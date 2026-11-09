# Scalar bit admission is also used by the device kernels. No FP zero comparison
# may erase a forbidden stored subnormal before its range classification.
@inline _ch_bits(x::Float64) = reinterpret(UInt64, x) & 0x7fffffffffffffff
@inline _ch_zero(x::Float64) = _ch_bits(x) == 0
@inline _ch_supported(x::Float64) =
    _ch_bits(x) == 0 || 0x0010000000000000 <= _ch_bits(x) < 0x7ff0000000000000

function _ch_positive(name, x::Float64)
    _ch_supported(x) && x > 0 ||
        throw(ArgumentError("CUDA hybrid $name must be positive normal Float64"))
    return x
end

function _ch_model(model::_CHModel)
    _ch_supported(model.closure.Te) && model.closure.Te >= 0 ||
        throw(ArgumentError("CUDA hybrid Te must be zero or positive normal Float64"))
    _ch_positive("nfloor", model.nfloor)
    all(_CH._particle_stored_zero, (model.η, model.ηH, model.de2)) ||
        throw(ArgumentError("CUDA hybrid requires original η = ηH = de2 = 0"))
    return nothing
end

function _ch_duration(dt::Real, fixed::Union{Nothing,Float64} = nothing)
    isfinite(dt) && _CH._stored_nonnegative(dt) ||
        throw(ArgumentError("CUDA hybrid dt must be finite and nonnegative"))
    value = Float64(dt)
    isfinite(value) || throw(ArgumentError("CUDA hybrid dt is not finite in Float64"))
    original_zero = _CH._particle_stored_zero(dt)
    original_zero ||
        !_ch_zero(value) ||
        throw(ArgumentError("CUDA hybrid nonzero dt converts to stored zero"))
    original_zero && return value
    _ch_positive("dt", value)
    fixed === nothing ||
        value == fixed ||
        throw(ArgumentError("CUDA hybrid requires its fixed saved dt"))
    return value
end

function _ch_fixed_step(dt::Real, NB::Integer)
    h = _ch_duration(dt)
    _ch_positive("dt", h)
    1 <= NB <= typemax(Int) || throw(ArgumentError("CUDA hybrid NB must be a positive Int"))
    count = Int(NB)
    for (name, value) in (
        ("dt/2", h / 2),
        ("dt/NB", h / count),
        ("dt/NB/2", (h / count) / 2),
        ("dt/NB/6", (h / count) / 6),
    )
        _ch_positive(name, value)
    end
    return h, count
end

function _ch_phase(initialized, primed, time, count, last, h; ready = true)
    initialized isa Bool &&
        primed isa Bool &&
        time isa Float64 &&
        count isa Int &&
        last isa Float64 || throw(ArgumentError("CUDA hybrid phase has incorrect scalar types"))
    ready && !initialized && throw(ArgumentError("CUDA hybrid is not initialized"))
    isfinite(time) && time >= 0 && count >= 0 ||
        throw(ArgumentError("CUDA hybrid time/count must be finite and nonnegative"))
    if count == 0
        !primed && _ch_zero(time) && _ch_zero(last) ||
            throw(ArgumentError("CUDA hybrid count zero requires unprimed time0/last0"))
    else
        initialized && primed && last == h && h <= time ||
            throw(ArgumentError("CUDA hybrid advanced phase requires primed last_dt=h≤time"))
    end
    return nothing
end

function _ch_acquire(st::_ResidentCUDAHybrid)
    Threads.atomic_cas!(st.busy, false, true) == false ||
        throw(ArgumentError("CUDA hybrid already has an active operation"))
    return nothing
end
_ch_release(st::_ResidentCUDAHybrid) = (Threads.atomic_xchg!(st.busy, false); nothing)

function _ch_scope(f, ctx, stream)
    # Materialize an actual caller scope, including on a previously unused task.
    CUDA.context()
    CUDA.stream()
    return CUDA.context!(ctx) do
        CUDA.stream!(stream) do
            f()
        end
    end
end
_ch_scope(f, st::_ResidentCUDAHybrid) = _ch_scope(f, st.context, st.stream)

function _ch_leaf(a, n, ctx, dev)
    a isa Union{_CHF,_CHI,_CHT,_CHC,_CHStatus} ||
        throw(ArgumentError("CUDA hybrid requires ordinary contiguous DeviceMemory CuVectors"))
    length(a) == n && axes(a, 1) == Base.OneTo(n) ||
        throw(DimensionMismatch("CUDA hybrid device leaf has a mismatched count"))
    CUDA.context(a) == ctx && CUDA.device(a) == dev ||
        throw(ArgumentError("CUDA hybrid leaves must share the bound context/device"))
    # Accessing the manager also rejects a freed reference before publication.
    a.data[]
    return nothing
end

function _ch_alias(a, b)
    return a === b || a.data[] === b.data[] || Base.mightalias(a, b)
end

function _ch_independent(arrays)
    for i in eachindex(arrays), j = 1:(i-1)
        _ch_alias(arrays[i], arrays[j]) &&
            throw(ArgumentError("CUDA hybrid public/private/particle storage roles overlap"))
    end
    return nothing
end

function _ch_owned(st::_ResidentCUDAHybrid)
    arrays = Any[
        _ch_field_arrays(st.fields)...,
        _ch_field_arrays(st.candidate)...,
        _ch_work_arrays(st.work)...,
        st.fourier.input,
        st.fourier.accumulator,
        st.g.ik,
        st.status,
    ]
    for w in st.motion
        append!(arrays, _ch_motion_arrays(w))
    end
    for ps in st.species
        append!(arrays, _ch_particle_leaves(ps))
    end
    return arrays
end

function _ch_storage(st::_ResidentCUDAHybrid, species; counts = true, binding = true)
    species isa AbstractVector &&
        axes(species, 1) == Base.OneTo(length(species)) &&
        length(species) == length(st.species) ||
        throw(ArgumentError("CUDA hybrid requires its bound ordered populations"))
    _ch_model(st.model)
    isequal(_ch_configuration(st.g, st.model, st.dt, st.NB, st.scalars), st.configuration) ||
        throw(ArgumentError("CUDA hybrid fixed configuration changed"))
    all(
        length(a) == length(st.species) for
        a in (st.scalars, st.factors, st.counts, st.binding, st.motion)
    ) || throw(ArgumentError("CUDA hybrid population workspace/binding layout changed"))
    _ch_fixed_step(st.dt, st.NB)
    st.shape isa _CH.CIC || throw(ArgumentError("CUDA hybrid requires CIC"))
    _CH._boris_accuracy_mode() ||
        throw(ArgumentError("CUDA hybrid requires host nearest rounding and gradual underflow"))
    for (s, ps) in enumerate(species)
        ps === st.species[s] && ps isa _CHParticle ||
            throw(ArgumentError("CUDA hybrid population object/order or storage type changed"))
        isequal((ps.q, ps.m), st.scalars[s]) ||
            throw(ArgumentError("CUDA hybrid original q/m must match construction"))
        n = length(ps.weight)
        0 <= n <= _CH_INDEX_LIMIT ||
            throw(ArgumentError("CUDA hybrid particle count exceeds its index domain"))
        counts &&
            n != st.counts[s] &&
            throw(ArgumentError("CUDA hybrid count changed; call init! explicitly"))
        leaves = _ch_binding(ps)
        binding &&
            !all(i -> leaves[i] === st.binding[s][i], eachindex(leaves)) &&
            throw(ArgumentError("CUDA hybrid particle leaves changed; call init! explicitly"))
        foreach(a -> _ch_leaf(a, counts ? n : length(a), st.context, st.device), leaves)
        factor = st.factors[s]
        isequal((factor.q, factor.m, factor.qm), (ps.q, ps.m, ps.q / ps.m)) ||
            throw(ArgumentError("CUDA hybrid original-factor cache changed"))
    end
    n = st.g.n[1]
    for f in (st.fields, st.candidate), a in _ch_field_arrays(f)
        _ch_leaf(a, n, st.context, st.device)
    end
    for a in _ch_work_arrays(st.work)
        _ch_leaf(a, n, st.context, st.device)
    end
    for (s, w) in enumerate(st.motion), a in _ch_motion_arrays(w)
        _ch_leaf(a, counts ? st.counts[s] : length(a), st.context, st.device)
    end
    for a in (st.fourier.input, st.fourier.accumulator, st.g.ik)
        _ch_leaf(a, n, st.context, st.device)
    end
    for plan in (st.fourier.forward, st.fourier.inverse)
        plan.ctx == st.context &&
            plan.input_size == st.g.n &&
            plan.output_size == st.g.n &&
            plan.region == (1,) &&
            plan.handle != C_NULL ||
            throw(ArgumentError("CUDA hybrid FFT plan context/shape/handle changed"))
    end
    _ch_leaf(st.status, 3, st.context, st.device)
    _ch_independent(_ch_owned(st))
    return nothing
end

function _ch_status!(status, stream, label)
    CUDA.synchronize(stream)
    record = Array(status) # Exactly three Int64 control values, never field/particle data.
    record[1] == 0 || throw(
        ArgumentError(
            "CUDA hybrid $label rejected category $(record[1]) at index $(record[2]), component $(record[3])",
        ),
    )
    return nothing
end

function _ch_launch!(st, label, kernel::F, n::Int, args...) where {F}
    fill!(st.status, Int64(0))
    if n > 0
        CUDA.@cuda fastmath = false threads = 256 blocks = cld(n, 256) kernel(st.status, args..., n)
    end
    _ch_status!(st.status, st.stream, label)
    return nothing
end

function _ch_values(st::_ResidentCUDAHybrid, species = st.species)
    for (s, ps) in enumerate(species)
        # init! may admit a changed count; validate every newly loaded leaf before
        # the unchecked device reads use the weight count as their loop bound.
        foreach(a -> _ch_leaf(a, length(ps.weight), st.context, st.device), _ch_particle_leaves(ps))
        _ch_launch!(
            st,
            "population $s",
            _ch_particle_values!,
            length(ps.weight),
            ps.x[1],
            ps.v,
            ps.weight,
            st.g.L[1],
        )
    end
    return nothing
end

function _ch_fields_values(st, f; density = true)
    _ch_launch!(
        st,
        "physical fields",
        _ch_field_values!,
        st.g.n[1],
        _ch_field_arrays(f),
        f.n,
        f.pe,
        f.ninv,
        st.model.nfloor,
        density,
    )
    minimum(f.B[1]) == maximum(f.B[1]) ||
        throw(ArgumentError("CUDA hybrid Bx must be spatially constant"))
    return nothing
end

function _ch_boundary(st, species = st.species)
    _ch_storage(st, species)
    _ch_phase(st.initialized, st.primed, st.time, st.step, st.last_dt, st.dt)
    _ch_values(st, species)
    _ch_fields_values(st, st.fields)
    return nothing
end
