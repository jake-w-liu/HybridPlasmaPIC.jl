function _ch_populations(species)
    species isa AbstractVector &&
        !isempty(species) &&
        axes(species, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("CUDA hybrid requires a nonempty ordered population vector"))
    all(ps -> ps isa _CHParticle, species) ||
        throw(ArgumentError("CUDA hybrid requires ordinary Float64 D1 DeviceMemory particles"))
    return _CHParticle[ps for ps in species]
end

function _CH.cuda_hybrid_stepper(n::Tuple, L::Tuple, model, species; kwargs...)
    throw(
        ArgumentError(
            "CUDA hybrid requires Float64 isothermal HybridModel and ordinary D1 DeviceMemory particles",
        ),
    )
end

function _CH.cuda_hybrid_stepper(
    n::Tuple,
    L::Tuple,
    model::_CHModel,
    species::AbstractVector;
    dt::Real,
    NB::Integer = 1,
    B = nothing,
)
    Sys.WORD_SIZE == 64 || throw(ArgumentError("CUDA hybrid requires a 64-bit Julia host"))
    length(n) == length(L) == 1 && n[1] isa Integer && L[1] isa Float64 ||
        throw(ArgumentError("CUDA hybrid requires one integer grid count and Float64 length"))
    2 <= n[1] <= _CH_INDEX_LIMIT ||
        throw(ArgumentError("CUDA hybrid grid exceeds its D1 index domain"))
    _ch_positive("L", L[1])
    _ch_model(model)
    h, nb = _ch_fixed_step(dt, NB)
    _CH._boris_accuracy_mode() ||
        throw(ArgumentError("CUDA hybrid requires host nearest rounding and gradual underflow"))
    CUDA.functional() && CUDA.cuFFT.functional() ||
        throw(ArgumentError("CUDA hybrid requires a functional CUDA device and cuFFT"))
    CUDA.allowscalar(false)
    owners = _ch_populations(species)
    scalars = Tuple{Float64,Float64}[]
    factors = _CH._BorisSpeciesFactors{Float64}[]
    for ps in owners
        _ch_positive("q", ps.q)
        _ch_positive("m", ps.m)
        qm = _ch_positive("q/m", ps.q / ps.m)
        push!(scalars, (ps.q, ps.m))
        push!(factors, _CH._BorisSpeciesFactors(ps.q, ps.m, qm))
    end
    metadata = _CH.SpectralOperators._fourier_metadata((Int(n[1]),), L)
    _ch_positive("dx", metadata.dx[1])
    ctx, dev = CUDA.context(owners[1].x[1]), CUDA.device(owners[1].x[1])
    counts = Int[length(ps.weight) for ps in owners]
    for (ps, k) in zip(owners, counts)
        0 <= k <= _CH_INDEX_LIMIT ||
            throw(ArgumentError("CUDA hybrid particle count exceeds its index domain"))
        foreach(a -> _ch_leaf(a, k, ctx, dev), _ch_particle_leaves(ps))
    end
    _ch_independent(Any[a for ps in owners for a in _ch_particle_leaves(ps)])
    CUDA.context()
    CUDA.stream()
    return CUDA.context!(ctx) do
        stream = CUDA.CuStream()
        _ch_scope(ctx, stream) do
            count = Int(n[1])
            ik = CUDA.CuArray{ComplexF64,1,CUDA.DeviceMemory}(metadata.ik[1])
            g = _CUDAHybridGrid(
                (count,),
                L,
                metadata.dx,
                metadata,
                ik,
                maximum(abs, metadata.kvec[1]),
            )
            fields, candidate = _ch_fields(count), _ch_fields(count)
            for f in (fields, candidate), a in _ch_field_arrays(f)
                fill!(a, 0.0)
            end
            motion = _CUDAHybridMotion[_ch_motion(k) for k in counts]
            binding = [_ch_binding(ps) for ps in owners]
            state = _ResidentCUDAHybrid(
                g,
                model,
                _CH.CIC(),
                h,
                nb,
                _ch_configuration(g, model, h, nb, scalars),
                owners,
                scalars,
                factors,
                counts,
                binding,
                fields,
                candidate,
                motion,
                _ch_work(count),
                _ch_fourier(count),
                CUDA.CuArray{Int64,1,CUDA.DeviceMemory}(undef, 3),
                ctx,
                dev,
                stream,
                false,
                false,
                0.0,
                0,
                0.0,
                Threads.Atomic{Bool}(false),
            )
            _ch_storage(state, owners)
            _ch_values(state, owners)
            if B !== nothing
                B isa NTuple{3,_CHF} || throw(
                    ArgumentError("CUDA hybrid B seed requires three ordinary device vectors"),
                )
                foreach(a -> _ch_leaf(a, count, ctx, dev), B)
                _ch_independent(Any[_ch_owned(state)..., B...])
                foreach(copyto!, fields.B, B)
            end
            _ch_launch!(state, "seed B", _ch_vector_values!, count, fields.B)
            minimum(fields.B[1]) == maximum(fields.B[1]) ||
                throw(ArgumentError("CUDA hybrid seed Bx must be constant"))
            CUDA.synchronize(stream)
            state
        end
    end
end

_CH.cuda_hybrid_stepper(n::Tuple, L::Tuple, model::_CHModel, ps::_CHParticle; kwargs...) =
    _CH.cuda_hybrid_stepper(n, L, model, _CHParticle[ps]; kwargs...)

function _ch_copy_fields!(destination, source)
    foreach(copyto!, _ch_field_arrays(destination), _ch_field_arrays(source))
    return nothing
end

function _ch_gather_fields!(st, w, x, E, B)
    n = length(x)
    _ch_launch!(st, "electric gather", _ch_gather!, n, w.Ep, E, x, st.g.L[1], st.g.dx[1], st.g.n[1])
    _ch_launch!(st, "magnetic gather", _ch_gather!, n, w.Bp, B, x, st.g.L[1], st.g.dx[1], st.g.n[1])
    return nothing
end

function _ch_scales(st, w, v, coefficient, mode)
    n = length(v[1])
    n == 0 && return (0.0, 0.0, 0.0, Inf)
    _ch_launch!(st, "motion scales", _ch_scales!, n, w.scales, v, w.Ep, w.Bp, coefficient, mode)
    return (maximum(w.scales[1]), maximum(w.scales[2]), maximum(w.scales[3]), minimum(w.scales[4]))
end

function _ch_euler_stage!(st, w, out, v, factor, mode)
    half = st.dt / 2
    hq = half * factor.qm
    scales = _ch_scales(st, w, v, hq, Val(:euler))
    _CH._euler_batch_certificate(half, factor, hq, scales...) ||
        throw(ArgumentError("CUDA hybrid Euler stage needs unavailable native completion"))
    _ch_launch!(st, "Euler half-kick", _ch_euler!, length(v[1]), out, v, w.Ep, w.Bp, hq, mode)
    return nothing
end

function _ch_moments!(st, n, u, positions, velocities, species)
    out = (n, u...)
    foreach(a -> fill!(a, 0.0), out)
    for (s, ps) in enumerate(species)
        foreach(a -> fill!(a, 0.0), st.work.raw)
        _ch_launch!(
            st,
            "population $s deposit",
            _ch_deposit!,
            length(ps.weight),
            st.work.raw,
            positions[s][1],
            velocities[s],
            ps.weight,
            st.g.L[1],
            st.g.dx[1],
            st.g.n[1],
        )
        _ch_launch!(
            st,
            "population $s moment finish",
            _ch_accumulate!,
            st.g.n[1],
            out,
            st.work.raw,
            ps.q,
            st.g.dx[1],
        )
    end
    _ch_launch!(st, "density/ion velocity", _ch_normalize!, st.g.n[1], n, u, st.model.nfloor)
    return nothing
end

function _ch_pressure_prepare!(st, pe, gradp, ninv, n)
    _ch_launch!(
        st,
        "pressure/reciprocal",
        _ch_pressure!,
        st.g.n[1],
        pe,
        ninv,
        n,
        st.model.closure.Te,
        st.model.nfloor,
    )
    _ch_gradient!(st, gradp[1], pe)
    return nothing
end

function _ch_constitutive!(st, E, J, B, u, gradp, ninv)
    _ch_curl!(st, J, B)
    _ch_launch!(
        st,
        "Ohm/source-risk",
        _ch_ohm!,
        st.g.n[1],
        E,
        u,
        B,
        J,
        gradp,
        ninv,
        st.model.η,
        st.model.ηH,
    )
    return nothing
end

function _ch_prepare_fields!(st, f)
    _ch_pressure_prepare!(st, f.pe, f.gradp, f.ninv, f.n)
    _ch_constitutive!(st, f.E, f.J, f.B, f.ui, f.gradp, f.ninv)
    return nothing
end

function _ch_rhs!(st, out, B)
    w = st.work
    _ch_constitutive!(st, w.E, w.J, B, w.frozen_ui, w.gradp, w.ninv)
    _ch_curl!(st, out, w.E)
    _ch_launch!(st, "Faraday sign", _ch_negate!, st.g.n[1], out)
    return nothing
end

function _ch_advance_B!(st, B)
    w = st.work
    hb = st.dt / st.NB
    half, whole, finish = _CH._RK4Affine(hb, 2), _CH._RK4Affine(hb, 1), _CH._RK4Affine(hb, 6)
    for _ = 1:st.NB
        _ch_rhs!(st, w.rates[1], B)
        _ch_launch!(st, "RK4 stage2", _ch_affine!, st.g.n[1], w.trial, B, (w.rates[1],), half)
        _ch_rhs!(st, w.rates[2], w.trial)
        _ch_launch!(st, "RK4 stage3", _ch_affine!, st.g.n[1], w.trial, B, (w.rates[2],), half)
        _ch_rhs!(st, w.rates[3], w.trial)
        _ch_launch!(st, "RK4 stage4", _ch_affine!, st.g.n[1], w.trial, B, (w.rates[3],), whole)
        _ch_rhs!(st, w.rates[4], w.trial)
        _ch_launch!(st, "RK4 finish", _ch_affine!, st.g.n[1], B, B, w.rates, finish)
    end
    return nothing
end

function _ch_recycled_motion(species, motion)
    return _CUDAHybridMotion[
        _CUDAHybridMotion(ps.x, ps.v, w.xmid, w.vpred, w.Ep, w.Bp, w.scales) for
        (ps, w) in zip(species, motion)
    ]
end

function _ch_replacements(species, motion)
    return _CHParticle[
        _CH.ParticleSet{1,Float64}(w.x, w.v, ps.weight, ps.id, ps.tag, ps.q, ps.m) for
        (ps, w) in zip(species, motion)
    ]
end

function _ch_commit!(st, particles, fields, counts, binding, recycled, time, count, last, primed)
    # Every argument/typed replacement is prepared, all device work synchronized,
    # and caller context/stream restored before this first public store.
    for s in eachindex(st.species)
        st.species[s].x = particles[s].x
        st.species[s].v = particles[s].v
        st.species[s].weight = particles[s].weight
        st.species[s].id = particles[s].id
        st.species[s].tag = particles[s].tag
    end
    old = st.fields
    st.fields = fields
    st.candidate = old
    st.motion = recycled
    st.counts = counts
    st.binding = binding
    st.time = time
    st.step = count
    st.last_dt = last
    st.primed = primed
    st.initialized = true
    return st
end

function _CH.init!(st::_ResidentCUDAHybrid, species::AbstractVector = st.species)
    _ch_acquire(st)
    try
        prepared = _ch_scope(st) do
            _ch_storage(st, species; counts = false, binding = false)
            _ch_values(st, species)
            _ch_launch!(st, "initial B", _ch_vector_values!, st.g.n[1], st.fields.B)
            minimum(st.fields.B[1]) == maximum(st.fields.B[1]) ||
                throw(ArgumentError("CUDA hybrid initial Bx must be constant"))
            counts = Int[length(ps.weight) for ps in species]
            motion = _CUDAHybridMotion[_ch_motion(n) for n in counts]
            _ch_independent(
                Any[_ch_owned(st)..., (a for w in motion for a in _ch_motion_arrays(w))...],
            )
            for (ps, w) in zip(species, motion)
                foreach(copyto!, w.x, ps.x)
                foreach(copyto!, w.v, ps.v)
            end
            foreach(copyto!, st.candidate.B, st.fields.B)
            _ch_moments!(
                st,
                st.candidate.n,
                st.candidate.ui,
                [w.x for w in motion],
                [w.v for w in motion],
                species,
            )
            _ch_prepare_fields!(st, st.candidate)
            _ch_fields_values(st, st.candidate)
            binding =
                [(w.x[1], w.v..., ps.weight, ps.id, ps.tag) for (ps, w) in zip(species, motion)]
            replacements = _ch_replacements(species, motion)
            recycled = _ch_recycled_motion(species, motion)
            CUDA.synchronize(st.stream)
            (replacements, st.candidate, counts, binding, recycled)
        end
        return _ch_commit!(st, prepared..., 0.0, 0, 0.0, false)
    finally
        _ch_release(st)
    end
end

function _CH.step!(st::_ResidentCUDAHybrid, species::AbstractVector = st.species; dt::Real = st.dt)
    _ch_acquire(st)
    try
        h = _ch_duration(dt, st.dt)
        prepared = _ch_scope(st) do
            _ch_boundary(st, species)
            _ch_zero(h) && return nothing
            time, count = _CH._validated_step_clock(st.time, st.step, h)
            time > st.time ||
                throw(ArgumentError("CUDA hybrid dt does not advance stored time"))
            _ch_copy_fields!(st.candidate, st.fields)
            for (s, ps) in enumerate(species)
                w = st.motion[s]
                foreach(copyto!, w.x, ps.x)
                foreach(copyto!, w.v, ps.v)
                _ch_gather_fields!(st, w, ps.x[1], st.fields.E, st.fields.B)
                st.primed || _ch_euler_stage!(st, w, w.v, w.v, st.factors[s], Val(:prime))
                hk = st.factors[s].qm * h / 2
                scales = _ch_scales(st, w, w.v, hk, Val(:boris))
                admitted, vbound = _CH._boris_batch_admission(hk, st.factors[s], h, scales...)
                xmax = isempty(ps.x[1]) ? 0.0 : maximum(abs, ps.x[1])
                admitted && _CH._coordinate_batch_certified(xmax, h, vbound) || throw(
                    ArgumentError("CUDA hybrid Boris/drift needs unavailable native completion"),
                )
                _ch_launch!(
                    st,
                    "Boris/drift/wrap",
                    _ch_boris_drift!,
                    length(ps.weight),
                    w.x[1],
                    w.v,
                    w.xmid[1],
                    ps.x[1],
                    w.v,
                    w.Ep,
                    w.Bp,
                    st.factors[s].qm,
                    h,
                    st.g.L[1],
                )
            end
            w = st.work
            _ch_moments!(
                st,
                w.frozen_n,
                w.frozen_ui,
                [m.xmid for m in st.motion],
                [m.v for m in st.motion],
                species,
            )
            _ch_pressure_prepare!(st, w.pe, w.gradp, w.ninv, w.frozen_n)
            _ch_advance_B!(st, st.candidate.B)
            positions = [m.x for m in st.motion]
            _ch_moments!(
                st,
                st.candidate.n,
                st.candidate.ui,
                positions,
                [m.v for m in st.motion],
                species,
            )
            _ch_prepare_fields!(st, st.candidate)
            for (s, m) in enumerate(st.motion)
                _ch_gather_fields!(st, m, m.x[1], st.candidate.E, st.candidate.B)
                _ch_euler_stage!(st, m, m.vpred, m.v, st.factors[s], Val(:predict))
            end
            _ch_moments!(
                st,
                st.candidate.n,
                st.candidate.ui,
                positions,
                [m.vpred for m in st.motion],
                species,
            )
            _ch_prepare_fields!(st, st.candidate)
            _ch_fields_values(st, st.candidate)
            for (ps, m) in zip(species, st.motion)
                _ch_launch!(
                    st,
                    "candidate particles",
                    _ch_particle_values!,
                    length(ps.weight),
                    m.x[1],
                    m.v,
                    ps.weight,
                    st.g.L[1],
                )
            end
            _ch_phase(true, true, time, count, h, h)
            binding = [
                (m.x[1], m.v..., ps.weight, ps.id, ps.tag) for (ps, m) in zip(species, st.motion)
            ]
            replacements = _ch_replacements(species, st.motion)
            recycled = _ch_recycled_motion(species, st.motion)
            CUDA.synchronize(st.stream)
            (replacements, st.candidate, copy(st.counts), binding, recycled, time, count, h, true)
        end
        prepared === nothing && return st
        return _ch_commit!(st, prepared...)
    finally
        _ch_release(st)
    end
end
