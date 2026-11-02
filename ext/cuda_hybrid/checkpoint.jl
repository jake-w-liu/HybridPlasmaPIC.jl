import Serialization

const _CH_SNAPSHOT_KEYS = (
    :integrator,
    :schema,
    :policy,
    :configuration,
    :metadata,
    :initialized,
    :primed,
    :time,
    :step,
    :last_dt,
    :velocity_phase,
    :velocity_offset,
    :species,
    :fields,
    :floor_count,
    :runtime,
)
const _CH_PARTICLE_KEYS = (:x, :v, :weight, :id, :tag, :q, :m)
const _CH_FIELD_KEYS = (:n, :ui, :B, :E, :J, :pe, :gradp, :ninv)
const _CH_RUNTIME_KEYS =
    (:julia, :cuda, :cufft, :spectraloperators, :device_name, :device_capability)

function _ch_keys(value, keys)
    value isa NamedTuple && Base.keys(value) == keys ||
        throw(ArgumentError("CUDA hybrid payload has unsupported keys or record type"))
    return nothing
end

function _ch_host_vector(a, ::Type{T}, n) where {T}
    a isa Vector{T} && length(a) == n ||
        throw(ArgumentError("CUDA hybrid payload requires ordinary matching host Vectors"))
    return nothing
end

function _ch_host_supported(arrays)
    all(a -> all(_ch_supported, a), arrays) || throw(
        ArgumentError("CUDA hybrid payload has nonfinite or nonzero subnormal physical values"),
    )
    return nothing
end

function _ch_host_state(s, st)
    _ch_keys(s, _CH_SNAPSHOT_KEYS)
    s.integrator === :cuda_hybrid_d1 && s.schema === 1 && s.policy === _CH_POLICY ||
        throw(ArgumentError("unsupported CUDA hybrid checkpoint family/schema/policy"))
    typeof(s.configuration) === typeof(st.configuration) &&
        isequal(s.configuration, st.configuration) ||
        throw(ArgumentError("CUDA hybrid checkpoint fixed configuration must match"))
    metadata = _CH.SpectralOperators._fourier_metadata(st.g.n, st.g.L)
    typeof(s.metadata) === typeof(metadata) && isequal(s.metadata, metadata) ||
        throw(ArgumentError("CUDA hybrid checkpoint Fourier metadata must match exactly"))
    _ch_phase(s.initialized, s.primed, s.time, s.step, s.last_dt, st.dt)
    s.velocity_offset isa Float64 &&
        isequal(s.velocity_offset, s.step == 0 ? 0.0 : -st.dt / 2) &&
        s.velocity_phase === (s.step == 0 ? :loaded_integer : :integer_fields_half_velocity) ||
        throw(ArgumentError("CUDA hybrid checkpoint velocity phase is inconsistent"))
    s.floor_count === 0 || throw(ArgumentError("CUDA hybrid requires inactive density floor"))
    _ch_keys(s.runtime, _CH_RUNTIME_KEYS)
    all(v -> v isa String && !isempty(v), values(s.runtime)) ||
        throw(ArgumentError("CUDA hybrid checkpoint version/device provenance must be strings"))
    s.species isa Vector && length(s.species) == length(st.species) ||
        throw(ArgumentError("CUDA hybrid checkpoint ordered population count differs"))
    arrays = AbstractArray[]
    positive_weight = false
    for (index, r) in enumerate(s.species)
        _ch_keys(r, _CH_PARTICLE_KEYS)
        r.q isa Float64 && r.m isa Float64 && isequal((r.q, r.m), st.scalars[index]) ||
            throw(ArgumentError("CUDA hybrid checkpoint original charges/masses differ"))
        r.weight isa Vector{Float64} ||
            throw(ArgumentError("CUDA hybrid checkpoint weight must be an ordinary Float64 Vector"))
        n = length(r.weight)
        0 <= n <= _CH_INDEX_LIMIT ||
            throw(ArgumentError("CUDA hybrid checkpoint particle count exceeds index domain"))
        r.x isa Tuple{Vector{Float64}} && r.v isa NTuple{3,Vector{Float64}} ||
            throw(ArgumentError("CUDA hybrid checkpoint particle dimension/storage is unsupported"))
        for a in (r.x..., r.v..., r.weight)
            _ch_host_vector(a, Float64, n)
        end
        _ch_host_vector(r.id, UInt64, n)
        _ch_host_vector(r.tag, UInt32, n)
        physical = (r.x..., r.v..., r.weight)
        _ch_host_supported(physical)
        all(x -> 0 <= x < st.g.L[1], r.x[1]) && all(w -> w >= 0, r.weight) || throw(
            ArgumentError(
                "CUDA hybrid checkpoint positions/weights are outside their admitted domain",
            ),
        )
        positive_weight |= any(w -> w > 0, r.weight)
        append!(arrays, (physical..., r.id, r.tag))
    end
    positive_weight ||
        throw(ArgumentError("CUDA hybrid initialized checkpoint needs positive ion inventory"))
    _ch_keys(s.fields, _CH_FIELD_KEYS)
    n = st.g.n[1]
    for name in (:ui, :B, :E, :J)
        getproperty(s.fields, name) isa NTuple{3,Vector{Float64}} || throw(
            ArgumentError(
                "CUDA hybrid checkpoint vector fields require three ordinary host Vectors",
            ),
        )
    end
    s.fields.gradp isa Tuple{Vector{Float64}} ||
        throw(ArgumentError("CUDA hybrid checkpoint gradp must have one spatial component"))
    physical = (
        s.fields.n,
        s.fields.ui...,
        s.fields.B...,
        s.fields.E...,
        s.fields.J...,
        s.fields.pe,
        s.fields.gradp...,
        s.fields.ninv,
    )
    for a in physical
        _ch_host_vector(a, Float64, n)
    end
    _ch_host_supported(physical)
    all(x -> x >= st.model.nfloor, s.fields.n) &&
        all(x -> x > 0, s.fields.ninv) &&
        all(x -> x >= 0, s.fields.pe) ||
        throw(ArgumentError("CUDA hybrid checkpoint density/pressure/reciprocal is invalid"))
    all(x -> x == s.fields.B[1][1], s.fields.B[1]) ||
        throw(ArgumentError("CUDA hybrid checkpoint Bx must be spatially constant"))
    append!(arrays, physical)
    append!(arrays, (s.metadata.midx..., s.metadata.kfull..., s.metadata.ik..., s.metadata.kvec...))
    destination_metadata =
        (st.g.metadata.midx..., st.g.metadata.kfull..., st.g.metadata.ik..., st.g.metadata.kvec...)
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _CH._particle_arrays_mightalias(arrays[i], arrays[j]) && throw(
                ArgumentError("CUDA hybrid payload physical/metadata arrays must be independent"),
            )
        end
        for a in destination_metadata
            _CH._particle_arrays_mightalias(arrays[i], a) && throw(
                ArgumentError("CUDA hybrid payload must be detached from destination metadata"),
            )
        end
    end
    return nothing
end

_CH._validate_checkpoint_state(s, st::_ResidentCUDAHybrid) = _ch_host_state(s, st)
_CH._checkpoint_state(st::_ResidentCUDAHybrid, species::AbstractVector) =
    _ch_snapshot_locked(st, species)

_ch_device_vector(a::Vector{T}) where {T} = CUDA.CuArray{T,1,CUDA.DeviceMemory}(a)

function _ch_restore_particles(records)
    return _CHParticle[
        _CH.ParticleSet{1,Float64}(
            map(_ch_device_vector, r.x),
            map(_ch_device_vector, r.v),
            _ch_device_vector(r.weight),
            _ch_device_vector(r.id),
            _ch_device_vector(r.tag),
            r.q,
            r.m,
        ) for r in records
    ]
end

function _ch_restore_fields(f)
    return _CUDAHybridFields(
        _ch_device_vector(f.n),
        map(_ch_device_vector, f.ui),
        map(_ch_device_vector, f.B),
        map(_ch_device_vector, f.E),
        map(_ch_device_vector, f.J),
        _ch_device_vector(f.pe),
        map(_ch_device_vector, f.gradp),
        _ch_device_vector(f.ninv),
    )
end

function _CH._restore_checkpoint_state!(st::_ResidentCUDAHybrid, species::AbstractVector, s)
    prepared = _ch_scope(st) do
        _ch_storage(st, species; counts = false, binding = false)
        _ch_host_state(s, st) # Original links must be rejected before detachment.
        saved = deepcopy(s)
        _ch_host_state(saved, st)
        _ch_device_metadata(st)
        particles = _ch_restore_particles(saved.species)
        fields = _ch_restore_fields(saved.fields)
        counts = Int[length(ps.weight) for ps in particles]
        motion = _CUDAHybridMotion[_ch_motion(n) for n in counts]
        arrays = Any[_ch_owned(st)..., _ch_field_arrays(fields)...]
        for ps in particles
            append!(arrays, _ch_particle_leaves(ps))
        end
        for w in motion
            append!(arrays, _ch_motion_arrays(w))
        end
        _ch_independent(arrays)
        _ch_values(st, particles)
        _ch_fields_values(st, fields)
        binding = [_ch_binding(ps) for ps in particles]
        CUDA.synchronize(st.stream)
        (
            particles,
            fields,
            counts,
            binding,
            motion,
            saved.time,
            saved.step,
            saved.last_dt,
            saved.primed,
        )
    end
    return _ch_commit!(st, prepared...)
end

"Save the complete detached CUDA physical state at an initialized ready boundary."
function _CH.save_checkpoint(
    path::AbstractString,
    st::_ResidentCUDAHybrid,
    species::AbstractVector = st.species,
)
    _ch_acquire(st)
    try
        return _CH._save_checkpoint(path, st, _ch_snapshot_locked(st, species))
    finally
        _ch_release(st)
    end
end

"Restore all saved CUDA physical fields, particles and phase without reconstructing them."
function _CH.load_checkpoint!(st::_ResidentCUDAHybrid, path::AbstractString)
    return _CH.load_checkpoint!(st, st.species, path)
end
function _CH.load_checkpoint!(
    st::_ResidentCUDAHybrid,
    species::AbstractVector,
    path::AbstractString,
)
    _ch_acquire(st)
    try
        return _CH._restore_checkpoint_state!(st, species, Serialization.deserialize(path))
    finally
        _ch_release(st)
    end
end

"Archive the same complete CUDA payload through the canonical checksummed run container."
function _CH.archive_run(
    path::AbstractString,
    st::_ResidentCUDAHybrid,
    species::AbstractVector = st.species;
    kwargs...,
)
    _ch_acquire(st)
    try
        return _CH._archive_run(path, st, _ch_snapshot_locked(st, species); kwargs...)
    finally
        _ch_release(st)
    end
end

"Verify the run checksum and continue through the same complete CUDA restore owner."
function _CH.load_archive!(st::_ResidentCUDAHybrid, path::AbstractString)
    return _CH.load_archive!(st, st.species, path)
end
function _CH.load_archive!(st::_ResidentCUDAHybrid, species::AbstractVector, path::AbstractString)
    _ch_acquire(st)
    try
        loaded = _CH.load_run(path)
        return _CH._restore_checkpoint_state!(st, species, loaded.state)
    finally
        _ch_release(st)
    end
end
