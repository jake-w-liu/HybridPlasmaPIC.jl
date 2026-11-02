function _ch_runtime(st)
    modules = (CUDA, CUDA.cuFFT, _CH.SpectralOperators)
    versions = map(Base.pkgversion, modules)
    all(v -> v isa VersionNumber, versions) ||
        throw(ArgumentError("CUDA hybrid dependency version provenance is unavailable"))
    return (
        julia = string(VERSION),
        cuda = string(versions[1]),
        cufft = string(versions[2]),
        spectraloperators = string(versions[3]),
        device_name = CUDA.name(st.device),
        device_capability = string(CUDA.capability(st.device)),
    )
end

function _ch_host_particle(ps)
    return (
        x = map(Array, ps.x),
        v = map(Array, ps.v),
        weight = Array(ps.weight),
        id = Array(ps.id),
        tag = Array(ps.tag),
        q = ps.q,
        m = ps.m,
    )
end

function _ch_host_fields(f)
    return (
        n = Array(f.n),
        ui = map(Array, f.ui),
        B = map(Array, f.B),
        E = map(Array, f.E),
        J = map(Array, f.J),
        pe = Array(f.pe),
        gradp = map(Array, f.gradp),
        ninv = Array(f.ninv),
    )
end

function _ch_device_metadata(st)
    expected = _CH.SpectralOperators._fourier_metadata(st.g.n, st.g.L)
    isequal(st.g.metadata, expected) &&
        isequal(Array(st.g.ik), expected.ik[1]) &&
        st.g.gain == maximum(abs, expected.kvec[1]) ||
        throw(ArgumentError("CUDA hybrid spectral metadata changed"))
    return deepcopy(expected)
end

function _ch_snapshot_locked(st, species)
    return _ch_scope(st) do
        _ch_boundary(st, species)
        CUDA.synchronize(st.stream)
        metadata = _ch_device_metadata(st) # Explicit snapshot boundary permits full transfer.
        saved = (
            integrator = :cuda_hybrid_d1,
            schema = 1,
            policy = _CH_POLICY,
            configuration = st.configuration,
            metadata = metadata,
            initialized = st.initialized,
            primed = st.primed,
            time = st.time,
            step = st.step,
            last_dt = st.last_dt,
            velocity_phase = st.step == 0 ? :loaded_integer : :integer_fields_half_velocity,
            velocity_offset = st.step == 0 ? 0.0 : -st.dt / 2,
            species = [_ch_host_particle(ps) for ps in species],
            fields = _ch_host_fields(st.fields),
            floor_count = 0,
            runtime = _ch_runtime(st),
        )
        _ch_host_state(saved, st)
        saved
    end
end

"""
    cuda_hybrid_host_snapshot(st, species=st.species)

Synchronize an initialized completed CUDA state and return detached ordinary
host particles, all sixteen physical fields, fixed metadata and actual phase.
This explicit full transfer does not prime, recenter or recompute saved fields.
Refetch live device leaves after mutating operations; snapshots remain detached.
"""
function _CH.cuda_hybrid_host_snapshot(
    st::_ResidentCUDAHybrid,
    species::AbstractVector = st.species,
)
    _ch_acquire(st)
    try
        return _ch_snapshot_locked(st, species)
    finally
        _ch_release(st)
    end
end

"""
    cuda_hybrid_host_diagnostics(st, species=st.species)

Take one explicit host snapshot and use the canonical host energy, momentum and
Fourier divergence owners. Energy uses carried half-time velocity after the
first step, alongside integer fields/predictor density; phase is reported.
"""
function _CH.cuda_hybrid_host_diagnostics(
    st::_ResidentCUDAHybrid,
    species::AbstractVector = st.species,
)
    s = _CH.cuda_hybrid_host_snapshot(st, species)
    populations = _CH.ParticleSet{1,Float64}[
        _CH.ParticleSet{1,Float64}(r.x, r.v, r.weight, r.id, r.tag, r.q, r.m) for r in s.species
    ]
    g = _CH.FourierGrid(s.configuration.n, s.configuration.L)
    divergence = zeros(Float64, g.n)
    _CH.divergence!(divergence, s.fields.B, g)
    return (
        energy = _CH.energy_budget(
            populations,
            s.fields.B,
            s.fields.n,
            _CH.IsothermalElectrons(s.configuration.Te),
            g,
        ),
        momentum = _CH.momentum_budget(populations, s.fields.B, g),
        divergence_l2 = LinearAlgebra.norm(divergence),
        divergence_max = maximum(abs, divergence),
        time = s.time,
        step = s.step,
        velocity_phase = s.velocity_phase,
        velocity_offset = s.velocity_offset,
        floor_count = s.floor_count,
    )
end
