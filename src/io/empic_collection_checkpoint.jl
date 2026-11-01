# Collection payloads extend the existing checkpoint/archive owners. Workspaces,
# gathers and FFT scratch remain recomputable and are not serialized.
const _EMPIC_COLLECTION_CHECKPOINT_FIELDS = (
    :integrator, :D, :T, :ncell, :L, :n0, :c, :shape, :n_sub,
    :particle_layout, :schema_version, :relativistic, :mobile,
    :neutrality, :charge_rtol, :initialized, :species, :fields,
    :time, :step, :last_dt, :last_sub_dt,
)

function _validate_empic_collection_checkpoint(s, st::EMPIC{D,T}) where {D,T}
    _validate_checkpoint_container(s, _EMPIC_COLLECTION_CHECKPOINT_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    s.integrator === :empic_collection && s.particle_layout === :collection &&
        s.schema_version === 1 ||
        throw(ArgumentError("unsupported electromagnetic collection checkpoint schema"))
    work = _empic_collection(st)
    s.n0 isa T && isequal(s.n0, st.n0) ||
        throw(ArgumentError("checkpoint background_charge does not match destination"))
    s.neutrality === work.neutrality ||
        throw(ArgumentError("checkpoint neutrality policy does not match destination"))
    s.charge_rtol isa T && isequal(s.charge_rtol, work.charge_rtol) ||
        throw(ArgumentError("checkpoint charge_rtol does not match destination"))
    s.mobile isa AbstractVector{Bool} &&
        axes(s.mobile, 1) == Base.OneTo(length(work.mobile)) && s.mobile == work.mobile ||
        throw(ArgumentError("checkpoint per-population mobility does not match destination"))
    s.relativistic isa Bool ||
        throw(ArgumentError("checkpoint relativistic flag must be Bool"))
    typeof(s.shape) === typeof(st.shape) && s.shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("checkpoint shape does not match destination"))
    s.c isa T || throw(ArgumentError("checkpoint c precision does not match destination"))
    _require_finite_positive_real("checkpoint c", s.c, T)
    s.n_sub isa Int && s.n_sub >= 1 ||
        throw(ArgumentError("checkpoint n_sub must be a positive Int"))
    s.initialized isa Bool || throw(ArgumentError("checkpoint initialized flag must be Bool"))
    _validate_checkpoint_time_state(s, T)
    s.last_sub_dt isa T && isfinite(s.last_sub_dt) &&
        zero(T) <= s.last_sub_dt <= s.last_dt &&
        (s.step == 0 ? iszero(s.last_sub_dt) : s.last_sub_dt > zero(T)) ||
        throw(ArgumentError("checkpoint fast timestep history is inconsistent"))
    !s.initialized && s.step != 0 &&
        throw(ArgumentError("advanced electromagnetic checkpoint must be initialized"))
    s.species isa AbstractVector && axes(s.species, 1) == Base.OneTo(length(s.species)) ||
        throw(ArgumentError("checkpoint species must be an ordered one-based vector"))
    populations = ParticleSet{D,T}[_checkpoint_species(record, Val(D), T) for record in s.species]
    _empic_collection_validate_populations(st, populations, s.c, s.relativistic)
    _validate_checkpoint_fields(s.fields, _empic_checkpoint_fields(st))
    arrays = (s.fields.E..., s.fields.B..., s.fields.J..., s.fields.rho_n, s.fields.rho_np1)
    for (i, a) in enumerate(arrays)
        all(isfinite, a) || throw(ArgumentError("checkpoint electromagnetic fields must be finite"))
        for j = 1:(i-1)
            Base.mightalias(a, arrays[j]) &&
                throw(ArgumentError("checkpoint electromagnetic fields must be independent"))
        end
        for ps in populations
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint particles must not alias fields"))
        end
    end
    return nothing
end

function _checkpoint_state(st::EMPIC{D,T}, species::AbstractVector{<:ParticleSet}) where {D,T}
    work = _empic_collection(st)
    _empic_collection_validate_species(st, species)
    _empic_collection_validate_phase(st)
    work.initialized[] && _empic_collection_bound(work, species)
    return (
        integrator = :empic_collection, D = D, T = T, ncell = st.g.n, L = st.g.L,
        n0 = st.n0, c = st.c, shape = st.shape, n_sub = st.n_sub,
        particle_layout = :collection, schema_version = 1, relativistic = st.relativistic,
        mobile = copy(work.mobile), neutrality = work.neutrality, charge_rtol = work.charge_rtol,
        initialized = work.initialized[], species = [_checkpoint_particle_state(ps) for ps in species],
        fields = _empic_checkpoint_fields(st), time = st.time[], step = st.step[],
        last_dt = st.last_dt[], last_sub_dt = st.last_sub_dt[],
    )
end

function _restore_checkpoint_state!(
    st::EMPIC{D,T}, species::AbstractVector{<:ParticleSet}, s,
) where {D,T}
    work = _empic_collection(st)
    _validate_checkpoint_state(s, st)
    _empic_collection_validate_storage(st, species)
    length(species) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    all(ps -> ps isa eltype(work.owners), species) ||
        throw(ArgumentError("population storage type changed; construct a new EMPIC state"))
    saved_populations = [_checkpoint_species(record, Val(D), T) for record in s.species]
    saved_fields = (s.fields.E..., s.fields.B..., s.fields.J..., s.fields.rho_n, s.fields.rho_np1)
    for (ps, record) in zip(species, s.species)
        _validate_checkpoint_destination(ps, length(record.weight))
        for saved in saved_populations
            _particle_sets_mightalias(ps, saved) &&
                throw(ArgumentError("checkpoint source and destination particles must be independent"))
        end
        for a in saved_fields
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint fields must not alias destination particles"))
        end
    end
    # As in the collisional checkpoint owner, detach every incoming leaf before
    # publication. In-memory payload fields can be a permutation of live fields;
    # copying one directly would overwrite a later source leaf.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    # Reserve ordinary Vector capacity before any values or lengths publish.
    # Fixed strided views have already passed the exact-count preflight above.
    for (ps, record) in zip(species, saved.species)
        count = length(record.weight)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            if length(a) != count
                a isa Vector || throw(ArgumentError("only Vector checkpoint destinations may resize"))
                sizehint!(a, count)
            end
        end
    end
    sizehint!(work.owners, length(species))
    for (ps, record) in zip(species, saved.species)
        _restore_checkpoint_particles!(ps, record)
    end
    _restore_checkpoint_fields!(_empic_checkpoint_fields(st), saved.fields)
    st.c, st.shape, st.n_sub = saved.c, saved.shape, saved.n_sub
    st.relativistic = saved.relativistic
    st.time[], st.step[] = saved.time, saved.step
    st.last_dt[], st.last_sub_dt[] = saved.last_dt, saved.last_sub_dt
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = saved.initialized
    return st
end
