# Ordered synchronous populations share the ordinary EMPIC Maxwell owner. The
# private ordinary state stages field publication; it does not define a solver.

"""
    EMPIC(g, species::AbstractVector; mobile=trues(length(species)),
          background_charge=0, neutrality=:require, charge_rtol=64eps(T),
          c=5, shape=CIC(), n_sub=1, relativistic=false)

Periodic CPU Float32/Float64 electromagnetic PIC with arbitrary finite species
charges and positive masses. `mobile` selects populations that kick and drift at
every Maxwell substep `dt/n_sub`. Fixed populations require zero stored velocity,
deposit charge only and remain unchanged. NGP, CIC and TSC are supported.
Relativistic collections and independent species timesteps are not implemented.

`background_charge` is a fixed uniform charge density. `:require` uses the same
exact physical-charge tolerance as electrostatic collections; `:neutralize`
explicitly supplies a uniform compensator, reported by diagnostics. Neither
changes particle weights. Initialize with `init_empic!(state, species)` to bind
ordered population identities and physical full-time velocities. Replacing or
reordering populations requires initialization. Growth requires new velocities
at the existing common fast half level. Zero dt preserves that phase.

Particles, physical fields and clock publish after the complete operation
succeeds. Scratch and shared FFT work may change on rejection. Concurrent calls
require independent states, populations and grids. The integer-count constructor
retains its electron/ion contract and separate checkpoint layout.
"""
function EMPIC(
    g::FourierGrid{D,T}, species::AbstractVector;
    mobile::AbstractVector{Bool} = trues(length(species)),
    background_charge = zero(T), neutrality::Symbol = :require,
    charge_rtol = 64eps(T), c = 5.0, shape::ShapeFunction = CIC(),
    n_sub::Integer = 1, relativistic::Bool = false,
) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("EMPIC supports D = 1, 2, or 3"))
    T in (Float32, Float64) ||
        throw(ArgumentError("electromagnetic collections require Float32 or Float64"))
    relativistic && throw(ArgumentError("relativistic EMPIC collections are not implemented"))
    shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("electromagnetic collections require NGP, CIC or TSC"))
    neutrality in (:require, :neutralize) ||
        throw(ArgumentError("neutrality must be :require or :neutralize"))
    _espic_validate_collection_type(species, g)
    axes(mobile, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("mobile must have one Boolean per species with one-based axes"))
    n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    ns = Int(n_sub)
    background = _require_finite_real("background_charge", background_charge, T)
    rtol = _require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
    cT = _require_finite_positive_real("c", c, T)
    stage = _new_empic_state(g, 0, background, cT, shape, false, false, one(T), ns, nothing)
    P = eltype(species)
    work = _EMPICCollectionWork{D,T,P,typeof(stage)}(
        neutrality, rtol, BitVector(mobile), P[], _ElectrostaticHostParticle{D,T}[],
        NTuple{D,Vector{T}}[], stage, Ref(false),
    )
    es = _new_empic_state(g, 0, background, cT, shape, false, false, one(T), ns, work)
    _empic_collection_validate_species(es, species)
    return es
end

_require_scalar_empic(::EMPIC1D) = nothing
_require_scalar_empic(es::EMPIC) = es.collection === nothing ? nothing :
    throw(ArgumentError("collection EMPIC requires its ordered species vector"))

# Physical arrays carried across completed boundaries. Persistence, diagnostics
# and collection stepping share this inventory; gathers/derivatives are scratch.
_empic_physical_fields(es::EMPIC) =
    (E = es.E, B = es.B, J = es.J, rho_n = es.rho_n, rho_np1 = es.rho_np1)

function _empic_collection_require_finite_fields(es::EMPIC)
    for components in values(_empic_physical_fields(es))
        finite = components isa Tuple ? all(a -> all(isfinite, a), components) :
            all(isfinite, components)
        finite || throw(ArgumentError("electromagnetic physical fields must be finite"))
    end
    return nothing
end

function _empic_collection(es::EMPIC)
    work = es.collection
    work === nothing && throw(ArgumentError("species vectors require EMPIC(g, species; ...)"))
    return work
end

function _empic_collection_bound(work, species)
    length(work.owners) == length(species) &&
        all(a === b for (a, b) in zip(work.owners, species)) ||
        throw(ArgumentError("electromagnetic population membership/order changed; call init_empic!"))
    return nothing
end

function _empic_collection_validate_storage(es::EMPIC{D,T}, species) where {D,T}
    work = _empic_collection(es)
    _espic_validate_collection_type(species, es.g)
    length(work.mobile) == length(species) ||
        throw(ArgumentError("electromagnetic population count must match mobility flags"))
    # The staged and public field inventories are deliberately identical to the
    # existing transaction inventory, including continuity diagnostic snapshots.
    fields = (_empic_transaction_fields(es)..., _empic_transaction_fields(work.stage)...)
    work.stage.g === es.g || throw(ArgumentError("electromagnetic stage must share its owning grid"))
    fft = (es.g.cbuf, es.g.tbuf, es.g.abuf, es.g.sbuf)
    metadata = (es.g.midx..., es.g.kfull..., es.g.ik..., es.g.kvec...)
    for (i, a) in enumerate(fields)
        _require_grid_array(:electromagnetic_field, a, es.g)
        for j = 1:(i-1)
            Base.mightalias(a, fields[j]) &&
                throw(ArgumentError("electromagnetic public and staged fields must not alias"))
        end
        for b in (fft..., metadata...)
            Base.mightalias(a, b) &&
                throw(ArgumentError("electromagnetic fields must not alias Fourier storage"))
        end
    end
    buffers = (fields..., fft..., metadata...,
        _empic_electron_gathers(es)..., _empic_ion_gathers(es)...,
        _empic_electron_gathers(work.stage)..., _empic_ion_gathers(work.stage)...)
    gathers = (_empic_electron_gathers(es)..., _empic_ion_gathers(es)...,
        _empic_electron_gathers(work.stage)..., _empic_ion_gathers(work.stage)...)
    for (i, a) in enumerate(gathers)
        for b in (fields..., fft..., metadata..., gathers[1:(i-1)]...)
            Base.mightalias(a, b) &&
                throw(ArgumentError("electromagnetic gather storage must be independent"))
        end
    end
    clocks = (es.time, es.step, es.last_dt, es.last_sub_dt,
        work.stage.time, work.stage.step, work.stage.last_dt, work.stage.last_sub_dt)
    for (i, a) in enumerate(clocks), j = 1:(i-1)
        a === clocks[j] && throw(ArgumentError("electromagnetic clock storage must be independent"))
    end
    for (i, ps) in enumerate(species)
        _check_particle_storage_axes(ps)
        all(a -> a isa StridedVector, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)) ||
            throw(ArgumentError("electromagnetic collections require CPU strided particle storage"))
        _particle_components_mightalias(ps) &&
            throw(ArgumentError("electromagnetic particle components must not alias"))
        for j = 1:(i-1)
            _particle_sets_mightalias(ps, species[j]) &&
                throw(ArgumentError("electromagnetic populations must have independent storage"))
        end
        for a in buffers
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("electromagnetic particles must not alias solver storage"))
        end
        for staged in work.staged
            _particle_sets_mightalias(ps, staged) &&
                throw(ArgumentError("electromagnetic particles must not alias staged particles"))
        end
        for mid in work.midpoints, a in mid
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("electromagnetic particles must not alias midpoint storage"))
        end
    end
    return nothing
end

function _empic_collection_validate_species(es::EMPIC{D,T}, species) where {D,T}
    _empic_collection_validate_storage(es, species)
    work = _empic_collection(es)
    es.relativistic && throw(ArgumentError("relativistic EMPIC collections are not implemented"))
    es.mobile && throw(ArgumentError("collection EMPIC mobility is specified per population"))
    _require_finite_real("background_charge", es.n0, T)
    _require_finite_positive_real("c", es.c, T)
    es.n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    for (i, ps) in enumerate(species)
        _validated_boris_factor(ps)
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "electromagnetic collection")
        all(w -> isfinite(w) && w >= zero(T), ps.weight) ||
            throw(ArgumentError("electromagnetic weights must be finite and nonnegative"))
        !work.mobile[i] && !all(a -> all(iszero, a), ps.v) &&
            throw(ArgumentError("fixed electromagnetic populations require zero velocity"))
    end
    charges = _pic_collection_charges(es.g, es.n0, species)
    _espic_check_neutrality(work, charges)
    return charges
end

function _empic_collection_prepare!(es::EMPIC{D,T}, species) where {D,T}
    work = _empic_collection(es)
    _pic_prepare_species_stages!(work.staged, species)
    while length(work.midpoints) < length(species)
        push!(work.midpoints, ntuple(_ -> T[], Val(D)))
    end
    resize!(work.midpoints, length(species))
    for (ps, mid) in zip(species, work.midpoints), a in mid
        resize!(a, nparticles(ps))
    end
    stage = work.stage
    stage.n0, stage.c, stage.shape, stage.n_sub = es.n0, es.c, es.shape, es.n_sub
    foreach(copyto!, _empic_transaction_fields(stage), _empic_transaction_fields(es))
    stage.time[], stage.step[] = es.time[], es.step[]
    stage.last_dt[], stage.last_sub_dt[] = es.last_dt[], es.last_sub_dt[]
    return work
end

function _empic_collection_charge_target(es, work, charges)
    return work.neutrality === :neutralize ? zero(charges.net) :
        charges.net / prod(_exact_rational, es.g.dx)
end

function _empic_collection_charge!(rho, stage, work, target)
    _charge_density_sum!(rho, stage.ne, work.staged, stage.g, stage.shape;
        background = stage.n0)
    all(isfinite, rho) ||
        throw(ArgumentError("electromagnetic deposited charge exceeds numeric range"))
    shift = (_exact_float_sum(rho) - target) / length(rho)
    _empic_complete_charge_shift!(rho, shift)
    all(isfinite, rho) ||
        throw(ArgumentError("electromagnetic reconciled charge exceeds numeric range"))
    return rho
end

function _empic_collection_publish!(es, work; particles::Bool)
    if particles
        for (ps, staged) in zip(work.owners, work.staged)
            foreach(copyto!, ps.x, staged.x)
            foreach(copyto!, ps.v, staged.v)
        end
    end
    stage = work.stage
    foreach(copyto!, _empic_transaction_fields(es), _empic_transaction_fields(stage))
    es.time[], es.step[] = stage.time[], stage.step[]
    es.last_dt[], es.last_sub_dt[] = stage.last_dt[], stage.last_sub_dt[]
    return es
end

function init_empic!(es::EMPIC{D,T}, species::AbstractVector) where {D,T}
    work = _empic_collection(es)
    all(ps -> ps isa eltype(work.owners), species) ||
        throw(ArgumentError("population storage type changed; construct a new EMPIC state"))
    charges = _empic_collection_validate_species(es, species)
    all(a -> all(isfinite, a), (es.E..., es.B...)) ||
        throw(ArgumentError("electromagnetic fields must be finite"))
    _empic_collection_prepare!(es, species)
    sizehint!(work.owners, length(species))
    stage = work.stage
    _empic_collection_charge!(stage.rho_n, stage, work,
        _empic_collection_charge_target(es, work, charges))
    _enforce_gauss!(stage)
    all(a -> all(isfinite, a), stage.E) ||
        throw(ArgumentError("electromagnetic initialized field exceeds numeric range"))
    copyto!(stage.rho_np1, stage.rho_n)
    foreach(a -> fill!(a, zero(T)), (stage.J..., stage.Jraw...))
    foreach(copyto!, stage.Bc, stage.B)
    stage.time[], stage.step[] = zero(T), 0
    stage.last_dt[], stage.last_sub_dt[] = zero(T), zero(T)
    _empic_collection_publish!(es, work; particles = false)
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = true
    return es
end

function _empic_collection_current!(stage::EMPIC{D,T}, work) where {D,T}
    foreach(a -> fill!(a, zero(T)), stage.Jraw)
    eligible, lost = true, false
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_empic_electron_buffers!(stage, nparticles(ps))
        mid = work.midpoints[i]
        _swap_positions_with_midpoints!(ps, mid)
        states = try
            _add_species_current!(stage.Jraw, stage, ps, stage.worke, stage.ne, true)
        finally
            _swap_positions_with_midpoints!(ps, mid)
        end
        eligible &= all(s -> s.eligible, states)
        lost |= any(s -> s.lost, states)
        lost |= any(a -> any(_deposition_range_lost, a), stage.Jraw)
    end
    if eligible && lost
        for c = 1:3
            exact = fill(zero(Rational{BigInt}), stage.g.n)
            for (i, ps) in enumerate(work.staged)
                work.mobile[i] || continue
                contribution = _exact_midpoint_current(ps, work.midpoints[i], stage.g, stage.shape, c)
                for I in eachindex(exact)
                    exact[I] += contribution[I]
                end
            end
            for I in eachindex(exact)
                stage.Jraw[c][I] = _from_exact(T, exact[I])
            end
        end
    end
    all(a -> all(isfinite, a), stage.Jraw) ||
        throw(ArgumentError("electromagnetic deposited current exceeds numeric range"))
    return stage.Jraw
end

function _empic_collection_prime!(stage::EMPIC{D,T}, work, sub_dt::T) where {D,T}
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_empic_electron_buffers!(stage, nparticles(ps))
        gather_vector!(stage.Ep, stage.E, ps, stage.g, stage.shape)
        gather_vector!(stage.Bp, stage.B, ps, stage.g, stage.shape)
        _prime_leapfrog!(ps.v, stage.Ep, stage.Bp, _validated_boris_factor(ps),
            sub_dt / 2, nparticles(ps), false, stage.c)
    end
    _empic_prime_magnetic!(stage, stage.g, sub_dt)
    return stage
end

function _empic_collection_substep!(stage::EMPIC{D,T}, work, dt::T, previous::T, target) where {D,T}
    _empic_collection_charge!(stage.rho_n, stage, work, target)
    parts = _empic_center_magnetic!(stage, stage.g, previous)
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_empic_electron_buffers!(stage, nparticles(ps))
        gather_vector!(stage.Ep, stage.E, ps, stage.g, stage.shape)
        gather_vector!(stage.Bp, stage.Bc, ps, stage.g, stage.shape)
        _push_velocities!(stage, ps, stage.Ep, stage.Bp, previous / 2 + dt / 2)
        _drift_record_mid!(ps, work.midpoints[i], dt, stage.g.L)
        apply_periodic!(ps, ntuple(_ -> zero(T), Val(D)), stage.g.L)
        _require_finite_particle_velocities(ps, "electromagnetic collection output")
    end
    _empic_collection_current!(stage, work)
    _empic_collection_charge!(stage.rho_np1, stage, work, target)
    _correct_current_continuity!(stage, dt)
    _empic_update_fields!(stage, stage.g, stage.c * stage.c, parts, previous, dt)
    _empic_collection_require_finite_fields(stage)
    all(a -> all(isfinite, a), stage.Bc) ||
        throw(ArgumentError("electromagnetic centered magnetic field exceeds numeric range"))
    return stage
end

function _empic_collection_validate_phase(es::EMPIC{D,T}) where {D,T}
    isfinite(es.time[]) && es.time[] >= zero(T) && es.step[] >= 0 ||
        throw(ArgumentError("electromagnetic clock must be finite and nonnegative"))
    previous = es.last_sub_dt[]
    full = es.last_dt[]
    isfinite(full) && isfinite(previous) && zero(T) <= previous <= full &&
        (es.step[] == 0 ? iszero(full) && iszero(previous) : previous > zero(T)) ||
        throw(ArgumentError("electromagnetic timestep history is inconsistent with its phase"))
    return nothing
end

function step_empic!(es::EMPIC{D,T}, species::AbstractVector, dt::Real) where {D,T}
    dtT = _require_resolved_dt("step_empic!",
        _validated_nonnegative_dt(T, dt; name = "step_empic!"))
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    work.initialized[] || throw(ArgumentError("initialize the electromagnetic collection before stepping"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    iszero(dtT) && return es
    sub_dt = _require_resolved_dt("step_empic! substep dt/n_sub", dtT / es.n_sub)
    sub_dt > zero(T) || throw(ArgumentError("dt/n_sub must be representable and positive"))
    next_time, next_step = _validated_step_clock(es.time[], es.step[], sub_dt, es.n_sub)
    _empic_collection_prepare!(es, species)
    stage = work.stage
    stage.step[] == 0 && _empic_collection_prime!(stage, work, sub_dt)
    previous = stage.step[] == 0 ? sub_dt : stage.last_sub_dt[]
    # q and weights are unchanged by canonical kicks/drifts. Reuse their exact
    # conserved total rather than re-summing all particles at both endpoints.
    target = _empic_collection_charge_target(es, work, charges)
    for _ = 1:es.n_sub
        _empic_collection_substep!(stage, work, sub_dt, previous, target)
        previous = sub_dt
    end
    stage.time[], stage.step[] = next_time, next_step
    stage.last_dt[], stage.last_sub_dt[] = dtT, sub_dt
    return _empic_collection_publish!(es, work; particles = true)
end

"""
    electromagnetic_diagnostics(state::EMPIC, species)

Exact represented particle/background/neutralizing charge, deposited charge,
ordinary particle momentum/kinetic energy, shared electromagnetic field energy
(counted once) and their total. Particle velocities and B use the preceding fast
half level; E is at full time. This explicit mover does not conserve discrete
total energy exactly. Fixed populations have no velocity time offset.
"""
function electromagnetic_diagnostics(es::EMPIC{D,T}, species::AbstractVector) where {D,T}
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    work.initialized[] || throw(ArgumentError("initialize the electromagnetic collection before diagnostics"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    kinetic = kinetic_energy(species)
    momentum = total_momentum(species)
    field = em_field_energy(es)
    neutralizing = work.neutrality === :neutralize ? -charges.net : zero(charges.net)
    offset = es.step[] == 0 ? zero(T) : -es.last_sub_dt[] / 2
    return (
        particle_charge = charges.particle, background_charge = charges.background,
        neutralizing_charge = neutralizing, net_charge = charges.net + neutralizing,
        deposited_charge = _exact_float_sum(es.step[] == 0 ? es.rho_n : es.rho_np1) *
                           prod(_exact_rational, es.g.dx),
        momentum = momentum, kinetic = kinetic, field = field, total = kinetic + field,
        velocity_time_offsets = [mobile ? offset : zero(T) for mobile in work.mobile],
        magnetic_time_offset = offset,
    )
end
