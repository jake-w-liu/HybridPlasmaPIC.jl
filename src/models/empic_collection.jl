# Ordered periodic populations share the ordinary EMPIC Maxwell owner. The
# private ordinary state stages field publication; it does not define a solver.

"""
    EMPIC(g, species::AbstractVector; mobile=trues(length(species)),
          background_charge=0, neutrality=:require, charge_rtol=64eps(T),
          c=5, shape=CIC(), n_sub=1, species_stride=nothing, relativistic=false)

Periodic CPU Float32/Float64 electromagnetic PIC with arbitrary finite species
charges and positive masses. By default each mobile population kicks and drifts
at every Maxwell substep `dt/n_sub`. Explicit `species_stride` opts into dividing
integer kick schedules: mobile rows divide `n_sub`, fixed rows are zero. All
mobile populations still drift and deposit current at every fast substep. Use
[`set_empic_schedule!`](@ref) to change an opted-in future schedule. Fixed
populations require zero stored velocity, deposit charge only and remain unchanged.
NGP, CIC and TSC are supported. `relativistic=true` retains the existing prime/kick
and supported stored-output domain. Nondividing/asynchronous timesteps are unsupported.

`background_charge` is a fixed uniform charge density. `:require` uses the same
exact physical-charge tolerance as electrostatic collections; `:neutralize`
explicitly supplies a uniform compensator, reported by diagnostics. Neither
changes particle weights. Initialize with `init_empic!(state, species)` to bind
ordered population identities and physical full-time velocities. Replacing or
reordering populations requires initialization. Growth requires new velocities
at each species' existing half level (the common fast level by default).
Zero dt preserves that phase. Opted-in histories/schedules use collection schema2;
default synchronous collections retain schema1.

Particles, physical fields and clock publish after the complete operation
succeeds. Scratch and shared FFT work may change on rejection. Concurrent calls
require independent states, populations and grids. The integer-count constructor
retains its electron/ion contract and separate checkpoint layout.
All particle leaves, including fixed populations and metadata, must support
ordinary writes and expose aliases. Known read-only `Base.Experimental.Const`
and atomic-only storage reject through standard wrapper ancestry. OS protection
changes and arbitrary custom setter failures are outside the publication guarantee.
"""
function EMPIC(
    g::FourierGrid{D,T},
    species::AbstractVector;
    mobile::AbstractVector{Bool} = trues(length(species)),
    species_stride = nothing,
    background_charge = zero(T),
    neutrality::Symbol = :require,
    charge_rtol = 64eps(T),
    c = 5.0,
    shape::ShapeFunction = CIC(),
    n_sub::Integer = 1,
    relativistic::Bool = false,
) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("EMPIC supports D = 1, 2, or 3"))
    T in (Float32, Float64) ||
        throw(ArgumentError("electromagnetic collections require Float32 or Float64"))
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
    stage = _new_empic_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, nothing)
    P = eltype(species)
    schedule = species_stride === nothing ? nothing :
               _empic_species_schedule(mobile, ns, species_stride, T)
    work = _EMPICCollectionWork{D,T,P,typeof(stage),typeof(schedule)}(
        neutrality,
        rtol,
        BitVector(mobile),
        P[],
        _ElectrostaticHostParticle{D,T}[],
        NTuple{D,Vector{T}}[],
        stage,
        Ref(false),
        schedule,
    )
    es = _new_empic_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, work)
    _empic_collection_validate_species(es, species)
    return es
end

_require_scalar_empic(es::Union{EMPIC,EMPIC1D}) =
    es.collection === nothing ? nothing :
    throw(ArgumentError("electromagnetic collection state requires its ordered species vector"))

# Physical arrays carried across completed boundaries. Persistence, diagnostics
# and collection stepping share this inventory; gathers/derivatives are scratch.
_empic_physical_fields(es::EMPIC) =
    (E = es.E, B = es.B, J = es.J, rho_n = es.rho_n, rho_np1 = es.rho_np1)
_empic_physical_fields(es::EMPIC1D) =
    (Ex = es.Ex, Ey = es.Ey, Bz = es.Bz, Jx = es.Jx, Jy = es.Jy,
     rho_n = es.rho_n, rho_np1 = es.rho_np1)

function _empic_collection_require_finite_fields(es::Union{EMPIC,EMPIC1D})
    for components in values(_empic_physical_fields(es))
        finite =
            components isa Tuple ? all(a -> all(isfinite, a), components) :
            all(isfinite, components)
        finite || throw(ArgumentError("electromagnetic physical fields must be finite"))
    end
    return nothing
end

function _empic_collection(es::Union{EMPIC,EMPIC1D})
    work = es.collection
    work === nothing && throw(ArgumentError("species vectors require an electromagnetic collection state"))
    return work
end

function _empic_collection_bound(work, species)
    length(work.owners) == length(species) &&
        all(a === b for (a, b) in zip(work.owners, species)) || throw(
        ArgumentError("electromagnetic population membership/order changed; call init_empic!"),
    )
    return nothing
end

function _empic_collection_validate_storage(es::Union{EMPIC,EMPIC1D}, species)
    work = _empic_collection(es)
    es.mobile && throw(ArgumentError("collection EMPIC mobility is specified per population"))
    _espic_validate_collection_type(species, es.g)
    length(work.mobile) == length(species) ||
        throw(ArgumentError("electromagnetic population count must match mobility flags"))
    # The staged and public field inventories are deliberately identical to the
    # existing transaction inventory, including continuity diagnostic snapshots.
    fields = (_empic_transaction_fields(es)..., _empic_transaction_fields(work.stage)...)
    work.stage.g === es.g ||
        throw(ArgumentError("electromagnetic stage must share its owning grid"))
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
    schedule_arrays = _empic_collection_schedule_arrays(work.schedule)
    buffers = (
        schedule_arrays...,
        fields...,
        fft...,
        metadata...,
        _empic_electron_gathers(es)...,
        _empic_ion_gathers(es)...,
        _empic_electron_gathers(work.stage)...,
        _empic_ion_gathers(work.stage)...,
    )
    gathers = (
        _empic_electron_gathers(es)...,
        _empic_ion_gathers(es)...,
        _empic_electron_gathers(work.stage)...,
        _empic_ion_gathers(work.stage)...,
    )
    for (i, a) in enumerate(gathers)
        for b in (fields..., fft..., metadata..., gathers[1:(i-1)]...)
            Base.mightalias(a, b) &&
                throw(ArgumentError("electromagnetic gather storage must be independent"))
        end
    end
    for (i, a) in enumerate(schedule_arrays)
        axes(a, 1) == Base.OneTo(length(work.mobile)) ||
            throw(ArgumentError("electromagnetic schedule storage must have one row per population"))
        for b in (fields..., fft..., metadata..., gathers..., schedule_arrays[1:(i-1)]...)
            Base.mightalias(a, b) &&
                throw(ArgumentError("electromagnetic schedule storage must be independent"))
        end
    end
    clocks = (
        es.time,
        es.step,
        es.last_dt,
        es.last_sub_dt,
        work.stage.time,
        work.stage.step,
        work.stage.last_dt,
        work.stage.last_sub_dt,
    )
    for (i, a) in enumerate(clocks), j = 1:(i-1)
        a === clocks[j] && throw(ArgumentError("electromagnetic clock storage must be independent"))
    end
    for (i, ps) in enumerate(species)
        _check_particle_storage_axes(ps)
        _require_particle_ordinary_stores(ps, "electromagnetic collection")
        all(a -> a isa StridedVector && _particle_batch_injective(a),
            (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)) ||
            throw(ArgumentError("electromagnetic collections require injective CPU strided particle storage"))
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

function _empic_collection_validate_populations(
    es::Union{EMPIC,EMPIC1D},
    species,
    c::T,
    relativistic::Bool,
) where {T}
    _empic_collection_validate_storage(es, species)
    work = _empic_collection(es)
    _require_finite_real("background_charge", es.n0, T)
    for (i, ps) in enumerate(species)
        _validated_boris_factor(ps)
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "electromagnetic collection")
        _empic_collection_validate_polarization(es, ps, work.mobile[i])
        relativistic && _require_subluminal_velocities(ps, c, "electromagnetic collection")
        all(w -> isfinite(w) && _stored_nonnegative(w), ps.weight) ||
            throw(ArgumentError("electromagnetic weights must be finite and nonnegative"))
        !work.mobile[i] &&
            !all(a -> all(_particle_stored_zero, a), ps.v) &&
            throw(
                ArgumentError("fixed electromagnetic populations require original stored zero velocity"),
            )
    end
    charges = _pic_collection_charges(es.g, es.n0, species)
    _espic_check_neutrality(work, charges)
    return charges
end

function _empic_collection_validate_species(es::Union{EMPIC,EMPIC1D}, species)
    T = eltype(es.rho_n)
    _require_finite_positive_real("c", es.c, T)
    es.n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    _empic_collection_validate_schedule(es, _empic_collection(es).schedule)
    return _empic_collection_validate_populations(es, species, es.c, es.relativistic)
end

function _empic_collection_prepare!(es::Union{EMPIC,EMPIC1D}, species)
    T = eltype(es.rho_n)
    work = _empic_collection(es)
    _pic_prepare_species_stages!(work.staged, species)
    while length(work.midpoints) < length(species)
        push!(work.midpoints, ntuple(_ -> T[], length(es.g.n)))
    end
    resize!(work.midpoints, length(species))
    for (ps, mid) in zip(species, work.midpoints), a in mid
        resize!(a, nparticles(ps))
    end
    stage = work.stage
    stage.n0, stage.c, stage.shape, stage.n_sub = es.n0, es.c, es.shape, es.n_sub
    stage.relativistic = es.relativistic
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
    _charge_density_sum!(rho, stage.ne, work.staged, stage.g, stage.shape; background = stage.n0)
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
    _empic_collection_publish_intervals!(work.schedule)
    return es
end

_empic_collection_seed_fields(es::EMPIC) = (es.E..., es.B...)
_empic_collection_seed_fields(es::EMPIC1D) = (es.Ex, es.Ey, es.Bz)

function init_empic!(es::Union{EMPIC,EMPIC1D}, species::AbstractVector)
    T = eltype(es.rho_n)
    work = _empic_collection(es)
    all(ps -> ps isa eltype(work.owners), species) ||
        throw(ArgumentError("population storage type changed; construct a new electromagnetic collection state"))
    charges = _empic_collection_validate_species(es, species)
    all(a -> all(isfinite, a), _empic_collection_seed_fields(es)) ||
        throw(ArgumentError("electromagnetic fields must be finite"))
    _empic_collection_prepare!(es, species)
    sizehint!(work.owners, length(species))
    stage = work.stage
    _empic_collection_charge!(
        stage.rho_n,
        stage,
        work,
        _empic_collection_charge_target(es, work, charges),
    )
    _empic_collection_initialize_fields!(stage)
    stage.time[], stage.step[] = zero(T), 0
    stage.last_dt[], stage.last_sub_dt[] = zero(T), zero(T)
    _empic_collection_reset_intervals!(work.schedule)
    _empic_collection_publish!(es, work; particles = false)
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = true
    return es
end

function _empic_collection_initialize_fields!(stage::EMPIC{D,T}) where {D,T}
    _enforce_gauss!(stage)
    all(a -> all(isfinite, a), stage.E) ||
        throw(ArgumentError("electromagnetic initialized field exceeds numeric range"))
    copyto!(stage.rho_np1, stage.rho_n)
    foreach(a -> fill!(a, zero(T)), (stage.J..., stage.Jraw...))
    foreach(copyto!, stage.Bc, stage.B)
    return stage
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
            _empic_collection_exact_current!(stage.Jraw[c], stage, work, c)
        end
    end
    all(a -> all(isfinite, a), stage.Jraw) ||
        throw(ArgumentError("electromagnetic deposited current exceeds numeric range"))
    return stage.Jraw
end

function _empic_collection_exact_current!(out::Array{T}, stage, work, c) where {T}
    exact = fill(zero(Rational{BigInt}), stage.g.n)
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        contribution =
            _exact_midpoint_current(ps, work.midpoints[i], stage.g, stage.shape, c)
        for I in eachindex(exact)
            exact[I] += contribution[I]
        end
    end
    for I in eachindex(exact)
        out[I] = _from_exact(T, exact[I])
    end
    return out
end

function _empic_collection_prime!(stage::EMPIC{D,T}, work, sub_dt::T) where {D,T}
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_empic_electron_buffers!(stage, nparticles(ps))
        gather_vector!(stage.Ep, stage.E, ps, stage.g, stage.shape)
        gather_vector!(stage.Bp, stage.B, ps, stage.g, stage.shape)
        _prime_leapfrog!(
            ps.v,
            stage.Ep,
            stage.Bp,
            _validated_boris_factor(ps),
            _empic_collection_prime_interval(work.schedule, i, sub_dt) / 2,
            nparticles(ps),
            stage.relativistic,
            stage.c,
        )
    end
    _empic_prime_magnetic!(stage, stage.g, sub_dt)
    return stage
end

function _empic_collection_substep!(stage::EMPIC{D,T}, work, dt::T, previous::T, target, index) where {D,T}
    _empic_collection_charge!(stage.rho_n, stage, work, target)
    parts = _empic_center_magnetic!(stage, stage.g, previous)
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        if _empic_collection_kick_due(work.schedule, i, index)
            _ensure_empic_electron_buffers!(stage, nparticles(ps))
            gather_vector!(stage.Ep, stage.E, ps, stage.g, stage.shape)
            gather_vector!(stage.Bp, stage.Bc, ps, stage.g, stage.shape)
            interval = _empic_collection_kick_interval(work.schedule, i, previous, dt)
            _push_velocities!(stage, ps, stage.Ep, stage.Bp, interval)
            _empic_collection_record_kick!(work.schedule, i)
        end
        _drift_record_mid!(ps, work.midpoints[i], dt, stage.g.L)
        apply_periodic!(ps, ntuple(_ -> zero(T), Val(D)), stage.g.L)
        _require_finite_particle_velocities(ps, "electromagnetic collection output")
        stage.relativistic &&
            _require_subluminal_velocities(ps, stage.c, "electromagnetic collection output")
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

function _empic_collection_validate_phase(es::Union{EMPIC,EMPIC1D})
    T = eltype(es.rho_n)
    isfinite(es.time[]) && es.time[] >= zero(T) && es.step[] >= 0 ||
        throw(ArgumentError("electromagnetic clock must be finite and nonnegative"))
    previous = es.last_sub_dt[]
    full = es.last_dt[]
    isfinite(full) &&
        isfinite(previous) &&
        zero(T) <= previous <= full &&
        (es.step[] == 0 ? iszero(full) && iszero(previous) : previous > zero(T)) ||
        throw(ArgumentError("electromagnetic timestep history is inconsistent with its phase"))
    _empic_collection_validate_species_history(es, _empic_collection(es).schedule)
    return nothing
end

function step_empic!(es::Union{EMPIC,EMPIC1D}, species::AbstractVector, dt::Real)
    T = eltype(es.rho_n)
    dtT =
        _require_resolved_dt("step_empic!", _validated_nonnegative_dt(T, dt; name = "step_empic!"))
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    work.initialized[] ||
        throw(ArgumentError("initialize the electromagnetic collection before stepping"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    iszero(dtT) && return es
    sub_dt = _require_resolved_dt("step_empic! substep dt/n_sub", dtT / es.n_sub)
    sub_dt > zero(T) || throw(ArgumentError("dt/n_sub must be representable and positive"))
    next_time, next_step = _validated_step_clock(es.time[], es.step[], sub_dt, es.n_sub)
    _empic_collection_advance_candidate!(es, species, work, charges, dtT, sub_dt, next_time, next_step)
    return _empic_collection_publish!(es, work; particles = true)
end

function _empic_collection_advance_candidate!(es, species, work, charges, dtT, sub_dt, next_time, next_step)
    _empic_collection_prepare_intervals!(work.schedule, work.mobile, es, dtT, sub_dt)
    _empic_collection_prepare!(es, species)
    stage = work.stage
    stage.step[] == 0 && _empic_collection_prime!(stage, work, sub_dt)
    previous = stage.step[] == 0 ? sub_dt : stage.last_sub_dt[]
    # q and weights are unchanged by canonical kicks/drifts. Reuse their exact
    # conserved total rather than re-summing all particles at both endpoints.
    target = _empic_collection_charge_target(es, work, charges)
    for index = 1:es.n_sub
        _empic_collection_substep!(stage, work, sub_dt, previous, target, index)
        previous = sub_dt
    end
    stage.time[], stage.step[] = next_time, next_step
    stage.last_dt[], stage.last_sub_dt[] = dtT, sub_dt
    return stage
end

"""
    electromagnetic_diagnostics(state::EMPIC, species)

Exact represented particle/background/neutralizing charge, deposited charge,
particle momentum/kinetic energy (relativistic when configured), shared electromagnetic field energy
(counted once) and their total. Particle velocities use each preceding declared
species half interval and B uses the fast half interval; E is at full time.
Offsets are nominal stage intervals, not exact timestamps relative to the rounded
public clock. This explicit mover does not conserve discrete total energy exactly. Fixed populations have no velocity time offset.
"""
function electromagnetic_diagnostics(es::Union{EMPIC,EMPIC1D}, species::AbstractVector)
    T = eltype(es.rho_n)
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    work.initialized[] ||
        throw(ArgumentError("initialize the electromagnetic collection before diagnostics"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    kinetic, momentum = es.relativistic ? _relativistic_particle_moments(species, es.c) :
        (kinetic_energy(species), total_momentum(species))
    field = em_field_energy(es)
    neutralizing = work.neutrality === :neutralize ? -charges.net : zero(charges.net)
    offset = es.step[] == 0 ? zero(T) : -es.last_sub_dt[] / 2
    result = (
        particle_charge = charges.particle,
        background_charge = charges.background,
        neutralizing_charge = neutralizing,
        net_charge = charges.net + neutralizing,
        deposited_charge = _exact_float_sum(es.step[] == 0 ? es.rho_n : es.rho_np1) *
                           prod(_exact_rational, es.g.dx),
        momentum = momentum,
        kinetic = kinetic,
        field = field,
        total = kinetic + field,
        velocity_time_offsets = _empic_collection_velocity_offsets(work.schedule, work.mobile, offset),
        magnetic_time_offset = offset,
    )
    return _empic_collection_schedule_diagnostics(work.schedule, result)
end

"""
    EMPIC1D(g, species::AbstractVector; mobile=trues(length(species)),
            background_charge=0, neutrality=:require, charge_rtol=64eps(T),
            c=5, shape=CIC(), n_sub=1, species_stride=nothing, relativistic=false)

Specialized periodic CPU Float32/64 PIC with fields `Ex/Ey/Bz` and ordered
arbitrary finite species charges and positive masses. Every original stored
`v_z` must be zero: the omitted `Ez/By` polarization is not evolved. Fixed
populations require zero velocity and deposit charge only. Every mobile
population drifts and contributes current on every `h=dt/n_sub`; its integer
`species_stride` divides `n_sub` and selects kicks at fast indices `0,r,2r,...`.
Default mobile strides are 1 and fixed strides are 0. The species kick interval
is `H=dt/(n_sub÷r)`, with its own previous-half history. Use
[`set_empic1d_schedule!`](@ref) between completed external operations. NGP, CIC
and TSC are supported.

The background/neutrality, identity binding, variable fast phase and physical
moment diagnostics follow the collection API of [`EMPIC`](@ref). Initialize
with physical full-time velocities and seeded `Ey/Bz`; initialization solves
`Ex` with the specialized zero-DC/Nyquist convention. Expected input/arithmetic
failures publish no particles, physical fields or clock. Scratch may change.
Calls require independent states, species and grids for concurrency. The
integer-count constructor and its separate electron/ion payload are unchanged.
"""
function EMPIC1D(
    g::FourierGrid{1,T},
    species::AbstractVector;
    mobile::AbstractVector{Bool} = trues(length(species)),
    species_stride = nothing,
    background_charge = zero(T),
    neutrality::Symbol = :require,
    charge_rtol = 64eps(T),
    c = 5.0,
    shape::ShapeFunction = CIC(),
    n_sub::Integer = 1,
    relativistic::Bool = false,
) where {T}
    T in (Float32, Float64) ||
        throw(ArgumentError("electromagnetic collections require Float32 or Float64"))
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
    stage = _new_empic1d_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, nothing)
    P = eltype(species)
    schedule = _empic_species_schedule(mobile, ns,
        species_stride === nothing ? [active ? 1 : 0 for active in mobile] : species_stride, T)
    work = _EMPICCollectionWork{1,T,P,typeof(stage),typeof(schedule)}(
        neutrality, rtol, BitVector(mobile), P[], _ElectrostaticHostParticle{1,T}[],
        NTuple{1,Vector{T}}[], stage, Ref(false), schedule,
    )
    es = _new_empic1d_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, work)
    _empic_collection_validate_species(es, species)
    return es
end

_empic_collection_validate_polarization(::EMPIC, ps, mobile) = nothing
function _empic_collection_validate_polarization(::EMPIC1D, ps, mobile)
    mobile && !all(_particle_stored_zero, ps.v[3]) &&
        throw(ArgumentError("EMPIC1D collections require original stored v_z = 0"))
    return nothing
end

function _empic_collection_initialize_fields!(stage::EMPIC1D{T}) where {T}
    g = stage.g
    _empic1d_gauss!(stage, g, g.n[1], g.cbuf)
    all(a -> all(isfinite, a), (stage.Ex, stage.Ey)) ||
        throw(ArgumentError("electromagnetic initialized field exceeds numeric range"))
    copyto!(stage.rho_np1, stage.rho_n)
    fill!(stage.Jx, zero(T))
    fill!(stage.Jy, zero(T))
    copyto!(stage.Bzc, stage.Bz)
    return stage
end

function _empic_collection_prime!(stage::EMPIC1D{T}, work, sub_dt::T) where {T}
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_electron_buffers!(stage, nparticles(ps))
        gather_scalar!(stage.Exp, stage.Ex, ps, stage.g, stage.shape)
        gather_scalar!(stage.Eyp, stage.Ey, ps, stage.g, stage.shape)
        gather_scalar!(stage.Bzp, stage.Bz, ps, stage.g, stage.shape)
        _prime_empic1d!(ps.v, stage.Exp, stage.Eyp, stage.Bzp,
            _validated_boris_factor(ps), _empic_collection_prime_interval(work.schedule, i, sub_dt) / 2,
            nparticles(ps), stage.relativistic, stage.c)
    end
    _empic1d_prime_magnetic!(stage, stage.g, sub_dt)
    return stage
end

function _empic1d_collection_Jx0(stage::EMPIC1D{T}, work) where {T}
    value = zero(T)
    eligible, lost = true, false
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        state = _species_Jx0_with_range(ps, stage.g.L[1])
        eligible &= state.eligible
        value += state.value
        lost |= state.lost || _deposition_range_lost(value)
    end
    if eligible && lost
        exact = zero(Rational{BigInt})
        for (i, ps) in enumerate(work.staged)
            work.mobile[i] || continue
            exact += _exact_species_Jx0(ps, stage.g.L[1])
        end
        value = _from_exact(T, exact)
    end
    isfinite(value) || throw(ArgumentError("electromagnetic mean current exceeds numeric range"))
    return value
end

function _empic_collection_current!(stage::EMPIC1D{T}, work) where {T}
    fill!(stage.Jy, zero(T))
    eligible, lost, accumulate = true, false, false
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        _ensure_electron_buffers!(stage, nparticles(ps))
        mid = work.midpoints[i]
        _swap_positions_with_midpoints!(ps, mid)
        state = try
            _add_transverse_Jy!(stage, ps, stage.Eyp, stage.dBz, accumulate)
        finally
            _swap_positions_with_midpoints!(ps, mid)
        end
        accumulate = true
        eligible &= state.eligible
        lost |= state.lost || any(_deposition_range_lost, stage.Jy)
    end
    if eligible && lost
        _empic_collection_exact_current!(stage.Jy, stage, work, 2)
    end
    all(isfinite, stage.Jy) ||
        throw(ArgumentError("electromagnetic deposited current exceeds numeric range"))
    _zero_nyquist_mode!(stage.Jy)
    return stage.Jy
end

function _empic_collection_substep!(stage::EMPIC1D{T}, work, dt::T, previous::T, target, index) where {T}
    g = stage.g
    _empic_collection_charge!(stage.rho_n, stage, work, target)
    parts = _empic1d_center_magnetic!(stage, g, previous, zero(T))
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        if _empic_collection_kick_due(work.schedule, i, index)
            _ensure_electron_buffers!(stage, nparticles(ps))
            gather_scalar!(stage.Exp, stage.Ex, ps, g, stage.shape)
            gather_scalar!(stage.Eyp, stage.Ey, ps, g, stage.shape)
            gather_scalar!(stage.Bzp, stage.Bzc, ps, g, stage.shape)
            interval = _empic_collection_kick_interval(work.schedule, i, previous, dt)
            _push_velocities!(stage, ps, stage.Exp, stage.Eyp, stage.Bzp, interval)
            _empic_collection_record_kick!(work.schedule, i)
        end
        _drift_record_mid!(ps, work.midpoints[i][1], dt, g.L[1])
        apply_periodic!(ps, (zero(T),), g.L)
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "electromagnetic collection output")
        stage.relativistic &&
            _require_subluminal_velocities(ps, stage.c, "electromagnetic collection output")
    end
    _empic_collection_current!(stage, work)
    _empic_collection_charge!(stage.rho_np1, stage, work, target)
    _esirkepov_Jx!(stage, _empic1d_collection_Jx0(stage, work), dt)
    _empic1d_update_fields!(stage, g, g.n[1], stage.c * stage.c, parts, previous, dt)
    _empic_collection_require_finite_fields(stage)
    all(isfinite, stage.Bzc) ||
        throw(ArgumentError("electromagnetic centered magnetic field exceeds numeric range"))
    return stage
end


_empic_collection_schedule_arrays(::Nothing) = ()
_empic_collection_schedule_arrays(schedule::_EMPICSpeciesSchedule) =
    (schedule.stride, schedule.last_species_dt, schedule.new_species_dt, schedule.previous_species_dt)

function _empic_validate_species_schedule(mobile, M, values)
    M isa Int && M >= 1 || throw(ArgumentError("n_sub must be a positive Int"))
    values isa AbstractVector && axes(values, 1) == Base.OneTo(length(mobile)) ||
        throw(ArgumentError("species_stride must have one integer row per population"))
    for (active, value) in zip(mobile, values)
        value isa Integer && !(value isa Bool) && 0 <= value <= typemax(Int) ||
            throw(ArgumentError("species strides must be nonnegative integers fitting Int"))
        if active
            1 <= value <= M && M % value == 0 ||
                throw(ArgumentError("mobile species strides must be positive divisors of n_sub"))
        else
            iszero(value) || throw(ArgumentError("fixed species strides must be zero"))
        end
    end
    return nothing
end

function _empic_species_stride_values(mobile, M, values)
    _empic_validate_species_schedule(mobile, M, values)
    return Int[value for value in values]
end

function _empic_species_schedule(mobile, M, values, ::Type{T}) where {T}
    stride = _empic_species_stride_values(mobile, M, values)
    return _EMPICSpeciesSchedule{T}(stride, zeros(T, length(mobile)),
        zeros(T, length(mobile)), zeros(T, length(mobile)))
end

_empic_collection_prime_interval(::Nothing, i, fast) = fast
_empic_collection_prime_interval(schedule::_EMPICSpeciesSchedule, i, fast) = schedule.new_species_dt[i]
_empic_collection_kick_due(::Nothing, i, index) = true
_empic_collection_kick_due(schedule::_EMPICSpeciesSchedule, i, index) =
    (index - 1) % schedule.stride[i] == 0
_empic_collection_kick_interval(::Nothing, i, previous, fast) = previous / 2 + fast / 2
_empic_collection_kick_interval(schedule::_EMPICSpeciesSchedule, i, previous, fast) =
    schedule.previous_species_dt[i] / 2 + schedule.new_species_dt[i] / 2
_empic_collection_record_kick!(::Nothing, i) = nothing
function _empic_collection_record_kick!(schedule::_EMPICSpeciesSchedule, i)
    schedule.previous_species_dt[i] = schedule.new_species_dt[i]
    return nothing
end

function _empic_validate_species_history(mobile, history, step, full, fast, ::Type{T}) where {T}
    history isa AbstractVector{T} && axes(history, 1) == Base.OneTo(length(mobile)) ||
        throw(ArgumentError("species timestep history must match population count and precision"))
    if step > 0
        isfinite(fast) && fast >= 2 * floatmin(T) ||
            throw(ArgumentError("previous fast interval must be finite and resolved"))
    end
    for (active, h) in zip(mobile, history)
        if step == 0 || !active
            _particle_stored_zero(h) ||
                throw(ArgumentError("step-zero and fixed species histories must be stored zero"))
        else
            isfinite(h) && h >= 2 * floatmin(T) && fast <= h <= full ||
                throw(ArgumentError("species timestep history must lie between previous fast and full intervals"))
        end
    end
    return nothing
end

_empic_collection_validate_schedule(es, ::Nothing) = nothing
_empic_collection_validate_schedule(es, schedule::_EMPICSpeciesSchedule) =
    _empic_validate_species_schedule(es.collection.mobile, es.n_sub, schedule.stride)
_empic_collection_validate_species_history(es, ::Nothing) = nothing
_empic_collection_validate_species_history(es, schedule::_EMPICSpeciesSchedule{T}) where {T} =
    _empic_validate_species_history(es.collection.mobile, schedule.last_species_dt,
        es.step[], es.last_dt[], es.last_sub_dt[], T)

_empic_collection_prepare_intervals!(::Nothing, mobile, es, dt, fast) = nothing
function _empic_collection_prepare_intervals!(schedule::_EMPICSpeciesSchedule{T}, mobile, es, dt, fast) where {T}
    for (i, active) in enumerate(mobile)
        h = active ? dt / (es.n_sub ÷ schedule.stride[i]) : zero(T)
        if active
            isfinite(h) && fast <= h <= dt && h > zero(T) ||
                throw(ArgumentError("species drift interval must lie between fast and full intervals"))
            _require_resolved_dt("step_empic! species interval", h)
        end
        schedule.new_species_dt[i] = h
        schedule.previous_species_dt[i] = es.step[] == 0 ? h : schedule.last_species_dt[i]
    end
    return nothing
end

_empic_collection_reset_intervals!(::Nothing) = nothing
function _empic_collection_reset_intervals!(schedule::_EMPICSpeciesSchedule{T}) where {T}
    fill!(schedule.new_species_dt, zero(T))
    fill!(schedule.previous_species_dt, zero(T))
    return nothing
end
_empic_collection_publish_intervals!(::Nothing) = nothing
_empic_collection_publish_intervals!(schedule::_EMPICSpeciesSchedule) =
    copyto!(schedule.last_species_dt, schedule.previous_species_dt)

_empic_collection_velocity_offsets(::Nothing, mobile, offset) =
    [active ? offset : zero(offset) for active in mobile]
_empic_collection_velocity_offsets(schedule::_EMPICSpeciesSchedule{T}, mobile, offset) where {T} =
    [active && !iszero(h) ? -h / 2 : zero(T) for (active, h) in zip(mobile, schedule.last_species_dt)]
_empic_collection_schedule_diagnostics(::Nothing, result) = result
_empic_collection_schedule_diagnostics(schedule::_EMPICSpeciesSchedule, result) =
    merge(result, (species_stride = copy(schedule.stride), last_species_dt = copy(schedule.last_species_dt)))

"""
    set_empic1d_schedule!(state; n_sub=state.n_sub, species_stride=nothing)

Set a specialized collection's future dividing integer kick strides at a completed
external boundary. `nothing` keeps its current strides. Mobile rows divide
`n_sub`; fixed rows are zero. The complete proposal is copied and validated
before configuration changes. Previous species and magnetic half intervals are
preserved; the next positive step joins the old and new species half levels.
The default all-fast state can be changed before or after initialization.
"""
set_empic1d_schedule!(es::EMPIC1D; n_sub::Integer = es.n_sub, species_stride = nothing) =
    _set_empic_collection_schedule!(es, n_sub, species_stride)

"""
    set_empic_schedule!(state::EMPIC; n_sub=state.n_sub, species_stride=nothing)

Set future dividing species kick strides on an opted-in generic collection.
Construct with explicit `species_stride` to opt in; a default synchronous state
keeps its concrete Nothing schedule owner. The complete proposal is copied and
validated before either configuration value changes. Prior species/fast histories
are preserved, and `nothing` keeps current strides. Call only at completed
external boundaries with no concurrent use.
"""
set_empic_schedule!(es::EMPIC; n_sub::Integer = es.n_sub, species_stride = nothing) =
    _set_empic_collection_schedule!(es, n_sub, species_stride)

function _set_empic_collection_schedule!(es, n_sub, species_stride)
    work = _empic_collection(es)
    work.schedule isa _EMPICSpeciesSchedule || throw(
        ArgumentError("construct EMPIC with explicit species_stride before setting a species schedule"),
    )
    1 <= n_sub <= typemax(Int) || throw(ArgumentError("n_sub must be positive and fit Int"))
    ns = Int(n_sub)
    values = species_stride === nothing ? work.schedule.stride : species_stride
    stride = _empic_species_stride_values(work.mobile, ns, values)
    _empic_collection_validate_phase(es)
    es.n_sub = ns
    work.schedule.stride = stride
    return es
end
