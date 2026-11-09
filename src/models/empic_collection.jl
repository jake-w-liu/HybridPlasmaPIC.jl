# Ordered periodic populations share the ordinary EMPIC Maxwell owner. The
# private ordinary state stages field publication; it does not define a solver.

const _EMPICEventQ = Rational{BigInt}
struct _EMPICEventPhase{T}
    time::_EMPICEventQ
    primed::Bool
    interval::Vector{T}
    last::Vector{_EMPICEventQ}
    next::Vector{_EMPICEventQ}
    report::Any
end
mutable struct _EMPICEventSchedule{T}
    proposal::NamedTuple
    phase::_EMPICEventPhase{T}
    candidate::_EMPICEventPhase{T}
    shadow::Any
    candidate_shadow::Any
end

function _empic_event_zero_phase(::Type{T}, count) where {T}
    return _EMPICEventPhase{T}(
        zero(_EMPICEventQ),
        false,
        zeros(T, count),
        fill(zero(_EMPICEventQ), count),
        fill(zero(_EMPICEventQ), count),
        nothing,
    )
end
function _empic_event_scalar(name, value, ::Type{T}) where {T}
    value isa Real && isfinite(value) || throw(ArgumentError("$name must be finite and real"))
    result = _require_finite_real(name, value, T)
    (iszero(value) || !iszero(result)) || throw(ArgumentError("$name disappears during conversion"))
    return result
end
function _empic_event_proposal(mobile, cap, adaptive, minimum, limits, budget, ::Type{T}) where {T}
    count = length(mobile)
    cap isa AbstractVector && axes(cap, 1) == Base.OneTo(count) ||
        throw(ArgumentError("species_dt must have one row per population"))
    active = adaptive === nothing ? falses(count) : adaptive
    active isa AbstractVector{Bool} && axes(active, 1) == Base.OneTo(count) ||
        throw(ArgumentError("species_adaptive must have one Boolean row per population"))
    low = minimum === nothing ? zeros(T, count) : minimum
    low isa AbstractVector && axes(low, 1) == Base.OneTo(count) ||
        throw(ArgumentError("species_dt_min must have one row per population"))
    budget isa Integer && !(budget isa Bool) && 1 <= budget <= typemax(Int) ||
        throw(ArgumentError("persistent schedules require a positive Int event_budget"))
    upper = T[_empic_event_scalar("species_dt", value, T) for value in cap]
    lower = T[_empic_event_scalar("species_dt_min", value, T) for value in low]
    flags = BitVector(active)
    for i in eachindex(mobile)
        if mobile[i]
            upper[i] > zero(T) || throw(ArgumentError("mobile species_dt must be positive"))
            _require_resolved_dt("species_dt", upper[i])
            if flags[i]
                minimum === nothing && throw(ArgumentError("adaptive rows require species_dt_min"))
                zero(T) < lower[i] <= upper[i] || throw(
                    ArgumentError("adaptive minimum must be positive and no larger than its cap"),
                )
                _require_resolved_dt("species_dt_min", lower[i])
            else
                _particle_stored_zero(lower[i]) ||
                    throw(ArgumentError("nonadaptive species_dt_min rows must be stored zero"))
            end
        else
            _particle_stored_zero(upper[i]) && _particle_stored_zero(lower[i]) && !flags[i] ||
                throw(ArgumentError("fixed species require zero intervals and no adaptation"))
        end
    end
    policy = if limits === nothing
        any(flags) && throw(ArgumentError("adaptive rows require explicit schedule_limits"))
        (rotation = zero(T), electric = zero(T), travel = zero(T), growth = one(T))
    else
        limits isa NamedTuple && keys(limits) == (:rotation, :electric, :travel, :growth) ||
            throw(ArgumentError("schedule_limits requires rotation/electric/travel/growth"))
        result = NamedTuple{keys(limits)}(map(values(limits)) do value
            _empic_event_scalar("schedule limit", value, T)
        end)
        (
            (all(>(zero(T)), values(result)) && result.growth >= one(T)) || (
                !any(flags) && isequal(
                    result,
                    (rotation = zero(T), electric = zero(T), travel = zero(T), growth = one(T)),
                )
            )
        ) || throw(ArgumentError("schedule limits must be positive, with growth at least one"))
        result
    end
    return (cap = upper, adaptive = flags, minimum = lower, limits = policy, budget = Int(budget))
end
function _empic_collection_schedule(
    mobile,
    ns,
    stride,
    cap,
    adaptive,
    minimum,
    limits,
    budget,
    ::Type{T},
    error_control = nothing,
) where {T}
    error_control === nothing || return _empic_coupled_schedule(
        mobile,
        ns,
        stride,
        cap,
        adaptive,
        minimum,
        limits,
        budget,
        T,
        error_control,
    )
    if cap === nothing
        all(isnothing, (adaptive, minimum, limits, budget)) ||
            throw(ArgumentError("persistent schedule options require species_dt"))
        return stride === nothing ? nothing : _empic_species_schedule(mobile, ns, stride, T)
    end
    stride === nothing ||
        throw(ArgumentError("species_stride and species_dt are mutually exclusive"))
    proposal = _empic_event_proposal(mobile, cap, adaptive, minimum, limits, budget, T)
    ns <= proposal.budget ||
        throw(ArgumentError("base field partition exceeds persistent event budget"))
    return _EMPICEventSchedule{T}(
        proposal,
        _empic_event_zero_phase(T, length(mobile)),
        _empic_event_zero_phase(T, length(mobile)),
        nothing,
        nothing,
    )
end

"""
    EMPIC(g, species::AbstractVector; mobile=trues(length(species)),
          background_charge=0, neutrality=:require, charge_rtol=64eps(T),
          c=5, shape=CIC(), n_sub=1, species_stride=nothing,
          species_dt=nothing, event_budget=nothing, error_control=nothing,
          relativistic=false)

Periodic CPU Float32/Float64 electromagnetic PIC with arbitrary finite species
charges and positive masses. By default each mobile population kicks and drifts
at every Maxwell substep `dt/n_sub`. Explicit `species_stride` opts into dividing
integer kick schedules: mobile rows divide `n_sub`, fixed rows are zero. All
mobile populations still drift and deposit current at every fast substep. Use
[`set_empic_schedule!`](@ref) to change an opted-in future schedule. Fixed
populations require zero stored velocity, deposit charge only and remain unchanged.
NGP, CIC and TSC are supported. `relativistic=true` retains the existing prime/kick
and supported stored-output domain. Explicit `species_dt` with `event_budget`
opts into persistent nondividing kick deadlines and is mutually exclusive with
`species_stride`. All mobile rows still drift and supply current on every actual
field segment. `species_adaptive`, `species_dt_min` and explicit `schedule_limits`
enable a sampled-state resolution heuristic; it is not coupled error control.

Explicit `error_control` with a persistent calendar selects the closed full-time
projected Euler/Heun controller described in `docs/EMPIC_COUPLED_ERROR_CONTROL.md`.
That policy supports CIC/TSC knot crossings, evaluates all mobile forces at both
trial stages and uses schema4/policy3. NGP rejects at this opt-in boundary. Its
measured embedded defect supplies no universal true-error or stability theorem.

`background_charge` is a fixed uniform charge density. `:require` uses the same
exact physical-charge tolerance as electrostatic collections; `:neutralize`
explicitly supplies a uniform compensator, reported by diagnostics. Neither
changes particle weights. Initialize with `init_empic!(state, species)` to bind
ordered population identities and physical full-time velocities. Replacing or
reordering populations requires initialization. Growth requires new velocities
at each species' existing half level (the common fast level by default).
Zero dt preserves that phase. Persistent schedules fix counts, q/m, weights,
provenance, component identities and committed physical values between
initializations; prepare physical full-time inputs before reinitializing after
external edits. They use schema3/policy2 and nominal velocity offsets
`(last_kick+next_kick)/2-source_time`. Dividing schedules use schema2; default
synchronous collections retain schema1. No implicit phase migration is supplied.

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
    species_dt = nothing,
    species_adaptive = nothing,
    species_dt_min = nothing,
    schedule_limits = nothing,
    event_budget = nothing,
    error_control = nothing,
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
    error_control === nothing ||
        shape isa Union{CIC,TSC} ||
        throw(
            ArgumentError(
                "coupled error control requires CIC or TSC; NGP jump control is not implemented",
            ),
        )
    neutrality in (:require, :neutralize) ||
        throw(ArgumentError("neutrality must be :require or :neutralize"))
    _espic_validate_collection_type(species, g)
    axes(mobile, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("mobile must have one Boolean per species with one-based axes"))
    n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    species_dt === nothing ||
        (!(n_sub isa Bool) && n_sub <= typemax(Int)) ||
        throw(ArgumentError("persistent n_sub must be a positive fitting Int"))
    ns = Int(n_sub)
    background = _require_finite_real("background_charge", background_charge, T)
    rtol = _require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
    cT = _require_finite_positive_real("c", c, T)
    stage = _new_empic_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, nothing)
    P = eltype(species)
    schedule = _empic_collection_schedule(
        mobile,
        ns,
        species_stride,
        species_dt,
        species_adaptive,
        species_dt_min,
        schedule_limits,
        event_budget,
        T,
        error_control,
    )
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
    _empic_coupled_bind!(es, schedule)
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
_empic_physical_fields(es::EMPIC1D) = (
    Ex = es.Ex,
    Ey = es.Ey,
    Bz = es.Bz,
    Jx = es.Jx,
    Jy = es.Jy,
    rho_n = es.rho_n,
    rho_np1 = es.rho_np1,
)

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
    work === nothing &&
        throw(ArgumentError("species vectors require an electromagnetic collection state"))
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
        _empic_collection_schedule_axes(work.schedule, a, length(work.mobile))
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
        all(
            a -> (a isa StridedVector && _particle_batch_injective(a))::Bool,
            (ps.x..., ps.v..., ps.weight, ps.id, ps.tag),
        ) || throw(
            ArgumentError(
                "electromagnetic collections require injective CPU strided particle storage",
            ),
        )
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
    return _empic_collection_validate_population_values(es, species, c, relativistic)
end

function _empic_collection_validate_population_values(
    es::Union{EMPIC,EMPIC1D},
    species,
    c::T,
    relativistic::Bool,
) where {T}
    work = _empic_collection(es)
    _require_finite_real("background_charge", es.n0, T)
    for (i, ps) in enumerate(species)
        _validated_boris_factor(ps)
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "electromagnetic collection")
        _empic_collection_validate_polarization(es, ps, work.mobile[i])
        relativistic && _require_subluminal_velocities(ps, c, "electromagnetic collection")
        all(w -> (isfinite(w) && _stored_nonnegative(w))::Bool, ps.weight) ||
            throw(ArgumentError("electromagnetic weights must be finite and nonnegative"))
        !work.mobile[i] &&
            !all(a -> all(_particle_stored_zero, a)::Bool, ps.v) &&
            throw(
                ArgumentError(
                    "fixed electromagnetic populations require original stored zero velocity",
                ),
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
    schedule = _empic_collection(es).schedule
    schedule isa _EMPICEventSchedule &&
        es.n_sub > schedule.proposal.budget &&
        throw(ArgumentError("base field partition exceeds persistent event budget"))
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
    return _empic_charge_sum!(rho, stage.ne, work.staged, stage.g, stage.shape, stage.n0, target)
end

function _empic_charge_sum!(rho, density, species, g, shape, background, target)
    _charge_density_sum!(rho, density, species, g, shape; background = background)
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
    work.schedule isa _EMPICCoupledSchedule && return _empic_coupled_init!(es, species)
    all(ps -> ps isa eltype(work.owners), species) || throw(
        ArgumentError(
            "population storage type changed; construct a new electromagnetic collection state",
        ),
    )
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
    _empic_collection_finish_candidate!(es, work, species)
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
    return _empic_current_sum!(
        stage.Jraw,
        stage,
        work.staged,
        work.midpoints,
        work.mobile,
        stage.g,
        stage.shape,
        Val(:recorded_midpoints),
    )
end

_empic_current_prepare!(::Val{:recorded_midpoints}, stage, ps) =
    _ensure_empic_electron_buffers!(stage, nparticles(ps))
_empic_current_prepare!(::Val{:stored_positions}, data, ps) = resize!(data.work, nparticles(ps))

function _empic_current_population!(::Val{:recorded_midpoints}, Jraw, stage, ps, mid, g, shape)
    _swap_positions_with_midpoints!(ps, mid)
    return try
        _add_species_current!(Jraw, stage, ps, stage.worke, stage.ne, true)
    finally
        _swap_positions_with_midpoints!(ps, mid)
    end
end
function _empic_current_population!(::Val{:stored_positions}, Jraw, data, ps, mid, g, shape)
    return _empic_add_species_current!(Jraw, ps, data.work, data.density, g, shape, true)
end
_empic_current_exact(::Val{:recorded_midpoints}, ps, mid, g, shape, c) =
    _exact_midpoint_current(ps, mid, g, shape, c)
_empic_current_exact(::Val{:stored_positions}, ps, mid, g, shape, c) =
    _exact_moment(ps, (ps.weight, ps.v[c]), g, shape, ps.q)

function _empic_current_sum!(
    Jraw::NTuple{3,Array{T,D}},
    data,
    species,
    midpoints,
    mobile,
    g,
    shape,
    placement,
) where {D,T}
    foreach(a -> fill!(a, zero(T)), Jraw)
    eligible, lost = true, false
    for (i, ps) in enumerate(species)
        mobile[i] || continue
        _empic_current_prepare!(placement, data, ps)
        states = _empic_current_population!(placement, Jraw, data, ps, midpoints[i], g, shape)
        eligible &= all(s -> s.eligible, states)
        lost |= any(s -> s.lost, states)
        lost |= any(a -> any(_deposition_range_lost, a), Jraw)
    end
    if eligible && lost
        for c = 1:3
            _empic_exact_current_sum!(Jraw[c], species, midpoints, mobile, g, shape, c, placement)
        end
    end
    all(a -> all(isfinite, a), Jraw) ||
        throw(ArgumentError("electromagnetic deposited current exceeds numeric range"))
    return Jraw
end

function _empic_collection_exact_current!(out::Array{T}, stage, work, c) where {T}
    return _empic_exact_current_sum!(
        out,
        work.staged,
        work.midpoints,
        work.mobile,
        stage.g,
        stage.shape,
        c,
        Val(:recorded_midpoints),
    )
end
function _empic_exact_current_sum!(
    out::Array{T},
    species,
    midpoints,
    mobile,
    g,
    shape,
    c,
    placement,
) where {T}
    exact = fill(zero(Rational{BigInt}), g.n)
    for (i, ps) in enumerate(species)
        mobile[i] || continue
        contribution = _empic_current_exact(placement, ps, midpoints[i], g, shape, c)
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

function _empic_collection_substep!(
    stage::EMPIC{D,T},
    work,
    dt::T,
    previous::T,
    target,
    index,
) where {D,T}
    parts = _empic_collection_begin_interval!(stage, work, target, previous)
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        if _empic_collection_kick_due(work.schedule, i, index)
            _empic_collection_gather!(stage, ps)
            interval = _empic_collection_kick_interval(work.schedule, i, previous, dt)
            _empic_collection_push_gathered!(stage, ps, interval)
            _empic_collection_record_kick!(work.schedule, i)
        end
        _empic_collection_drift!(stage, work, i, dt)
    end
    return _empic_collection_finish_interval!(stage, work, parts, previous, dt, target)
end

function _empic_collection_begin_interval!(stage::EMPIC, work, target, previous)
    _empic_collection_charge!(stage.rho_n, stage, work, target)
    return _empic_center_magnetic!(stage, stage.g, previous)
end
function _empic_collection_gather!(
    stage::EMPIC{D,T},
    ps::_ElectrostaticHostParticle{D,T};
    initial = false,
) where {D,T}
    _ensure_empic_electron_buffers!(stage, nparticles(ps))
    gather_vector!(stage.Ep, stage.E, ps, stage.g, stage.shape)
    gather_vector!(stage.Bp, initial ? stage.B : stage.Bc, ps, stage.g, stage.shape)
    return nothing
end
_empic_collection_push_gathered!(stage::EMPIC, ps, interval) =
    _push_velocities!(stage, ps, stage.Ep, stage.Bp, interval)
function _empic_collection_drift!(stage::EMPIC{D,T}, work, i, dt) where {D,T}
    ps = work.staged[i]
    _drift_record_mid!(ps, work.midpoints[i], dt, stage.g.L)
    apply_periodic!(ps, ntuple(_ -> zero(T), Val(D)), stage.g.L)
    _require_finite_particle_velocities(ps, "electromagnetic collection output")
    stage.relativistic &&
        _require_subluminal_velocities(ps, stage.c, "electromagnetic collection output")
    return nothing
end
function _empic_collection_finish_interval!(stage::EMPIC, work, parts, previous, dt, target)
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
    _validate_empic_clock_phase(es.time[], es.step[], es.last_dt[], es.last_sub_dt[])
    _empic_collection_validate_species_history(es, _empic_collection(es).schedule)
    return nothing
end

function step_empic!(es::Union{EMPIC,EMPIC1D}, species::AbstractVector, dt::Real)
    es.collection !== nothing &&
        es.collection.schedule isa _EMPICCoupledSchedule &&
        return _empic_coupled_step!(es, species, dt)
    T = eltype(es.rho_n)
    dtT =
        _require_resolved_dt("step_empic!", _validated_nonnegative_dt(T, dt; name = "step_empic!"))
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    if work.schedule isa _EMPICEventSchedule
        isfinite(dt) && dt >= 0 && (iszero(dt) || !iszero(dtT)) || throw(
            ArgumentError(
                "persistent duration must be nonnegative and nonzero inputs must survive conversion",
            ),
        )
    end
    work.initialized[] ||
        throw(ArgumentError("initialize the electromagnetic collection before stepping"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    iszero(dtT) && return es
    if work.schedule isa _EMPICEventSchedule
        _empic_event_advance!(es, species, work, charges, dtT)
        return _empic_collection_publish!(es, work; particles = true)
    end
    sub_dt = _require_resolved_dt("step_empic! substep dt/n_sub", dtT / es.n_sub)
    sub_dt > zero(T) || throw(ArgumentError("dt/n_sub must be representable and positive"))
    next_time, next_step = _validated_step_clock(es.time[], es.step[], sub_dt, es.n_sub)
    _empic_collection_advance_candidate!(
        es,
        species,
        work,
        charges,
        dtT,
        sub_dt,
        next_time,
        next_step,
    )
    return _empic_collection_publish!(es, work; particles = true)
end

function _empic_collection_advance_candidate!(
    es,
    species,
    work,
    charges,
    dtT,
    sub_dt,
    next_time,
    next_step,
)
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
    es.collection !== nothing &&
        es.collection.schedule isa _EMPICCoupledSchedule &&
        return _empic_coupled_diagnostics(es, species)
    T = eltype(es.rho_n)
    charges = _empic_collection_validate_species(es, species)
    work = _empic_collection(es)
    work.initialized[] ||
        throw(ArgumentError("initialize the electromagnetic collection before diagnostics"))
    _empic_collection_bound(work, species)
    _empic_collection_validate_phase(es)
    _empic_collection_require_finite_fields(es)
    return _empic_collection_moment_diagnostics(es, species, charges, work)
end
function _empic_collection_moment_diagnostics(es, species, charges, work)
    T = eltype(es.rho_n)
    kinetic, momentum =
        es.relativistic ? _relativistic_particle_moments(species, es.c) :
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
        velocity_time_offsets = _empic_collection_velocity_offsets(
            work.schedule,
            work.mobile,
            offset,
        ),
        magnetic_time_offset = offset,
    )
    return _empic_collection_schedule_diagnostics(work.schedule, result)
end

"""
    EMPIC1D(g, species::AbstractVector; mobile=trues(length(species)),
            background_charge=0, neutrality=:require, charge_rtol=64eps(T),
            c=5, shape=CIC(), n_sub=1, species_stride=nothing,
            species_dt=nothing, event_budget=nothing, error_control=nothing,
            relativistic=false)

Specialized periodic CPU Float32/64 PIC with fields `Ex/Ey/Bz` and ordered
arbitrary finite species charges and positive masses. Every original stored
`v_z` must be zero: the omitted `Ez/By` polarization is not evolved. Fixed
populations require zero velocity and deposit charge only. Every mobile
population drifts and contributes current on every actual field segment; its integer
`species_stride` divides `n_sub` and selects kicks at fast indices `0,r,2r,...`.
Default mobile strides are 1 and fixed strides are 0. The species kick interval
is `H=dt/(n_sub÷r)`, with its own previous-half history. Use
[`set_empic1d_schedule!`](@ref) between completed external operations. NGP, CIC
and TSC are supported.
Alternatively, `species_dt` plus `event_budget` opts into the persistent
nondividing schedule described by [`EMPIC`](@ref), including its fixed-source
policy, optional sampled-state adaptation and schema3/policy2 restart.
Adding the closed `error_control` proposal selects full-time CIC/TSC projected
Euler/Heun evolution, complete signed work and schema4/policy3 restart, as
specified in `docs/EMPIC_COUPLED_ERROR_CONTROL.md`. This constructor selection
keeps the planar `Ex/Ey/Bz` polarization and explicitly rejects NGP.

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
    species_dt = nothing,
    species_adaptive = nothing,
    species_dt_min = nothing,
    schedule_limits = nothing,
    event_budget = nothing,
    error_control = nothing,
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
    error_control === nothing ||
        shape isa Union{CIC,TSC} ||
        throw(
            ArgumentError(
                "coupled error control requires CIC or TSC; NGP jump control is not implemented",
            ),
        )
    neutrality in (:require, :neutralize) ||
        throw(ArgumentError("neutrality must be :require or :neutralize"))
    _espic_validate_collection_type(species, g)
    axes(mobile, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("mobile must have one Boolean per species with one-based axes"))
    n_sub >= 1 || throw(ArgumentError("n_sub must be ≥ 1"))
    species_dt === nothing ||
        (!(n_sub isa Bool) && n_sub <= typemax(Int)) ||
        throw(ArgumentError("persistent n_sub must be a positive fitting Int"))
    ns = Int(n_sub)
    background = _require_finite_real("background_charge", background_charge, T)
    rtol = _require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
    cT = _require_finite_positive_real("c", c, T)
    stage =
        _new_empic1d_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, nothing)
    P = eltype(species)
    rows =
        species_dt === nothing && species_stride === nothing ?
        [active ? 1 : 0 for active in mobile] : species_stride
    schedule = _empic_collection_schedule(
        mobile,
        ns,
        rows,
        species_dt,
        species_adaptive,
        species_dt_min,
        schedule_limits,
        event_budget,
        T,
        error_control,
    )
    work = _EMPICCollectionWork{1,T,P,typeof(stage),typeof(schedule)}(
        neutrality,
        rtol,
        BitVector(mobile),
        P[],
        _ElectrostaticHostParticle{1,T}[],
        NTuple{1,Vector{T}}[],
        stage,
        Ref(false),
        schedule,
    )
    es = _new_empic1d_state(g, 0, background, cT, shape, relativistic, false, one(T), ns, work)
    _empic_coupled_bind!(es, schedule)
    _empic_collection_validate_species(es, species)
    return es
end

_empic_collection_validate_polarization(::EMPIC, ps, mobile) = nothing
function _empic_collection_validate_polarization(::EMPIC1D, ps, mobile)
    mobile &&
        !all(_particle_stored_zero, ps.v[3]) &&
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
        _prime_empic1d!(
            ps.v,
            stage.Exp,
            stage.Eyp,
            stage.Bzp,
            _validated_boris_factor(ps),
            _empic_collection_prime_interval(work.schedule, i, sub_dt) / 2,
            nparticles(ps),
            stage.relativistic,
            stage.c,
        )
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

function _empic_collection_substep!(
    stage::EMPIC1D{T},
    work,
    dt::T,
    previous::T,
    target,
    index,
) where {T}
    g = stage.g
    parts = _empic_collection_begin_interval!(stage, work, target, previous)
    for (i, ps) in enumerate(work.staged)
        work.mobile[i] || continue
        if _empic_collection_kick_due(work.schedule, i, index)
            _empic_collection_gather!(stage, ps)
            interval = _empic_collection_kick_interval(work.schedule, i, previous, dt)
            _empic_collection_push_gathered!(stage, ps, interval)
            _empic_collection_record_kick!(work.schedule, i)
        end
        _empic_collection_drift!(stage, work, i, dt)
    end
    return _empic_collection_finish_interval!(stage, work, parts, previous, dt, target)
end

function _empic_collection_begin_interval!(stage::EMPIC1D{T}, work, target, previous) where {T}
    _empic_collection_charge!(stage.rho_n, stage, work, target)
    return _empic1d_center_magnetic!(stage, stage.g, previous, zero(T))
end
function _empic_collection_gather!(
    stage::EMPIC1D{T},
    ps::_ElectrostaticHostParticle{1,T};
    initial = false,
) where {T}
    _ensure_electron_buffers!(stage, nparticles(ps))
    gather_scalar!(stage.Exp, stage.Ex, ps, stage.g, stage.shape)
    gather_scalar!(stage.Eyp, stage.Ey, ps, stage.g, stage.shape)
    gather_scalar!(stage.Bzp, initial ? stage.Bz : stage.Bzc, ps, stage.g, stage.shape)
    return nothing
end
_empic_collection_push_gathered!(stage::EMPIC1D, ps, interval) =
    _push_velocities!(stage, ps, stage.Exp, stage.Eyp, stage.Bzp, interval)
function _empic_collection_drift!(stage::EMPIC1D{T}, work, i, dt) where {T}
    ps = work.staged[i]
    _drift_record_mid!(ps, work.midpoints[i][1], dt, stage.g.L[1])
    apply_periodic!(ps, (zero(T),), stage.g.L)
    _validate_particle_positions(ps)
    _require_finite_particle_velocities(ps, "electromagnetic collection output")
    stage.relativistic &&
        _require_subluminal_velocities(ps, stage.c, "electromagnetic collection output")
    return nothing
end
function _empic_collection_finish_interval!(
    stage::EMPIC1D{T},
    work,
    parts,
    previous,
    dt,
    target,
) where {T}
    g = stage.g
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
function _empic_collection_schedule_axes(schedule, a, count)
    axes(a, 1) == Base.OneTo(count) ||
        throw(ArgumentError("electromagnetic schedule storage must have one row per population"))
    return nothing
end
_empic_collection_finish_candidate!(es, work, species) = nothing
_empic_collection_schedule_arrays(schedule::_EMPICSpeciesSchedule) = (
    schedule.stride,
    schedule.last_species_dt,
    schedule.new_species_dt,
    schedule.previous_species_dt,
)

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
    return _EMPICSpeciesSchedule{T}(
        stride,
        zeros(T, length(mobile)),
        zeros(T, length(mobile)),
        zeros(T, length(mobile)),
    )
end

_empic_collection_prime_interval(::Nothing, i, fast) = fast
_empic_collection_prime_interval(schedule::_EMPICSpeciesSchedule, i, fast) =
    schedule.new_species_dt[i]
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
            isfinite(h) && h >= 2 * floatmin(T) && fast <= h <= full || throw(
                ArgumentError(
                    "species timestep history must lie between previous fast and full intervals",
                ),
            )
        end
    end
    return nothing
end

_empic_collection_validate_schedule(es, ::Nothing) = nothing
_empic_collection_validate_schedule(es, schedule::_EMPICSpeciesSchedule) =
    _empic_validate_species_schedule(es.collection.mobile, es.n_sub, schedule.stride)
_empic_collection_validate_species_history(es, ::Nothing) = nothing
_empic_collection_validate_species_history(es, schedule::_EMPICSpeciesSchedule{T}) where {T} =
    _empic_validate_species_history(
        es.collection.mobile,
        schedule.last_species_dt,
        es.step[],
        es.last_dt[],
        es.last_sub_dt[],
        T,
    )

_empic_collection_prepare_intervals!(::Nothing, mobile, es, dt, fast) = nothing
function _empic_collection_prepare_intervals!(
    schedule::_EMPICSpeciesSchedule{T},
    mobile,
    es,
    dt,
    fast,
) where {T}
    for (i, active) in enumerate(mobile)
        h = active ? dt / (es.n_sub ÷ schedule.stride[i]) : zero(T)
        if active
            isfinite(h) && fast <= h <= dt && h > zero(T) || throw(
                ArgumentError("species drift interval must lie between fast and full intervals"),
            )
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
_empic_collection_velocity_offsets(schedule::_EMPICSpeciesSchedule{T}, mobile, offset) where {T} = [
    active && !iszero(h) ? -h / 2 : zero(T) for (active, h) in zip(mobile, schedule.last_species_dt)
]
_empic_collection_schedule_diagnostics(::Nothing, result) = result
_empic_collection_schedule_diagnostics(schedule::_EMPICSpeciesSchedule, result) = merge(
    result,
    (species_stride = copy(schedule.stride), last_species_dt = copy(schedule.last_species_dt)),
)

"""
    set_empic1d_schedule!(state; n_sub=state.n_sub, species_stride=nothing,
                         species_dt=nothing, species_adaptive=nothing,
                         species_dt_min=nothing, schedule_limits=nothing,
                         event_budget=nothing)

Set a specialized collection's future dividing integer kick strides at a completed
external boundary. `nothing` keeps its current strides. Mobile rows divide
`n_sub`; fixed rows are zero. The complete proposal is copied and validated
before configuration changes. Previous species and magnetic half intervals are
preserved; the next positive step joins the old and new species half levels.
The default all-fast state can be changed before or after initialization.
For a state constructed with `species_dt`, the persistent keywords instead
replace its future proposal; `nothing` retains each current value. Pending
intervals/deadlines remain intact and new choices apply at the next due kick.
Dividing and persistent owner kinds cannot be switched by this setter.
On a constructor-selected full-time controller, `error_control` replaces its
complete future proposal. Historical admission, pending deadlines and work stay
owned by their accepted policy; no setter migrates a half phase.
"""
set_empic1d_schedule!(
    es::EMPIC1D;
    n_sub::Integer = es.n_sub,
    species_stride = nothing,
    species_dt = nothing,
    species_adaptive = nothing,
    species_dt_min = nothing,
    schedule_limits = nothing,
    event_budget = nothing,
    error_control = nothing,
) = _set_empic_collection_schedule!(
    es,
    n_sub,
    species_stride,
    species_dt,
    species_adaptive,
    species_dt_min,
    schedule_limits,
    event_budget,
    error_control,
)

"""
    set_empic_schedule!(state::EMPIC; n_sub=state.n_sub, species_stride=nothing,
                       species_dt=nothing, species_adaptive=nothing,
                       species_dt_min=nothing, schedule_limits=nothing,
                       event_budget=nothing)

Set future dividing species kick strides on an opted-in generic collection.
Construct with explicit `species_stride` to opt in; a default synchronous state
keeps its concrete Nothing schedule owner. The complete proposal is copied and
validated before either configuration value changes. Prior species/fast histories
are preserved, and `nothing` keeps current strides. Call only at completed
external boundaries with no concurrent use.
For an opted-in persistent owner, copied `species_dt`, `species_adaptive`,
`species_dt_min`, `schedule_limits` and `event_budget` replace future choices
without retiming the current pending intervals. `nothing` retains each value.
For a constructor-selected full-time controller, a complete `error_control`
proposal updates future choices while retaining historical admission and work.
"""
set_empic_schedule!(
    es::EMPIC;
    n_sub::Integer = es.n_sub,
    species_stride = nothing,
    species_dt = nothing,
    species_adaptive = nothing,
    species_dt_min = nothing,
    schedule_limits = nothing,
    event_budget = nothing,
    error_control = nothing,
) = _set_empic_collection_schedule!(
    es,
    n_sub,
    species_stride,
    species_dt,
    species_adaptive,
    species_dt_min,
    schedule_limits,
    event_budget,
    error_control,
)

function _set_empic_collection_schedule!(
    es,
    n_sub,
    stride,
    cap,
    adaptive,
    minimum,
    limits,
    budget,
    error_control = nothing,
)
    work = _empic_collection(es)
    work.schedule isa _EMPICCoupledSchedule && return _empic_coupled_set!(
        es,
        n_sub,
        stride,
        cap,
        adaptive,
        minimum,
        limits,
        budget,
        error_control,
    )
    error_control === nothing || throw(
        ArgumentError(
            "construct a coupled error-control owner before changing its control proposal",
        ),
    )
    if work.schedule isa _EMPICEventSchedule
        stride === nothing || throw(ArgumentError("persistent schedules do not use species_stride"))
        n_sub isa Integer && !(n_sub isa Bool) && 1 <= n_sub <= typemax(Int) ||
            throw(ArgumentError("n_sub must be a positive Int for persistent schedules"))
        old = work.schedule.proposal
        proposed = _empic_event_proposal(
            work.mobile,
            cap === nothing ? old.cap : cap,
            adaptive === nothing ? old.adaptive : adaptive,
            minimum === nothing ? old.minimum : minimum,
            limits === nothing ? old.limits : limits,
            budget === nothing ? old.budget : budget,
            eltype(es.rho_n),
        )
        n_sub <= proposed.budget ||
            throw(ArgumentError("base field partition exceeds persistent event budget"))
        _empic_collection_validate_schedule(es, work.schedule)
        _empic_collection_validate_phase(es)
        work.initialized[] && _empic_collection_validate_species(es, work.owners)
        es.n_sub = Int(n_sub)
        work.schedule.proposal = proposed
        return es
    end
    all(isnothing, (cap, adaptive, minimum, limits, budget)) ||
        throw(ArgumentError("persistent schedule changes require construction with species_dt"))
    return _set_empic_collection_schedule!(es, n_sub, stride)
end

function _set_empic_collection_schedule!(es, n_sub, species_stride)
    work = _empic_collection(es)
    work.schedule isa _EMPICSpeciesSchedule || throw(
        ArgumentError(
            "construct EMPIC with explicit species_stride before setting a species schedule",
        ),
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

# Persistent event phases use the same numerical interval owners above. Their
# rational source clock orders events; actual narrowed durations remain visible.
_empic_event_q(value) = value isa _EMPICEventQ && denominator(value) > 0
_empic_event_arrays(ps::ParticleSet) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
_empic_event_E(es::EMPIC) = es.E
_empic_event_E(es::EMPIC1D) = (es.Ex, es.Ey)
_empic_event_B(es::EMPIC) = es.B
_empic_event_B(es::EMPIC1D) = (es.Bz,)
_empic_event_J(es::EMPIC) = es.J
_empic_event_J(es::EMPIC1D) = (es.Jx, es.Jy)
_empic_event_raw(es::EMPIC) = es.Jraw
_empic_event_raw(es::EMPIC1D) = nothing
_empic_event_fields(es::EMPIC) = merge(_empic_physical_fields(es), (Jraw = es.Jraw,))
_empic_event_fields(es::EMPIC1D) = _empic_physical_fields(es)
_empic_event_energy(es) =
    _field_energy_rational(_empic_event_E(es), _empic_event_B(es), _exact_rational(es.c)^2, es.g.dx)
_empic_event_field_E(::EMPIC, fields) = fields.E
_empic_event_field_E(::EMPIC1D, fields) = (fields.Ex, fields.Ey)
_empic_event_field_B(::EMPIC, fields) = fields.B
_empic_event_field_B(::EMPIC1D, fields) = (fields.Bz,)
_empic_event_field_J(::EMPIC, fields) = fields.J
_empic_event_field_J(::EMPIC1D, fields) = (fields.Jx, fields.Jy)

function _empic_event_floor(::Type{T}, exact) where {T}
    value = _from_exact(T, exact)
    isfinite(value) || throw(ArgumentError("persistent schedule duration exceeds numeric range"))
    _exact_rational(value) > exact && (value = prevfloat(value))
    value > zero(T) || throw(ArgumentError("persistent schedule duration disappears"))
    return _require_resolved_dt("persistent species interval", value)
end
function _empic_event_interval(::Type{T}, exact) where {T}
    value = _from_exact(T, exact)
    isfinite(value) && value > zero(T) ||
        throw(ArgumentError("distinct persistent events require a representable positive interval"))
    return _require_resolved_dt("persistent field interval", value)
end
function _empic_event_bounds(stage::EMPIC, ps)
    Q = _exact_rational
    all(a -> all(isfinite, a), (stage.Ep..., stage.Bp...)) ||
        throw(ArgumentError("persistent sampled fields must be finite"))
    E = B = V = zero(_EMPICEventQ)
    for p in eachindex(ps.weight)
        E = max(E, sum(Q(abs(stage.Ep[c][p])) for c = 1:3))
        B = max(B, sum(Q(abs(stage.Bp[c][p])) for c = 1:3))
        V = max(V, sum(Q(abs(ps.v[c][p])) for c = 1:3))
    end
    return E, B, V
end
function _empic_event_bounds(stage::EMPIC1D, ps)
    Q = _exact_rational
    all(a -> all(isfinite, a), (stage.Exp, stage.Eyp, stage.Bzp)) ||
        throw(ArgumentError("persistent sampled fields must be finite"))
    E = B = V = zero(_EMPICEventQ)
    for p in eachindex(ps.weight)
        E = max(E, Q(abs(stage.Exp[p])) + Q(abs(stage.Eyp[p])))
        B = max(B, Q(abs(stage.Bzp[p])))
        V = max(V, sum(Q(abs(ps.v[c][p])) for c = 1:3))
    end
    return E, B, V
end
const _EMPIC_EVENT_KICK_KEYS = (
    :species,
    :time,
    :previous,
    :selected,
    :cap,
    :upper,
    :minimum,
    :adaptive,
    :limits,
    :q,
    :m,
    :E1,
    :B1,
    :V1,
    :halvings,
    :bootstrap,
    :duration,
    :duration_defect,
)
function _empic_event_within(h, q, m, E, B, V, limits, c, dx)
    Q = _exact_rational
    return Q(h) * abs(Q(q)) * B <= Q(limits.rotation) * Q(m) &&
           Q(h) * abs(Q(q)) * E <= Q(limits.electric) * Q(m) * Q(c) &&
           Q(h) * V <= Q(limits.travel) * minimum(Q, dx)
end
function _empic_event_interval_choice(stage, proposal, ps, i, previous, bootstrap)
    T, Q = eltype(ps), _exact_rational
    E, B, V = _empic_event_bounds(stage, ps)
    adaptive = proposal.adaptive[i]
    cap = Q(proposal.cap[i])
    adaptive && !bootstrap && (cap = min(cap, Q(proposal.limits.growth) * Q(previous)))
    selected, halvings = _empic_event_floor(T, cap), 0
    if adaptive
        lower = proposal.minimum[i]
        while true
            selected >= lower ||
                throw(ArgumentError("adaptive species interval fell below species_dt_min"))
            _empic_event_within(
                selected,
                ps.q,
                ps.m,
                E,
                B,
                V,
                proposal.limits,
                stage.c,
                stage.g.dx,
            ) && break
            halvings < proposal.budget ||
                throw(ArgumentError("adaptive interval halving budget exhausted"))
            selected = _require_resolved_dt("adaptive species interval", selected / T(2))
            halvings += 1
        end
    end
    return (
        selected = selected,
        cap = cap,
        upper = proposal.cap[i],
        minimum = proposal.minimum[i],
        adaptive = adaptive,
        limits = proposal.limits,
        q = ps.q,
        m = ps.m,
        E1 = E,
        B1 = B,
        V1 = V,
        halvings = halvings,
    )
end

function _empic_event_decision(stage, schedule, ps, i, time, previous, bootstrap)
    T, Q, proposal = eltype(ps), _exact_rational, schedule.proposal
    choice = _empic_event_interval_choice(stage, proposal, ps, i, previous, bootstrap)
    selected, cap, adaptive, E, B, V, halvings = choice.selected,
    choice.cap,
    choice.adaptive,
    choice.E1,
    choice.B1,
    choice.V1,
    choice.halvings
    prior = bootstrap ? selected : previous
    consumed = prior / T(2) + selected / T(2)
    return (
        species = i,
        time = time,
        previous = prior,
        selected = selected,
        cap = cap,
        upper = proposal.cap[i],
        minimum = proposal.minimum[i],
        adaptive = adaptive,
        limits = proposal.limits,
        q = ps.q,
        m = ps.m,
        E1 = E,
        B1 = B,
        V1 = V,
        halvings = halvings,
        bootstrap = bootstrap,
        duration = consumed,
        duration_defect = Q(consumed) - (Q(prior) + Q(selected)) / 2,
    )
end

function _empic_event_configuration(es; c = es.c, shape = es.shape, relativistic = es.relativistic)
    work = _empic_collection(es)
    return (
        c = c,
        shape = shape,
        relativistic = relativistic,
        background = es.n0,
        neutrality = work.neutrality,
        charge_rtol = work.charge_rtol,
        mobile = Tuple(work.mobile),
    )
end
function _empic_event_shadow_values(es, populations, fields, bindings, configuration, clock)
    D, T = length(es.g.n), eltype(es.rho_n)
    copies = _ElectrostaticHostParticle{D,T}[]
    _pic_prepare_species_stages!(copies, populations)
    detached = NamedTuple{keys(fields)}(map(values(fields)) do value
        value isa Tuple ? map(copy, value) : copy(value)
    end)
    return (
        species = copies,
        fields = detached,
        bindings = Tuple[_empic_event_arrays(ps) for ps in bindings],
        configuration = configuration,
        clock = clock,
    )
end
function _empic_event_shadow(es, work, species)
    stage = work.stage
    return _empic_event_shadow_values(
        es,
        work.staged,
        _empic_event_fields(stage),
        species,
        _empic_event_configuration(es),
        (stage.time[], stage.step[], stage.last_dt[], stage.last_sub_dt[]),
    )
end
function _empic_event_shadow_matches(es, schedule)
    shadow, work = schedule.shadow, _empic_collection(es)
    shadow isa NamedTuple &&
        keys(shadow) == (:species, :fields, :bindings, :configuration, :clock) ||
        throw(ArgumentError("persistent schedule has no valid committed shadow"))
    expected = _empic_event_configuration(es)
    isequal(shadow.configuration, expected) &&
        isequal(shadow.clock, (es.time[], es.step[], es.last_dt[], es.last_sub_dt[])) &&
        isequal(shadow.fields, _empic_event_fields(es)) &&
        length(shadow.species) == length(shadow.bindings) == length(work.owners) || throw(
        ArgumentError("persistent source/field/clock changed; initialize from full-time data"),
    )
    for (ps, saved, binding) in zip(work.owners, shadow.species, shadow.bindings)
        length(binding) == length(_empic_event_arrays(ps)) &&
            all(pair -> pair[1] === pair[2], zip(_empic_event_arrays(ps), binding)) &&
            isequal((ps.q, ps.m), (saved.q, saved.m)) &&
            all(
                pair -> isequal(pair[1], pair[2]),
                zip(_empic_event_arrays(ps), _empic_event_arrays(saved)),
            ) || throw(
            ArgumentError(
                "persistent particle source/phase/provenance changed; initialize from full-time data",
            ),
        )
    end
    return nothing
end
function _empic_event_shadow_arrays(shadow)
    shadow === nothing && return AbstractArray[]
    shadow isa NamedTuple &&
        keys(shadow) == (:species, :fields, :bindings, :configuration, :clock) ||
        throw(ArgumentError("persistent committed shadow layout changed"))
    # Bindings intentionally retain public component identities; their data is
    # not a second private storage role. The copied physical values are separate.
    arrays = AbstractArray[shadow.species, shadow.bindings]
    for ps in shadow.species
        append!(arrays, _empic_event_arrays(ps))
    end
    for value in values(shadow.fields)
        append!(arrays, value isa Tuple ? value : (value,))
    end
    return arrays
end
function _empic_event_report_arrays(report)
    report === nothing && return AbstractArray[]
    arrays = AbstractArray[
        report.initial.interval,
        report.initial.last,
        report.initial.next,
        report.intervals,
        report.edges,
        report.duration_defects,
        report.kicks,
        report.continuity,
        report.segments,
        report.rho_start,
        report.rho_end,
        report.impulse...,
    ]
    report.raw_impulse === nothing || append!(arrays, report.raw_impulse)
    for segment in report.segments
        append!(
            arrays,
            (
                segment.E_before...,
                segment.E_after...,
                segment.current...,
                segment.rho_before,
                segment.rho_after,
            ),
        )
        segment.raw_current === nothing || append!(arrays, segment.raw_current)
    end
    return arrays
end
function _empic_collection_schedule_arrays(schedule::_EMPICEventSchedule)
    p, a, b = schedule.proposal, schedule.phase, schedule.candidate
    arrays = AbstractArray[
        p.cap,
        p.adaptive,
        p.minimum,
        a.interval,
        a.last,
        a.next,
        b.interval,
        b.last,
        b.next,
    ]
    append!(arrays, _empic_event_report_arrays(a.report))
    append!(arrays, _empic_event_report_arrays(b.report))
    append!(arrays, _empic_event_shadow_arrays(schedule.shadow))
    append!(arrays, _empic_event_shadow_arrays(schedule.candidate_shadow))
    return arrays
end
_empic_collection_schedule_axes(::_EMPICEventSchedule, a, count) = nothing
function _empic_collection_validate_schedule(es, schedule::_EMPICEventSchedule{T}) where {T}
    rounding(T) === RoundNearest && !get_zero_subnormals() ||
        throw(ArgumentError("persistent schedules require nearest rounding and gradual underflow"))
    p = schedule.proposal
    p isa NamedTuple && keys(p) == (:cap, :adaptive, :minimum, :limits, :budget) ||
        throw(ArgumentError("persistent proposal layout changed"))
    _empic_event_proposal(es.collection.mobile, p.cap, p.adaptive, p.minimum, p.limits, p.budget, T)
    for phase in (schedule.phase, schedule.candidate)
        all(
            a -> axes(a, 1) == Base.OneTo(length(es.collection.mobile)),
            (phase.interval, phase.last, phase.next),
        ) || throw(ArgumentError("persistent phase row geometry changed"))
        phase.report === nothing || _empic_event_report_layout(es, phase.report)
    end
    return nothing
end

const _EMPIC_EVENT_REPORT_KEYS = (
    :kind,
    :start,
    :finish,
    :requested,
    :public_before,
    :public_after,
    :budget,
    :n_sub,
    :bootstrap,
    :initial,
    :intervals,
    :edges,
    :duration_defects,
    :kicks,
    :continuity,
    :charge_source_scale,
    :segments,
    :rho_start,
    :rho_end,
    :impulse,
    :raw_impulse,
    :work_used,
    :work_raw,
    :field_before,
    :field_after,
    :field_defect,
    :clock_defect,
)
function _empic_event_report_layout(es, report)
    T, D = eltype(es.rho_n), length(es.g.n)
    report isa NamedTuple && keys(report) == _EMPIC_EVENT_REPORT_KEYS ||
        throw(ArgumentError("persistent report layout differs"))
    report.kind === (es isa EMPIC ? :vector : :planar) &&
        report.bootstrap isa Bool &&
        report.budget isa Int &&
        report.budget > 0 &&
        report.n_sub isa Int &&
        1 <= report.n_sub <= report.budget &&
        report.requested isa T &&
        report.public_before isa T &&
        report.public_after isa T &&
        all(isfinite, (report.requested, report.public_before, report.public_after)) &&
        report.requested > zero(T) &&
        report.public_before >= zero(T) &&
        report.public_after > report.public_before &&
        all(
            _empic_event_q,
            (
                report.start,
                report.finish,
                report.work_used,
                report.charge_source_scale,
                report.field_before,
                report.field_after,
                report.field_defect,
                report.clock_defect,
            ),
        ) || throw(ArgumentError("persistent report scalar source/clock values differ"))
    report.initial isa NamedTuple &&
        keys(report.initial) == (:primed, :interval, :last, :next) &&
        report.initial.primed isa Bool &&
        report.initial.interval isa Vector{T} &&
        report.initial.last isa Vector{_EMPICEventQ} &&
        report.initial.next isa Vector{_EMPICEventQ} &&
        all(
            a -> axes(a, 1) == Base.OneTo(length(es.collection.mobile)),
            (report.initial.interval, report.initial.last, report.initial.next),
        ) &&
        all(isfinite, report.initial.interval) &&
        all(_empic_event_q, (report.initial.last..., report.initial.next...)) ||
        throw(ArgumentError("persistent report initial phase differs"))
    report.intervals isa Vector{T} &&
        report.edges isa Vector{_EMPICEventQ} &&
        report.duration_defects isa Vector{_EMPICEventQ} &&
        report.continuity isa Vector{T} &&
        report.kicks isa Vector &&
        report.segments isa Vector &&
        all(v -> all(i -> isassigned(v, i), eachindex(v)), (report.kicks, report.segments)) ||
        throw(ArgumentError("persistent report event storage differs"))
    count = length(report.intervals)
    count > 0 &&
        length(report.edges) == count + 1 &&
        length(report.duration_defects) ==
        length(report.continuity) ==
        length(report.segments) ==
        count &&
        length(report.kicks) <= report.budget &&
        count <= report.budget - length(report.kicks) &&
        all(h -> isfinite(h) && h >= 2floatmin(T), report.intervals) &&
        all(x -> isfinite(x) && x >= zero(T), report.continuity) &&
        all(_empic_event_q, (report.edges..., report.duration_defects...)) ||
        throw(ArgumentError("persistent report event counts/ranges differ"))
    all(
        a -> a isa Array{T,D} && size(a) == es.g.n && all(isfinite, a),
        (report.rho_start, report.rho_end),
    ) || throw(ArgumentError("persistent report charge snapshots differ"))
    components = es isa EMPIC ? 3 : 2
    report.impulse isa NTuple{components,Array{_EMPICEventQ,D}} &&
        all(a -> size(a) == es.g.n && all(_empic_event_q, a), report.impulse) ||
        throw(ArgumentError("persistent report current impulse differs"))
    if es isa EMPIC
        report.raw_impulse isa NTuple{3,Array{_EMPICEventQ,D}} &&
            all(a -> size(a) == es.g.n && all(_empic_event_q, a), report.raw_impulse) &&
            _empic_event_q(report.work_raw) ||
            throw(ArgumentError("persistent report raw current/work differs"))
    else
        report.raw_impulse === nothing && report.work_raw === nothing ||
            throw(ArgumentError("planar reports cannot invent a raw longitudinal current"))
    end
    for segment in report.segments
        segment isa NamedTuple &&
            keys(segment) ==
            (:E_before, :E_after, :current, :raw_current, :rho_before, :rho_after) ||
            throw(ArgumentError("persistent segment source layout differs"))
        all(
            v ->
                v isa NTuple{components,Array{T,D}} &&
                    all(a -> size(a) == es.g.n && all(isfinite, a), v),
            (segment.E_before, segment.E_after, segment.current),
        ) || throw(ArgumentError("persistent segment electric/current sources differ"))
        all(
            a -> a isa Array{T,D} && size(a) == es.g.n && all(isfinite, a),
            (segment.rho_before, segment.rho_after),
        ) || throw(ArgumentError("persistent segment deposited charge differs"))
        if es isa EMPIC
            segment.raw_current isa NTuple{3,Array{T,D}} &&
                all(a -> size(a) == es.g.n && all(isfinite, a), segment.raw_current) ||
                throw(ArgumentError("persistent segment raw current differs"))
        else
            segment.raw_current === nothing ||
                throw(ArgumentError("planar segment cannot invent a raw longitudinal current"))
        end
    end
    for kick in report.kicks
        kick isa NamedTuple && keys(kick) == _EMPIC_EVENT_KICK_KEYS ||
            throw(ArgumentError("persistent kick decision layout differs"))
        kick.species isa Int &&
            1 <= kick.species <= length(es.collection.mobile) &&
            es.collection.mobile[kick.species] &&
            kick.bootstrap isa Bool &&
            kick.adaptive isa Bool &&
            kick.halvings isa Int &&
            0 <= kick.halvings <= report.budget &&
            all(
                x -> x isa T && isfinite(x),
                (
                    kick.previous,
                    kick.selected,
                    kick.upper,
                    kick.minimum,
                    kick.q,
                    kick.m,
                    kick.duration,
                ),
            ) &&
            kick.previous >= 2floatmin(T) &&
            kick.selected >= 2floatmin(T) &&
            kick.upper >= kick.selected &&
            kick.minimum >= zero(T) &&
            kick.m > zero(T) &&
            all(
                _empic_event_q,
                (kick.time, kick.cap, kick.E1, kick.B1, kick.V1, kick.duration_defect),
            ) &&
            all(>=(zero(_EMPICEventQ)), (kick.E1, kick.B1, kick.V1)) &&
            kick.limits isa NamedTuple &&
            keys(kick.limits) == (:rotation, :electric, :travel, :growth) &&
            all(x -> x isa T && isfinite(x), values(kick.limits)) ||
            throw(ArgumentError("persistent kick decision source values differ"))
    end
    return nothing
end

function _empic_event_report_valid(
    es,
    phase,
    species;
    speed = es.c,
    fields = _empic_event_fields(es),
)
    r, Q = phase.report, _exact_rational
    _empic_event_report_layout(es, r)
    r.start >= 0 &&
        r.finish == phase.time &&
        r.finish - r.start == Q(r.requested) &&
        first(r.edges) == r.start &&
        last(r.edges) == r.finish &&
        r.bootstrap === !r.initial.primed &&
        r.field_before >= 0 &&
        r.field_after >= 0 &&
        r.charge_source_scale == _pic_collection_charges(es.g, es.n0, species).scale &&
        r.field_defect == r.field_after - r.field_before + r.work_used &&
        r.clock_defect == Q(r.public_after) - Q(r.public_before) - (r.finish - r.start) ||
        throw(ArgumentError("persistent report source/work algebra differs"))
    clock = r.public_before
    for i in eachindex(r.intervals)
        delta = r.edges[i+1] - r.edges[i]
        delta > 0 &&
            isequal(r.intervals[i], _empic_event_interval(eltype(r.intervals), delta)) &&
            Q(r.intervals[i]) - delta == r.duration_defects[i] ||
            throw(ArgumentError("persistent report event duration/materialization differs"))
        nextclock, _ = _validated_step_clock(clock, 0, r.intervals[i])
        nextclock > clock || throw(ArgumentError("persistent report clock did not advance"))
        clock = nextclock
    end
    isequal(clock, r.public_after) || throw(ArgumentError("persistent report public clock differs"))
    isequal(last(r.segments).E_after, _empic_event_field_E(es, fields)) &&
        isequal(last(r.segments).current, _empic_event_field_J(es, fields)) &&
        isequal(first(r.segments).rho_before, r.rho_start) &&
        isequal(last(r.segments).rho_before, fields.rho_n) &&
        isequal(last(r.segments).rho_after, r.rho_end) &&
        isequal(r.rho_end, fields.rho_np1) &&
        r.field_after == _field_energy_rational(
            _empic_event_field_E(es, fields),
            _empic_event_field_B(es, fields),
            Q(speed)^2,
            es.g.dx,
        ) || throw(ArgumentError("persistent report differs from its actual saved endpoint"))
    es isa EMPIC &&
        !isequal(last(r.segments).raw_current, fields.Jraw) &&
        throw(ArgumentError("persistent report raw current differs from its actual saved endpoint"))
    impulse = map(a -> fill(zero(_EMPICEventQ), size(a)), r.impulse)
    raw_impulse =
        r.raw_impulse === nothing ? nothing :
        map(a -> fill(zero(_EMPICEventQ), size(a)), r.raw_impulse)
    work_used, work_raw = zero(_EMPICEventQ), r.work_raw === nothing ? nothing : zero(_EMPICEventQ)
    volume = prod(Q, es.g.dx)
    for i in eachindex(r.segments)
        segment = r.segments[i]
        i == 1 ||
            isequal(r.segments[i-1].E_after, segment.E_before) ||
            throw(ArgumentError("persistent segment electric history is discontinuous"))
        i == 1 ||
            isequal(r.segments[i-1].rho_after, segment.rho_before) ||
            throw(ArgumentError("persistent segment charge history is discontinuous"))
        work_used += _empic_event_work!(
            impulse,
            segment.current,
            segment.E_before,
            segment.E_after,
            r.intervals[i],
            volume,
        )
        if raw_impulse !== nothing
            work_raw += _empic_event_work!(
                raw_impulse,
                segment.raw_current,
                segment.E_before,
                segment.E_after,
                r.intervals[i],
                volume,
            )
        end
    end
    isequal(impulse, r.impulse) &&
        isequal(raw_impulse, r.raw_impulse) &&
        work_used == r.work_used &&
        isequal(work_raw, r.work_raw) ||
        throw(ArgumentError("persistent current impulse/work differs from actual segment sources"))
    H, a, b = copy(r.initial.interval), copy(r.initial.last), copy(r.initial.next)
    mobile = es.collection.mobile
    for i in eachindex(mobile)
        if !mobile[i] || r.bootstrap
            _particle_stored_zero(H[i]) && iszero(a[i]) && iszero(b[i]) ||
                throw(ArgumentError("persistent report inactive initial phase differs"))
        else
            H[i] >= 2floatmin(eltype(H)) && a[i] <= r.start <= b[i] && b[i] - a[i] == Q(H[i]) ||
                throw(ArgumentError("persistent report preceding active phase differs"))
        end
    end
    previous_time, previous_row = r.start, 0
    for k in r.kicks
        i = k.species
        r.start <= k.time < r.finish &&
            any(t -> t == k.time, r.edges[1:(end-1)]) &&
            (k.time > previous_time || (k.time == previous_time && i > previous_row)) &&
            k.q == species[i].q &&
            k.m == species[i].m ||
            throw(ArgumentError("persistent report kick ordering/source differs"))
        if r.bootstrap && k.time == r.start
            k.bootstrap && _particle_stored_zero(H[i]) && k.previous == k.selected ||
                throw(ArgumentError("persistent bootstrap decision was not reused"))
        else
            !k.bootstrap && k.time == b[i] && k.previous == H[i] ||
                throw(ArgumentError("persistent kick missed its pending deadline"))
        end
        cap =
            k.adaptive && !k.bootstrap ? min(Q(k.upper), Q(k.limits.growth) * Q(k.previous)) :
            Q(k.upper)
        cap == k.cap || throw(ArgumentError("persistent historical growth cap differs"))
        isequal(k.duration, k.previous / eltype(H)(2) + k.selected / eltype(H)(2)) &&
            k.duration_defect == Q(k.duration) - (Q(k.previous) + Q(k.selected)) / 2 ||
            throw(ArgumentError("persistent historical kick materialization differs"))
        expected = _empic_event_floor(eltype(H), cap)
        for _ = 1:k.halvings
            k.adaptive &&
                !_empic_event_within(
                    expected,
                    k.q,
                    k.m,
                    k.E1,
                    k.B1,
                    k.V1,
                    k.limits,
                    speed,
                    es.g.dx,
                ) ||
                throw(ArgumentError("persistent historical halving skipped an admitted interval"))
            expected = _require_resolved_dt("historical adaptive interval", expected / eltype(H)(2))
        end
        isequal(expected, k.selected) ||
            throw(ArgumentError("persistent historical halving differs"))
        if k.adaptive
            all(>(0), values(k.limits)) &&
                k.limits.growth >= 1 &&
                2floatmin(eltype(H)) <= k.minimum <= k.selected &&
                _empic_event_within(
                    k.selected,
                    k.q,
                    k.m,
                    k.E1,
                    k.B1,
                    k.V1,
                    k.limits,
                    speed,
                    es.g.dx,
                ) || throw(ArgumentError("persistent historical adaptive bounds differ"))
        else
            _particle_stored_zero(k.minimum) &&
                k.halvings == 0 &&
                (
                    (all(>(0), values(k.limits)) && k.limits.growth >= 1) || isequal(
                        k.limits,
                        (
                            rotation = zero(eltype(H)),
                            electric = zero(eltype(H)),
                            travel = zero(eltype(H)),
                            growth = one(eltype(H)),
                        ),
                    )
                ) || throw(ArgumentError("persistent fixed decision carries adaptive history"))
        end
        H[i], a[i], b[i] = k.selected, k.time, k.time + Q(k.selected)
        previous_time, previous_row = k.time, i
    end
    isequal(H, phase.interval) && a == phase.last && b == phase.next ||
        throw(ArgumentError("persistent pending phase differs from its actual decisions"))
    # Reuse the actual calendar selector to verify the whole historical field
    # partition, including base edges and due rows that cannot be omitted.
    replay = _EMPICEventPhase{eltype(H)}(
        r.start,
        r.initial.primed,
        copy(r.initial.interval),
        deepcopy(r.initial.last),
        deepcopy(r.initial.next),
        nothing,
    )
    cursor, base = 1, 1
    for j in eachindex(r.intervals)
        time = r.edges[j]
        for i in eachindex(mobile)
            mobile[i] || continue
            if (r.bootstrap && j == 1) || replay.next[i] == time
                cursor <= length(r.kicks) &&
                    r.kicks[cursor].species == i &&
                    r.kicks[cursor].time == time ||
                    throw(ArgumentError("persistent partition omitted a due species kick"))
                k = r.kicks[cursor]
                replay.interval[i], replay.last[i], replay.next[i] =
                    k.selected, time, time + Q(k.selected)
                cursor += 1
            end
        end
        edge, base = _empic_event_next(replay, mobile, time, r.finish, base, r.n_sub)
        edge == r.edges[j+1] ||
            throw(ArgumentError("persistent field partition differs from actual deadlines"))
    end
    cursor == length(r.kicks) + 1 ||
        throw(ArgumentError("persistent partition carries an unused kick"))
    return nothing
end

function _empic_event_phase_values(
    es,
    phase,
    species,
    initialized,
    time,
    step,
    full,
    fast;
    speed = es.c,
    fields = _empic_event_fields(es),
)
    T = eltype(es.rho_n)
    phase isa _EMPICEventPhase{T} &&
        _empic_event_q(phase.time) &&
        phase.time >= 0 &&
        initialized isa Bool &&
        time isa T &&
        step isa Int &&
        full isa T &&
        fast isa T || throw(ArgumentError("persistent phase scalar types differ"))
    _validate_empic_clock_phase(time, step, full, fast)
    count = length(es.collection.mobile)
    (species === nothing ? step == 0 : length(species) == count) &&
        all(a -> axes(a, 1) == Base.OneTo(count), (phase.interval, phase.last, phase.next)) &&
        all(isfinite, phase.interval) &&
        all(_empic_event_q, (phase.last..., phase.next...)) ||
        throw(ArgumentError("persistent phase population geometry differs"))
    if step == 0
        !phase.primed &&
            iszero(phase.time) &&
            _particle_stored_zero(time) &&
            all(_particle_stored_zero, phase.interval) &&
            all(iszero, (phase.last..., phase.next...)) &&
            phase.report === nothing ||
            throw(ArgumentError("step-zero persistent state carries an advanced phase"))
    else
        initialized &&
            phase.primed &&
            phase.time >= _exact_rational(full) &&
            phase.report !== nothing ||
            throw(ArgumentError("advanced persistent state lacks its phase/report"))
        for i in eachindex(es.collection.mobile)
            if es.collection.mobile[i]
                phase.interval[i] >= 2floatmin(T) &&
                    phase.last[i] <= phase.time <= phase.next[i] &&
                    phase.next[i] - phase.last[i] == _exact_rational(phase.interval[i]) ||
                    throw(ArgumentError("persistent species deadline/half phase differs"))
            else
                _particle_stored_zero(phase.interval[i]) &&
                    iszero(phase.last[i]) &&
                    iszero(phase.next[i]) ||
                    throw(ArgumentError("fixed species carries an active persistent phase"))
            end
        end
        _empic_event_report_valid(es, phase, species; speed, fields)
        r = phase.report
        isequal(r.public_after, time) &&
            isequal(r.requested, full) &&
            isequal(last(r.intervals), fast) &&
            r.bootstrap === (step == 1) ||
            throw(ArgumentError("persistent report differs from completed public phase"))
        if step == 1
            iszero(r.start) && _particle_stored_zero(r.public_before) ||
                throw(ArgumentError("first persistent interval has an invalid origin"))
        end
    end
    return nothing
end
function _empic_collection_validate_species_history(es, schedule::_EMPICEventSchedule)
    work = _empic_collection(es)
    _empic_event_phase_values(
        es,
        schedule.phase,
        work.initialized[] ? work.owners : nothing,
        work.initialized[],
        es.time[],
        es.step[],
        es.last_dt[],
        es.last_sub_dt[],
    )
    work.initialized[] && _empic_event_shadow_matches(es, schedule)
    return nothing
end
function _empic_collection_reset_intervals!(schedule::_EMPICEventSchedule{T}) where {T}
    schedule.candidate = _empic_event_zero_phase(T, length(schedule.proposal.cap))
    schedule.candidate_shadow = nothing
    return nothing
end
function _empic_collection_finish_candidate!(
    es,
    work::_EMPICCollectionWork{D,T,P,S,R},
    species,
) where {D,T,P,S,R<:_EMPICEventSchedule}
    work.schedule.candidate_shadow = _empic_event_shadow(es, work, species)
    return nothing
end
function _empic_collection_publish_intervals!(schedule::_EMPICEventSchedule)
    old, oldshadow = schedule.phase, schedule.shadow
    schedule.phase, schedule.shadow = schedule.candidate, schedule.candidate_shadow
    schedule.candidate, schedule.candidate_shadow = old, oldshadow
    return nothing
end
_empic_collection_prime_interval(schedule::_EMPICEventSchedule, i, fast) =
    schedule.candidate.interval[i]
_empic_collection_velocity_offsets(schedule::_EMPICEventSchedule, mobile, offset) = [
    mobile[i] && schedule.phase.primed ?
    (schedule.phase.last[i] + schedule.phase.next[i]) / 2 - schedule.phase.time :
    zero(_EMPICEventQ) for i in eachindex(mobile)
]
function _empic_collection_schedule_diagnostics(schedule::_EMPICEventSchedule, result)
    return merge(
        result,
        (
            schedule_policy = 2,
            source_time = deepcopy(schedule.phase.time),
            public_clock_defect = _exact_rational(schedule.shadow.clock[1]) - schedule.phase.time,
            species_dt = copy(schedule.proposal.cap),
            last_species_dt = copy(schedule.phase.interval),
            last_kick = deepcopy(schedule.phase.last),
            next_kick = deepcopy(schedule.phase.next),
            last_schedule_report = deepcopy(schedule.phase.report),
        ),
    )
end

function _empic_event_budget(schedule, intervals, kicks)
    budget = schedule.proposal.budget
    length(kicks) <= budget && length(intervals) <= budget - length(kicks) ||
        throw(ArgumentError("persistent event/counter budget exhausted"))
    return nothing
end
function _empic_event_next(phase, mobile, time, finish, base, M)
    edge = finish
    while base <= M
        proposed = phase.time + (finish - phase.time) * base / M
        if proposed > time
            edge = min(edge, proposed)
            break
        end
        base += 1
    end
    for i in eachindex(mobile)
        mobile[i] || continue
        phase.next[i] > time || throw(ArgumentError("persistent pending event did not advance"))
        edge = min(edge, phase.next[i])
    end
    return edge, base
end
function _empic_event_work!(
    impulse::Vector{_EMPICEventQ},
    current::Tuple,
    oldE::Tuple,
    newE::Tuple,
    h::T,
    volume::_EMPICEventQ,
) where {T<:AbstractFloat}
    Q, work = _exact_rational, zero(_EMPICEventQ)
    for c in eachindex(impulse), I in eachindex(impulse[c])
        amount = Q(h) * Q(current[c][I])
        impulse[c][I] += amount
        work += volume * amount * (Q(oldE[c][I]) + Q(newE[c][I])) / 2
    end
    return work
end
function _empic_event_advance!(
    es::Union{EMPIC,EMPIC1D},
    species::AbstractVector,
    work::_EMPICCollectionWork,
    charges,
    duration::T,
) where {T}
    schedule, Q = work.schedule, _exact_rational
    es.n_sub <= schedule.proposal.budget ||
        throw(ArgumentError("base field partition exceeds persistent event budget"))
    _empic_collection_prepare!(es, species)
    stage = work.stage
    original = schedule.phase
    H, a, b = copy(original.interval), deepcopy(original.last), deepcopy(original.next)
    schedule.candidate =
        _EMPICEventPhase{T}(deepcopy(original.time), original.primed, H, a, b, nothing)
    initial = (
        primed = original.primed,
        interval = copy(original.interval),
        last = deepcopy(original.last),
        next = deepcopy(original.next),
    )
    start, finish, time = original.time, original.time + Q(duration), original.time
    intervals, edges, defects, continuity = T[], _EMPICEventQ[deepcopy(time)], _EMPICEventQ[], T[]
    kicks, segments = NamedTuple[], NamedTuple[]
    bootstrap = !original.primed
    field_before = _empic_event_energy(stage)
    if bootstrap
        for (i, ps) in enumerate(work.staged)
            work.mobile[i] || continue
            _empic_collection_gather!(stage, ps; initial = true)
            decision = _empic_event_decision(stage, schedule, ps, i, time, zero(T), true)
            push!(kicks, decision)
            _empic_event_budget(schedule, intervals, kicks)
            H[i], a[i], b[i] = decision.selected, time, time + Q(decision.selected)
        end
    end
    base = 1
    previous = stage.last_sub_dt[]
    if bootstrap
        edge, base =
            _empic_event_next(schedule.candidate, work.mobile, time, finish, base, es.n_sub)
        first_h = _empic_event_interval(T, edge - time)
        _empic_collection_prime!(stage, work, first_h)
        previous = first_h
    end
    target = _empic_collection_charge_target(es, work, charges)
    impulse = ntuple(_ -> fill(zero(_EMPICEventQ), es.g.n), length(_empic_event_J(stage)))
    raw = _empic_event_raw(stage)
    raw_impulse = raw === nothing ? nothing : map(a -> fill(zero(_EMPICEventQ), size(a)), raw)
    used_work = zero(_EMPICEventQ)
    raw_work = raw === nothing ? nothing : zero(_EMPICEventQ)
    volume = prod(Q, es.g.dx)
    rho_start = similar(stage.rho_n)
    public_before = stage.time[]
    first_interval = true
    while time < finish
        parts = _empic_collection_begin_interval!(stage, work, target, previous)
        first_interval && copyto!(rho_start, stage.rho_n)
        for (i, ps) in enumerate(work.staged)
            work.mobile[i] || continue
            if bootstrap && first_interval
                # The initial decision was made on full-time values exactly once.
                _empic_collection_gather!(stage, ps)
                _empic_collection_push_gathered!(stage, ps, H[i])
            elseif b[i] == time
                _empic_collection_gather!(stage, ps)
                decision = _empic_event_decision(stage, schedule, ps, i, time, H[i], false)
                push!(kicks, decision)
                _empic_event_budget(schedule, intervals, kicks)
                _empic_collection_push_gathered!(stage, ps, decision.duration)
                H[i], a[i], b[i] = decision.selected, time, time + Q(decision.selected)
            end
        end
        edge, base =
            _empic_event_next(schedule.candidate, work.mobile, time, finish, base, es.n_sub)
        h = _empic_event_interval(T, edge - time)
        next_time, _ = _validated_step_clock(stage.time[], es.step[], h)
        next_time > stage.time[] || throw(ArgumentError("persistent public clock cannot advance"))
        push!(intervals, h)
        _empic_event_budget(schedule, intervals, kicks)
        oldE = map(copy, _empic_event_E(stage))
        for i in eachindex(work.mobile)
            work.mobile[i] && _empic_collection_drift!(stage, work, i, h)
        end
        _empic_collection_finish_interval!(stage, work, parts, previous, h, target)
        push!(
            segments,
            (
                E_before = oldE,
                E_after = map(copy, _empic_event_E(stage)),
                current = map(copy, _empic_event_J(stage)),
                raw_current = raw === nothing ? nothing : map(copy, raw),
                rho_before = copy(stage.rho_n),
                rho_after = copy(stage.rho_np1),
            ),
        )
        used_work += _empic_event_work!(
            impulse,
            _empic_event_J(stage),
            oldE,
            _empic_event_E(stage),
            h,
            volume,
        )
        if raw !== nothing
            raw_work += _empic_event_work!(raw_impulse, raw, oldE, _empic_event_E(stage), h, volume)
        end
        residual = charge_conservation_residual(stage, h)
        isfinite(residual) && residual >= zero(T) ||
            throw(ArgumentError("persistent actual continuity observation is unresolved"))
        push!(continuity, residual)
        push!(defects, Q(h) - (edge - time))
        push!(edges, deepcopy(edge))
        time, previous, first_interval = edge, h, false
        stage.time[] = next_time
    end
    _, count = _validated_step_clock(es.time[], es.step[], zero(T))
    stage.step[], stage.last_dt[], stage.last_sub_dt[] = count, duration, previous
    field_after = _empic_event_energy(stage)
    report = (
        kind = es isa EMPIC ? :vector : :planar,
        start = deepcopy(start),
        finish = deepcopy(finish),
        requested = duration,
        public_before = public_before,
        public_after = stage.time[],
        budget = schedule.proposal.budget,
        n_sub = es.n_sub,
        bootstrap = bootstrap,
        initial = initial,
        intervals = intervals,
        edges = edges,
        duration_defects = defects,
        kicks = kicks,
        continuity = continuity,
        charge_source_scale = charges.scale,
        segments = segments,
        rho_start = rho_start,
        rho_end = copy(stage.rho_np1),
        impulse = impulse,
        raw_impulse = raw_impulse,
        work_used = used_work,
        work_raw = raw_work,
        field_before = field_before,
        field_after = field_after,
        field_defect = field_after - field_before + used_work,
        clock_defect = Q(stage.time[]) - Q(public_before) - Q(duration),
    )
    schedule.candidate = _EMPICEventPhase{T}(finish, true, H, a, b, report)
    _empic_event_phase_values(
        es,
        schedule.candidate,
        work.staged,
        true,
        stage.time[],
        stage.step[],
        stage.last_dt[],
        stage.last_sub_dt[];
        fields = _empic_event_fields(stage),
    )
    _empic_collection_finish_candidate!(es, work, species)
    _empic_collection_validate_storage(es, species)
    return stage
end
