# First-order fixed-step composition of the canonical hybrid and local-Coulomb owners.

const _CollisionalPopulation{D,T} =
    ParticleSet{D,T,Vector{T},Vector{T},Vector{T},Vector{UInt64},Vector{UInt32}}

mutable struct _CollisionalTransaction{D,T,R}
    fields::HybridFields{D,T}
    x::NTuple{D,Vector{T}}
    v::NTuple{3,Vector{T}}
    rng::R
    clock::Tuple{T,Int,T}
end

"""
    CollisionalHybridStepper(inner::HybridStepper, species::AbstractVector;
        dt, rng, coulomb_log, epsilon0, weight_to_number, u_floor=1e-3)

Fixed-step, first-order CPU composition of [`HybridStepper`](@ref) with cell-local
[`collide_coulomb!`](@ref). Each nonzero `step!(stepper, species, dt)` pushes the
particles and advances B to `t+dt`, scatters the velocities for the interval `dt`
at the new positions (velocities stay at the half level `t+dt/2`), and then forms
the integer-time moments and carried E from the scattered velocities. The result
equals `step!(inner, ...)`, `collide_coulomb!(..., dt; ...)` and a recomputation
of the moments, Ohm's law and carried E, without forming E twice. This is not a
symmetric or second-order collision split. Call `init!` before active stepping.

Supports Float32/Float64, D=1:3, standard Vector particle storage, positive ions,
NGP/CIC/TSC, and ideal massless `IsothermalElectrons`. The ordered populations,
counts, charges, masses, model and timestep are fixed; any other nonzero `dt` is
rejected before state changes and `dt = 0` is a no-op that consumes no random
draws. All-zero Coulomb logarithms reproduce plain `HybridStepper` stepping
exactly. The wrapper borrows the inner stepper exclusively and deep-copies an
explicit Xoshiro or MersenneTwister RNG; independent wrappers need independent
inner steppers and populations.

A failed operation restores particle positions/velocities, every public field,
the clock, source ledger and RNG in place, so a retry draws the same numbers.
Private scratch (inner workspaces, collision bins) may change. Initialization
resets the physical clock but preserves the current RNG stream.
Checkpoint/restart retains the full advanced RNG state for the same Julia version
and environment; a seed alone is only provenance.
"""
mutable struct CollisionalHybridStepper{D,T,S,R,L,M,SH,G<:FourierGrid{D,T}}
    const inner::S
    # Concrete, so the collision call dispatches statically and allocates nothing.
    const grid::G
    const species::Vector{_CollisionalPopulation{D,T}}
    const identities::Vector{Tuple{Int,T,T}}
    const model::M
    const shape::SH
    const dt::T
    const coulomb_log::L
    const log_binding::L
    const epsilon0::T
    const weight_to_number::T
    const u_floor::T
    const rng::R
    const workspace::CoulombCollisionWorkspace{T}
    const transaction::_CollisionalTransaction{D,T,R}
    const busy::Threads.Atomic{Bool}
    phase::Symbol
    source_time::T
    source_count::Int
end

function _validate_collisional_model(inner::HybridStepper{D,T}) where {D,T}
    T <: Union{Float32,Float64} ||
        throw(ArgumentError("CollisionalHybridStepper requires Float32 or Float64"))
    model = inner.model
    model.closure isa IsothermalElectrons ||
        throw(ArgumentError("CollisionalHybridStepper requires IsothermalElectrons"))
    (iszero(model.η) && iszero(model.ηH) && iszero(model.de2)) ||
        throw(ArgumentError("CollisionalHybridStepper requires η = ηH = de2 = 0"))
    _require_finite_nonnegative_real("Te", model.closure.Te, T)
    _require_finite_positive_real("nfloor", model.nfloor, T)
    inner.shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("CollisionalHybridStepper requires NGP, CIC or TSC"))
    return nothing
end

function _validate_collisional_rng(rng::Random.Xoshiro)
    iszero(rng.s0 | rng.s1 | rng.s2 | rng.s3) &&
        throw(ArgumentError("Xoshiro must have a nonzero recurrence state"))
    return nothing
end

function _validate_collisional_rng(rng::Random.MersenneTwister)
    # Invoke the installed stdlib storage/cache contracts, without copying or
    # drawing. The temporary wrapper is unused; its arrays remain caller-owned.
    Random.DSFMT.DSFMT_state(rng.state.val)
    Random.MersenneTwister(
        rng.seed,
        rng.state,
        rng.vals,
        rng.ints,
        rng.idxF,
        rng.idxI,
        rng.adv,
        rng.adv_jump,
        rng.adv_vals,
        rng.adv_ints,
    )
    # mt_pop! consumes only the unexhausted [1,2) Float64 cache. Empty caches may
    # contain reset zeros; validating those as draws would reject valid seeds.
    for i = (rng.idxF+1):length(rng.vals)
        1.0 <= rng.vals[i] < 2.0 ||
            throw(ArgumentError("MersenneTwister has an invalid live Float64 cache"))
    end
    return nothing
end

_validate_collisional_rng(rng) =
    throw(ArgumentError("CollisionalHybridStepper requires Xoshiro or MersenneTwister"))

function _validate_collisional_populations(species, g::FourierGrid{D,T}) where {D,T}
    species isa AbstractVector{<:ParticleSet{D,T}} ||
        throw(ArgumentError("collisional populations must match grid dimension and precision"))
    for ps in species
        ps isa _CollisionalPopulation{D,T} ||
            throw(ArgumentError("collisional populations require standard Vector storage"))
    end
    return _validate_local_coulomb_species(species, g)
end

function CollisionalHybridStepper(
    inner::HybridStepper{D,T},
    species::AbstractVector;
    dt::Real,
    rng,
    coulomb_log,
    epsilon0::Real,
    weight_to_number::Real,
    u_floor::Real = 1e-3,
) where {D,T}
    _validate_collisional_model(inner)
    _validate_collisional_populations(species, inner.g)
    # Charge, q/m and total count are bound here; repeated calls check their
    # identity and delegate changing storage/values to the local collision guard.
    count = _validate_hybrid_species(species)
    (iszero(inner.time[]) && iszero(inner.step[]) && iszero(inner.last_dt[])) ||
        throw(ArgumentError("construct collisional stepping from an unadvanced inner stepper"))
    h = _require_finite_positive_real("collisional dt", dt, T)
    epsilon = _require_finite_positive_real("epsilon0", epsilon0, T)
    number_scale = _require_finite_positive_real("weight_to_number", weight_to_number, T)
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    _validate_coulomb_logs(coulomb_log, length(species), T)
    logs = coulomb_log isa Real ? T(coulomb_log) : Matrix{T}(coulomb_log)
    _validate_collisional_rng(rng)
    _validate_checkpoint_configuration((;), inner, (rng,))
    owned_rng = deepcopy(rng)
    populations = Vector{_CollisionalPopulation{D,T}}(species)
    transaction = _CollisionalTransaction(
        HybridFields{D,T}(inner.g.n),
        ntuple(_ -> zeros(T, count), D),
        ntuple(_ -> zeros(T, count), 3),
        deepcopy(owned_rng),
        (zero(T), 0, zero(T)),
    )
    st = CollisionalHybridStepper(
        inner,
        inner.g,
        populations,
        [(nparticles(ps), ps.q, ps.m) for ps in populations],
        inner.model,
        inner.shape,
        h,
        logs,
        copy(logs),
        epsilon,
        number_scale,
        uf,
        owned_rng,
        CoulombCollisionWorkspace(populations, inner.g),
        transaction,
        Threads.Atomic{Bool}(false),
        :uninitialized,
        zero(T),
        0,
    )
    _validate_collisional_fields(st)
    return st
end

# Flatten the statically known public-field groups without a boxed iterator or
# temporary Vector. All leaves have the same Array{T,D} type, including empties.
@inline _collisional_field_arrays(array::Array) = (array,)
@inline _collisional_field_arrays(::Tuple{}) = ()
@inline _collisional_field_arrays(groups::Tuple) =
    (_collisional_field_arrays(first(groups))..., _collisional_field_arrays(Base.tail(groups))...)

function _validate_collisional_fields(st::CollisionalHybridStepper)
    f, saved = st.inner.fields, st.transaction.fields
    arrays = _collisional_field_arrays(values(_hybrid_field_state(f)))
    references = _collisional_field_arrays(values(_hybrid_field_state(saved)))
    for i in eachindex(arrays)
        size(arrays[i]) == size(references[i]) ||
            throw(DimensionMismatch("collisional field dimensions must match the grid"))
        _require_all_finite("collisional field", arrays[i], "composition")
    end
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("collisional public fields must have independent storage"))
        end
        for ps in st.species
            _mightalias_particle_storage(arrays[i], ps) &&
                throw(ArgumentError("collisional public fields must not alias particle storage"))
        end
    end
    f.floor_count[] >= 0 || throw(ArgumentError("density floor count must be nonnegative"))
    return nothing
end

function _validate_collisional_binding(st::CollisionalHybridStepper, species)
    st.phase in (:ready, :uninitialized) ||
        throw(ArgumentError("collisional state is in flight or has an invalid phase"))
    length(species) == length(st.species) ||
        throw(ArgumentError("collisional species count is fixed"))
    axes(species, 1) == Base.OneTo(length(species)) ||
        throw(ArgumentError("collisional species must use one-based axes"))
    for a in eachindex(st.species)
        ps = species[a]
        ps === st.species[a] ||
            throw(ArgumentError("collisional species identities and order are fixed"))
        isequal((nparticles(ps), ps.q, ps.m), st.identities[a]) ||
            throw(ArgumentError("collisional particle counts, charges and masses are fixed"))
    end
    st.inner.g === st.grid || throw(ArgumentError("collisional grid is fixed"))
    _validate_collisional_model(st.inner)
    isequal(st.inner.model, st.model) && isequal(st.inner.shape, st.shape) ||
        throw(ArgumentError("collisional model and shape are fixed"))
    isequal(st.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional Coulomb logarithms are fixed"))
    # The bound vector holds the same validated objects with a concrete eltype.
    _validate_collisional_populations(st.species, st.grid)
    _validate_collisional_fields(st)
    _validate_collisional_rng(st.rng)
    return nothing
end

function _validate_collisional_boundary(st::CollisionalHybridStepper)
    st.phase === :ready || throw(ArgumentError("initialize collisional stepping before use"))
    inner = st.inner
    isequal(st.source_time, inner.time[]) && st.source_count == inner.step[] ||
        throw(ArgumentError("collisional source ledger does not match the field clock"))
    isfinite(st.source_time) && st.source_time >= 0 && st.source_count >= 0 ||
        throw(ArgumentError("collisional source ledger must be finite and nonnegative"))
    (
        st.source_count == 0 ? iszero(inner.last_dt[]) && iszero(st.source_time) :
        isequal(inner.last_dt[], st.dt)
    ) || throw(ArgumentError("collisional timestep history does not match its fixed dt"))
    return nothing
end

function _save_collisional_transaction!(st::CollisionalHybridStepper{D}) where {D}
    w, inner = st.transaction, st.inner
    _copy_hybrid_fields!(w.fields, inner.fields)
    offset = 0
    for ps in st.species
        count = nparticles(ps)
        for d = 1:D
            copyto!(w.x[d], offset + 1, ps.x[d], 1, count)
        end
        for c = 1:3
            copyto!(w.v[c], offset + 1, ps.v[c], 1, count)
        end
        offset += count
    end
    copy!(w.rng, st.rng)
    w.clock = (inner.time[], inner.step[], inner.last_dt[])
    return nothing
end

function _restore_collisional_transaction!(st::CollisionalHybridStepper{D}) where {D}
    w, inner = st.transaction, st.inner
    _copy_hybrid_fields!(inner.fields, w.fields)
    offset = 0
    for ps in st.species
        count = nparticles(ps)
        for d = 1:D
            copyto!(ps.x[d], 1, w.x[d], offset + 1, count)
        end
        for c = 1:3
            copyto!(ps.v[c], 1, w.v[c], offset + 1, count)
        end
        offset += count
    end
    copy!(st.rng, w.rng)
    inner.time[], inner.step[], inner.last_dt[] = w.clock
    return nothing
end

function _acquire_collisional!(st::CollisionalHybridStepper)
    Threads.atomic_cas!(st.busy, false, true) &&
        throw(ArgumentError("concurrent or reentrant collisional operation"))
    return nothing
end

function init!(st::CollisionalHybridStepper, species::AbstractVector)
    _acquire_collisional!(st)
    try
        return _init_collisional!(st, species)
    finally
        st.busy[] = false
    end
end

function _init_collisional!(st::CollisionalHybridStepper, species::AbstractVector)
    _validate_collisional_binding(st, species)
    _save_collisional_transaction!(st)
    previous_phase = st.phase
    previous_ledger = (st.source_time, st.source_count)
    st.phase = :inflight
    completed = false
    try
        init!(st.inner, st.species)
        _validate_collisional_fields(st)
        st.source_time = st.inner.time[]
        st.source_count = st.inner.step[]
        completed = true
        return st
    finally
        if !completed
            _restore_collisional_transaction!(st)
            st.source_time, st.source_count = previous_ledger
        end
        st.phase = completed ? :ready : previous_phase
    end
end

@inline _collisional_active(logs::Real) = !iszero(logs)
@inline _collisional_active(logs::AbstractMatrix) = any(!iszero, logs)

# The only HybridStepper source: scatter the half-level velocities at the new
# positions for the full field interval before the carried E is formed.
function _apply_hybrid_source!(st::CollisionalHybridStepper, species, dtT)
    collide_coulomb!(
        species,
        st.grid,
        dtT;
        coulomb_log = st.coulomb_log,
        epsilon0 = st.epsilon0,
        weight_to_number = st.weight_to_number,
        u_floor = st.u_floor,
        rng = st.rng,
        workspace = st.workspace,
    )
    return nothing
end

function step!(
    st::CollisionalHybridStepper{D,T},
    species::AbstractVector,
    dt::Real;
    NB::Integer = 1,
) where {D,T}
    h = _validated_step_dt(T, dt, NB; min_NB = 1, name = "collisional step!")
    iszero(h) &&
        !iszero(dt) &&
        throw(ArgumentError("nonzero collisional dt must not underflow at grid precision"))
    iszero(h) && return st
    _acquire_collisional!(st)
    try
        return _step_collisional!(st, species, h, NB)
    finally
        st.busy[] = false
    end
end

function _step_collisional!(st::CollisionalHybridStepper, species, h, NB)
    isequal(h, st.dt) || throw(ArgumentError("collisional step! requires its fixed dt"))
    _validate_collisional_binding(st, species)
    _validate_collisional_boundary(st)
    _validated_step_clock(st.inner.time[], st.inner.step[], h)
    _save_collisional_transaction!(st)
    previous_ledger = (st.source_time, st.source_count)
    st.phase = :inflight
    completed = false
    try
        inner = st.inner
        if _collisional_active(st.coulomb_log)
            _step_hybrid!(inner, st.species, h, NB, st)
        else
            _step_hybrid!(inner, st.species, h, NB, nothing)
        end
        _validate_collisional_populations(st.species, st.grid)
        _validate_collisional_fields(st)
        st.source_time = inner.time[]
        st.source_count = inner.step[]
        completed = true
        return st
    finally
        if !completed
            _restore_collisional_transaction!(st)
            st.source_time, st.source_count = previous_ledger
        end
        st.phase = :ready
    end
end
