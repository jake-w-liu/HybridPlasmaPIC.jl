# First-order fixed-step composition of the canonical hybrid and local-Coulomb owners.

const _CollisionalPopulation{D,T} =
    ParticleSet{D,T,Vector{T},Vector{T},Vector{T},Vector{UInt64},Vector{UInt32}}

mutable struct _CollisionalTransaction{D,T,R,F}
    fields::HybridFields{D,T}
    Eforce::F
    x::NTuple{D,Vector{T}}
    v::NTuple{3,Vector{T}}
    rng::R
    clock::Tuple{T,Int,T}
end

"""
    CollisionalHybridStepper(inner::HybridStepper, species::AbstractVector;
        dt, rng, coulomb_log, epsilon0, weight_to_number, u_floor=1e-3,
        electron_ion=nothing)

Fixed-step, first-order CPU composition of [`HybridStepper`](@ref) with cell-local
[`collide_coulomb!`](@ref). Each nonzero `step!(stepper, species, dt)` pushes the
particles and advances B to `t+dt`, scatters the velocities for the interval `dt`
at the new positions (velocities stay at the half level `t+dt/2`), and then forms
the integer-time moments and carried E from the scattered velocities. The result
equals `step!(inner, ...)`, `collide_coulomb!(..., dt; ...)` and a recomputation
of the moments, Ohm's law and carried E, without forming E twice. Opt-in
`electron_ion=ElectronIonThermalExchange(...)` then exchanges stored electron
pressure and ion thermal energy before that recomputation. This is a fixed-step
first-order Lie composition. Call `init!` before active stepping.

Supports Float32/Float64, D=1:3, standard Vector particle storage, positive ions,
NGP/CIC/TSC, ideal massless `IsothermalElectrons`, `EvolvingElectrons` or `CGLElectrons`, and
`TransportElectrons`, `ConductingElectrons` or `FieldAlignedConductingElectrons`
with constant `η ≥ 0`, `ηH ≥ 0`
and `de2 = 0`. Conduction retains its original-quotient numerical policy and guard.
With stored pressure, initialize `inner.fields.pe` before `init!`; the field
substep evolves pressure and the ion-ion collision substep leaves it unchanged.
CGL keeps its algebraic global invariants and recomputes the full pressure-tensor
force at each trial-B stage. Its `pe` is scratch; no electron heat reservoir or
scalar thermal exchange is admitted. Full anisotropic fields are rolled back and
persisted under CGL schema 4. The existing force floor and raw-B energy limits
remain; this composition adds no global conservation or stability theorem.
Thermal exchange requires grid γ=T(5/3) and actual ne≥nfloor. NGP keeps its
simultaneous nodal star; CIC/TSC use ordered positive frame/node overlap edges;
it omits drag and finite electron momentum. Ohmic/hyperresistive heating and
conduction occur only in the field stages. Transport rolls back its independently
carried ion force as well as the full E. The ordered populations,
counts, charges, masses, model and timestep are fixed; any other nonzero `dt` is
rejected before state changes and `dt = 0` is a no-op that consumes no random
draws. With both sources inactive the physical evolution reproduces plain
`HybridStepper` stepping exactly. The wrapper borrows the inner stepper exclusively and deep-copies an
explicit Xoshiro or MersenneTwister RNG; independent wrappers need independent
inner steppers and populations.

A failed operation restores particle positions/velocities, every public field,
the clock, source ledger/report and RNG in place, so a retry draws the same numbers.
Private scratch (inner workspaces, collision bins) may change. Initialization
resets the physical clock but preserves the current RNG stream.
Checkpoint/restart retains the full advanced RNG state for the same Julia version
and environment; a seed alone is only provenance.
"""
mutable struct CollisionalHybridStepper{D,T,S,R,L,M,SH,G<:FourierGrid{D,T},F,X,W}
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
    const transaction::_CollisionalTransaction{D,T,R,F}
    const busy::Threads.Atomic{Bool}
    const electron_ion::X
    const thermal_work::W
    thermal_report::Union{Nothing,NamedTuple}
    phase::Symbol
    source_time::T
    source_count::Int
end

function _validate_collisional_model(inner::HybridStepper{D,T}) where {D,T}
    T <: Union{Float32,Float64} ||
        throw(ArgumentError("CollisionalHybridStepper requires Float32 or Float64"))
    model = inner.model
    model.closure isa
    Union{IsothermalElectrons,EvolvingElectrons,_ScalarTransportClosure,CGLElectrons} || throw(
        ArgumentError(
            "CollisionalHybridStepper requires IsothermalElectrons, EvolvingElectrons, TransportElectrons, ConductingElectrons, FieldAlignedConductingElectrons or CGLElectrons",
        ),
    )
    if _transports_electrons(model.closure)
        iszero(model.de2) || throw(ArgumentError("collisional scalar transport requires de2 = 0"))
    else
        (iszero(model.η) && iszero(model.ηH) && iszero(model.de2)) ||
            throw(ArgumentError("ideal collisional electrons require η = ηH = de2 = 0"))
    end
    if model.closure isa IsothermalElectrons
        _require_finite_nonnegative_real("Te", model.closure.Te, T)
    elseif model.closure isa CGLElectrons
        _validate_collisional_cgl_coefficients(model.closure)
    else
        _validate_evolving_model(model, T)
    end
    _require_finite_positive_real("nfloor", model.nfloor, T)
    inner.shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("CollisionalHybridStepper requires NGP, CIC or TSC"))
    return nothing
end

_collisional_closure_values(c::IsothermalElectrons) = (c.Te,)
_collisional_closure_values(c::Union{EvolvingElectrons,TransportElectrons}) = (c.γ,)
_collisional_closure_values(c::ConductingElectrons) = (c.γ, c.κ)
_collisional_closure_values(c::FieldAlignedConductingElectrons) =
    (c.γ, c.κ_parallel, c.κ_perpendicular)
_collisional_closure_values(c::CGLElectrons) = (c.Cperp, c.Cpar)

function _validate_collisional_cgl_coefficients(closure::CGLElectrons{C}) where {C}
    C in (Float16, Float32, Float64, BigFloat) ||
        throw(ArgumentError("collisional CGL requires standard floating invariant coefficients"))
    for coefficient in (closure.Cperp, closure.Cpar)
        isfinite(coefficient) && coefficient >= zero(C) ||
            throw(ArgumentError("collisional CGL invariants must be finite and nonnegative"))
    end
    return nothing
end

# Composite identity compares BigFloat backing storage, which changes on
# checkpoint copy. Bind the original scalar values instead of their addresses.
function _collisional_model_matches(a::HybridModel, b::HybridModel)
    return typeof(a) === typeof(b) && isequal(
        (a.η, a.ηH, a.nfloor, a.de2, _collisional_closure_values(a.closure)),
        (b.η, b.ηH, b.nfloor, b.de2, _collisional_closure_values(b.closure)),
    )
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

_collisional_rng_buffers(::Random.Xoshiro) = ()
_collisional_rng_buffers(rng::Random.MersenneTwister) = (rng.state.val, rng.vals, rng.ints)

function _validate_collisional_rng_pair(live, candidate, context::AbstractString)
    # Both stdlib owners validate fixed recurrence/cache sizes before any copy!.
    # copy!(MersenneTwister) shares seed metadata; draws never mutate that seed.
    _validate_collisional_rng(live)
    _validate_collisional_rng(candidate)
    live === candidate &&
        throw(ArgumentError("$context live and candidate RNGs must be independent"))
    for a in _collisional_rng_buffers(live), b in _collisional_rng_buffers(candidate)
        _particle_arrays_mightalias(a, b) &&
            throw(ArgumentError("$context RNG recurrence and caches must be independent"))
    end
    return nothing
end

# Stock metadata is read-only across draws; copy! may share it between streams.
_collisional_rng_metadata(::Random.Xoshiro) = ()
_collisional_rng_metadata(rng::Random.MersenneTwister) = (rng.seed, rng.adv_jump)

# Restore replaces stream values; only stock copy! destination storage is required.
_validate_collisional_rng_storage(::Random.Xoshiro) = nothing
function _validate_collisional_rng_storage(rng::Random.MersenneTwister)
    Random.DSFMT.DSFMT_state(rng.state.val)
    length(rng.vals) == Random.MT_CACHE_F && length(rng.ints) == (Random.MT_CACHE_I >> 4) ||
        throw(ArgumentError("MersenneTwister copy destination has invalid cache storage"))
    buffers = _collisional_rng_buffers(rng)
    for i in eachindex(buffers), j = 1:(i-1)
        _particle_arrays_mightalias(buffers[i], buffers[j]) &&
            throw(ArgumentError("MersenneTwister recurrence and cache buffers must be independent"))
    end
    return nothing
end

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
    electron_ion = nothing,
) where {D,T}
    _validate_collisional_model(inner)
    _validate_collisional_populations(species, inner.g)
    # Charge, mass and total count are bound here; repeated calls check their
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
    thermal = _bind_electron_ion(electron_ion, inner, species)
    _validate_collisional_rng(rng)
    _validate_checkpoint_configuration((;), inner, (rng,))
    owned_rng = deepcopy(rng)
    populations = Vector{_CollisionalPopulation{D,T}}(species)
    transaction = _CollisionalTransaction(
        HybridFields{D,T}(inner.g.n; anisotropic = is_anisotropic(inner.model.closure)),
        _transport_force_storage(inner.g, inner.model.closure),
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
        deepcopy(inner.model),
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
        thermal,
        _electron_ion_work(thermal, inner, populations),
        nothing,
        :uninitialized,
        zero(T),
        0,
    )
    _validate_collisional_fields(st)
    _validate_electron_ion_binding(st)
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
    # Pressure is persistent state for this closure. Reuse the inner owner's
    # positivity and independence checks before saving or mutating either state.
    _evolves_pressure(st.model.closure) && _validate_pressure_step_state(st.inner, st.species)
    _validate_collisional_force_storage(st, st.transaction.Eforce)
    _validate_collisional_cgl_storage(st, st.model.closure)
    return nothing
end

_validate_collisional_cgl_storage(st, ::ElectronClosure) = nothing
function _validate_collisional_cgl_storage(st, ::CGLElectrons)
    w = st.transaction
    saved = _collisional_field_arrays(values(_hybrid_field_state(w.fields)))
    live = _collisional_field_arrays(values(_hybrid_field_state(st.inner.fields)))
    for array in (saved..., live...)
        size(array) == st.grid.n || throw(
            DimensionMismatch("collisional CGL fields and rollback must have full grid layout"),
        )
    end
    count = sum(first, st.identities; init = 0)
    for array in (w.x..., w.v...)
        length(array) == count || throw(
            DimensionMismatch(
                "collisional CGL particle rollback length must match fixed populations",
            ),
        )
    end
    _validate_collisional_rng_pair(st.rng, w.rng, "collisional CGL rollback")
    # Enumerate logical snapshot roles: the graph deliberately permits the same
    # Array owner to recur, but two rollback roles cannot share that storage.
    buffers = (saved..., w.x..., w.v..., _collisional_rng_buffers(w.rng)...)
    for i in eachindex(buffers), j = 1:(i-1)
        _particle_arrays_mightalias(buffers[i], buffers[j]) &&
            throw(ArgumentError("collisional CGL rollback roles must have independent storage"))
    end
    # Stock RNG seed/advance metadata may be shared read-only by copy!. Only
    # recurrence/cache buffers enter the snapshot independence boundary.
    _validate_checkpoint_configuration_registered((
        w.fields,
        w.x,
        w.v,
        _collisional_rng_buffers(w.rng),
    ),) do remember, backing, owned
        _register_checkpoint_configuration_owned!(remember, backing, owned, (;), st.inner)
        _visit_checkpoint_objects(remember, (st.species, st.workspace, st.rng); backing)
    end
    return nothing
end

@inline _collisional_force(st::CollisionalHybridStepper) =
    _transports_electrons(st.model.closure) ? st.inner.pressure_work.Eforce : nothing

_validate_collisional_force_storage(st, ::Nothing) = nothing
function _validate_collisional_force_storage(st, saved::NTuple{3})
    _validate_transport_force_storage(st.inner, st.species, saved)
    for a in saved, b in _collisional_force(st)
        _particle_arrays_mightalias(a, b) &&
            throw(ArgumentError("collisional force rollback must have independent storage"))
    end
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
    _collisional_model_matches(st.inner.model, st.model) && isequal(st.inner.shape, st.shape) ||
        throw(ArgumentError("collisional model and shape are fixed"))
    isequal(st.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional Coulomb logarithms are fixed"))
    # The bound vector holds the same validated objects with a concrete eltype.
    _validate_collisional_populations(st.species, st.grid)
    _validate_collisional_fields(st)
    _validate_collisional_rng(st.rng)
    _validate_electron_ion_binding(st)
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
        isequal(inner.last_dt[], st.dt) && st.dt <= st.source_time
    ) || throw(ArgumentError("collisional timestep history does not match its fixed dt"))
    _validate_electron_ion_clock(st)
    return nothing
end

function _save_collisional_transaction!(st::CollisionalHybridStepper{D}) where {D}
    w, inner = st.transaction, st.inner
    _copy_hybrid_fields!(w.fields, inner.fields)
    _copy_transport_force!(w.Eforce, _collisional_force(st))
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
    _copy_transport_force!(_collisional_force(st), w.Eforce)
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
        st.thermal_report = nothing
        st.thermal_work === nothing || (st.thermal_work.pending = nothing)
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

# The canonical post-field source sees integer positions/pe and half-level v.
# Its default ion-only subflow leaves pe fixed; opt-in thermal exchange follows
# it before the existing owner recomputes moments, pressure force and carried E.
function _apply_hybrid_source!(st::CollisionalHybridStepper, species, dtT)
    if st.electron_ion === nothing || _collisional_active(st.coulomb_log)
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
    end
    _apply_electron_ion!(st.electron_ion, st, species, dtT)
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
        if _collisional_active(st.coulomb_log) || _electron_ion_active(st.electron_ion)
            _step_hybrid!(inner, st.species, h, NB, st)
        else
            _step_hybrid!(inner, st.species, h, NB, nothing)
        end
        _validate_collisional_populations(st.species, st.grid)
        _validate_collisional_fields(st)
        st.source_time = inner.time[]
        st.source_count = inner.step[]
        _commit_electron_ion_report!(st)
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

_apply_electron_ion!(::Nothing, st, species, h) = nothing
function _apply_electron_ion!(config::_ElectronIonThermalConfig, st, species, h)
    _electron_ion_active(config) || return nothing
    work = st.thermal_work
    work === nothing &&
        throw(ArgumentError("active electron-ion exchange requires thermal work storage"))
    config.policy == 2 && return _prepare_electron_ion_overlap!(
        config,
        work,
        st.inner,
        species,
        h,
        st.epsilon0,
        st.weight_to_number,
        st.rng,
    )
    return _prepare_electron_ion!(
        config,
        work,
        st.inner,
        species,
        h,
        st.epsilon0,
        st.weight_to_number,
        st.rng,
    )
end

function _validate_electron_ion_binding(st::CollisionalHybridStepper)
    config = st.electron_ion
    if config === nothing
        st.thermal_work === nothing && st.thermal_report === nothing ||
            throw(ArgumentError("disabled electron-ion exchange cannot own thermal state"))
        return nothing
    end
    _validate_electron_ion_config(config, st.inner, st.species)
    work = st.thermal_work
    _validate_coulomb_workspace(work.bins, st.species, st.grid, config.logarithm)
    arrays = (work.density, work.number_density, work.pressure)
    for i in eachindex(arrays)
        size(arrays[i]) == st.grid.n ||
            throw(DimensionMismatch("electron-ion scratch dimensions differ from the grid"))
        _validate_coulomb_work_array(arrays[i], st.species, config.logarithm)
        for j = 1:(i-1)
            _require_coulomb_independent(arrays[i], arrays[j])
        end
        for bins in work.bins.populations
            _require_coulomb_work_independent(arrays[i], bins)
        end
    end
    # Reuse the canonical native-resource/storage graph. Scratch must be separate
    # from fields, inner FFT/RK/predictor work, rollback, ion-ion work and RNG.
    _validate_checkpoint_configuration_registered((
        arrays...,
        work.bins,
        config.logarithm,
    )) do remember, backing, owned
        _register_checkpoint_configuration_owned!(remember, backing, owned, (;), st.inner)
        _visit_checkpoint_objects(
            remember,
            (st.transaction, st.workspace, st.rng, st.coulomb_log, st.thermal_report);
            backing,
        )
    end
    return nothing
end

function _validate_electron_ion_clock(st::CollisionalHybridStepper)
    st.electron_ion === nothing && return nothing
    report = st.thermal_report
    if st.source_count == 0
        report === nothing ||
            throw(ArgumentError("unadvanced electron-ion state cannot carry a report"))
    else
        report_keys =
            st.electron_ion.policy == 1 ? _ELECTRON_ION_REPORT_KEYS : _ELECTRON_ION_OVERLAP_KEYS
        _electron_ion_record(report, report_keys)
        report.policy === st.electron_ion.policy &&
            isequal(report.time, st.source_time) &&
            report.step === st.source_count &&
            isequal(report.dt, st.dt) ||
            throw(ArgumentError("electron-ion report does not match the source clock"))
    end
    return nothing
end

function _commit_electron_ion_report!(st::CollisionalHybridStepper)
    st.electron_ion === nothing && return nothing
    if _electron_ion_active(st.electron_ion)
        report = st.thermal_work.pending
        report !== nothing &&
            isequal(report.time, st.source_time) &&
            report.step == st.source_count ||
            throw(ArgumentError("electron-ion source did not prepare its completed report"))
        st.thermal_report = report
    else
        if st.electron_ion.policy == 2
            st.thermal_report = _electron_ion_overlap_empty_report(
                st.electron_ion,
                st.source_time,
                st.source_count,
                st.dt,
            )
            return nothing
        end
        st.thermal_report = (
            policy = 1,
            time = st.source_time,
            step = st.source_count,
            dt = st.dt,
            cells = NamedTuple[],
        )
    end
    return nothing
end
