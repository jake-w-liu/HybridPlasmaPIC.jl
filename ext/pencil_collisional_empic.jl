# Actual-pencil EM/Coulomb composition; collision physics stays in coulomb_local.jl.
mutable struct DistributedCollisionalEMPIC{T,S,P,R,L,G}
    const inner::S
    const g::G
    const species::Vector{P}
    const identities::Vector{Tuple{T,T}}
    const c::T
    const charge_rtol::T
    const dt::T
    const coulomb_log::L
    const log_binding::L
    const epsilon0::T
    const weight_to_number::T
    const u_floor::T
    const rng::R
    const candidate_rng::R
    const workspace::HybridPlasmaPIC.CoulombCollisionWorkspace{T}
    const busy::Base.Threads.Atomic{Bool}
    phase::Symbol
    source_time::T
    source_count::Int
end

function _pencil_cem_roles(st, extra = ())
    inner, g = st.inner, st.g
    writers = Any[parent(a) for a in (_pencil_empic_arrays(inner)..., inner._source_spectral,
        g.physical, g.spectral, g.spectral_work..., g.staging, g.staging_work...)]
    append!(writers, _pencil_empic_midpoints(inner))
    push!(writers, inner._population_counts)
    readonly = (g.midx..., g.kfull..., g.ik..., g.kvec..., inner._layout.by_rank,
        inner._layout.partitions..., inner._layout.starts..., inner._layout.lookup,
        inner._species_scalars, st.identities)
    return _pencil_collision_roles(st, writers, readonly, extra)
end

function _pencil_cem_binding(st::DistributedCollisionalEMPIC{T}, species; values = true, counts = true) where {T}
    inner = st.inner
    _pencil_empic_mode(T)
    st.g === inner.g || throw(ArgumentError("distributed collisional EMPIC grid is fixed"))
    (!values || st.phase in (:uninitialized, :ready)) ||
        throw(ArgumentError("invalid distributed collisional EMPIC phase"))
    _pencil_empic_storage(inner)
    _pencil_empic_species(inner, species; bound = true, values, factors = values, counts)
    length(species) == length(st.species) == length(st.identities) ||
        throw(ArgumentError("distributed collisional EMPIC ordered species are fixed"))
    for a in eachindex(species)
        species[a] === st.species[a] && isequal(inner._species_scalars[a], st.identities[a]) ||
            throw(ArgumentError("distributed collisional EMPIC owners and q/m are fixed"))
    end
    HybridPlasmaPIC._particle_stored_zero(inner.background_charge) &&
        inner.neutrality === :require && inner.grid_dc === :project &&
        isequal(inner.c, st.c) && isequal(inner.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("distributed collisional EMPIC requires fixed zero-background neutral/projected EMPIC"))
    HybridPlasmaPIC._require_resolved_dt("distributed collisional EMPIC dt",
        HybridPlasmaPIC._require_finite_positive_real("dt", st.dt, T))
    st.dt < inner._cfl || throw(ArgumentError("fixed collisional EMPIC dt must be strictly sub-CFL"))
    HybridPlasmaPIC._require_finite_positive_real("epsilon0", st.epsilon0, T)
    HybridPlasmaPIC._require_finite_positive_real("weight_to_number", st.weight_to_number, T)
    HybridPlasmaPIC._require_finite_positive_real("u_floor", st.u_floor, T)
    isequal(st.epsilon0, st.weight_to_number) ||
        throw(ArgumentError("distributed collisional EMPIC requires equal resolved epsilon0 and weight_to_number"))
    isequal(st.coulomb_log, st.log_binding) || throw(ArgumentError("distributed collisional EMPIC logs are fixed"))
    HybridPlasmaPIC._validate_coulomb_logs(st.coulomb_log, length(species), T)
    for width in st.g.dx
        HybridPlasmaPIC._require_finite_positive_real("grid cell width", width, T)
    end
    values && HybridPlasmaPIC._validate_local_coulomb_species(species, Val(3), T)
    HybridPlasmaPIC._validate_coulomb_workspace(st.workspace, species, T, _pencil_collision_cells(inner), st.coulomb_log)
    rng_roots = values ? (st.rng, st.candidate_rng) :
        (HybridPlasmaPIC._collisional_rng_buffers(st.rng), HybridPlasmaPIC._collisional_rng_buffers(st.candidate_rng))
    _pencil_collision_configuration(inner, (st.species, st.identities, st.busy, st.workspace, rng_roots),
        (st.coulomb_log, st.log_binding))
    _pencil_collision_configuration(inner, (st.species, st.identities, st.busy, st.coulomb_log, st.log_binding,
        rng_roots), (st.workspace,))
    _pencil_collision_rngs(st; values)
    _pencil_cem_roles(st)
    return nothing
end

function _pencil_cem_boundary(st)
    inner = st.inner
    st.phase === :ready && inner.initialized[] ||
        throw(ArgumentError("initialize distributed collisional EMPIC before use"))
    _pencil_empic_phase(inner)
    isequal(st.source_time, inner.time[]) && st.source_count == inner.step[] ||
        throw(ArgumentError("distributed collisional EMPIC source ledger differs from its inner clock"))
    (inner.step[] == 0 || isequal(inner.last_dt[], st.dt)) ||
        throw(ArgumentError("distributed collisional EMPIC carried interval differs from fixed dt"))
    _pencil_empic_finite_fields(_pencil_empic_physical(inner))
    return nothing
end

function _pencil_cem_config(st)
    return (_pencil_empic_bound_config(st.inner), st.dt, st.coulomb_log, st.epsilon0,
        st.weight_to_number, st.u_floor, typeof(st.rng), HybridPlasmaPIC._collisional_runtime())
end

function _pencil_cem_agreement(st, species, context; ready = true, values = true)
    return _agree_configuration(st.g.fft.comm, context) do
        _pencil_cem_binding(st, species; values, counts = ready)
        ready && _pencil_cem_boundary(st)
        (_pencil_cem_config(st), ready ? (st.phase, st.source_time, st.source_count) : ())
    end
end

"""
    distributed_collisional_empic(inner, species; dt, rng, coulomb_log,
        epsilon0, weight_to_number, u_floor=1e-3)

Collective fixed-step3D actual-pencil nonrelativistic EMPIC/Coulomb composition.
Uses all-mobile ordinary Vector species, zero background, physical neutrality and
explicit grid DC projection. Required positive resolved epsilon0 and
weight_to_number must be equal; dt is resolved and strictly sub-CFL. Each rank
owns detached stock live/candidate RNGs. The canonical EM candidate completes
before local scattering at new positions/carried-half velocities; its eleven
fields and preceding drift current remain unchanged. This is a first-order Lie
split, without an exact discrete energy law. Borrow inner/grid/comm exclusively;
see docs/DISTRIBUTED_COLLISIONAL_EMPIC.md for phases, limits and restart.
"""
function HybridPlasmaPIC.distributed_collisional_empic(
    inner::DistributedEMPIC{T}, species; dt, rng, coulomb_log,
    epsilon0, weight_to_number, u_floor = 1e-3,
) where {T}
    context = "Distributed collisional EMPIC construction"
    _pencil_empic_agreement(inner, species, context; ready = false)
    parameters = _agree_configuration(inner.g.fft.comm, context) do
        _pencil_empic_phase_values(inner.initialized[], inner.primed[], inner.last_dt[],
            inner.time[], inner.step[], inner._cfl)
        inner.step[] == 0 || throw(ArgumentError("construct collisional EMPIC before inner stepping"))
        HybridPlasmaPIC._particle_stored_zero(inner.background_charge) && inner.neutrality === :require ||
            throw(ArgumentError("distributed collisional EMPIC requires zero background and neutrality=:require"))
        h = HybridPlasmaPIC._require_resolved_dt(context, HybridPlasmaPIC._require_finite_positive_real("dt", dt, T))
        h < inner._cfl || throw(ArgumentError("fixed collisional EMPIC dt must be strictly sub-CFL"))
        ep = HybridPlasmaPIC._require_finite_positive_real("epsilon0", epsilon0, T)
        ns = HybridPlasmaPIC._require_finite_positive_real("weight_to_number", weight_to_number, T)
        isequal(ep, ns) || throw(ArgumentError("epsilon0 and weight_to_number must resolve equally"))
        uf = HybridPlasmaPIC._require_finite_positive_real("u_floor", u_floor, T)
        HybridPlasmaPIC._validate_coulomb_logs(coulomb_log, length(species), T)
        rng isa Union{HybridPlasmaPIC.Random.Xoshiro,HybridPlasmaPIC.Random.MersenneTwister} ||
            throw(ArgumentError("distributed collisional EMPIC requires stock Xoshiro or MersenneTwister"))
        HybridPlasmaPIC._validate_collisional_rng(rng)
        _pencil_collision_configuration(inner, (species, coulomb_log), (rng,))
        _pencil_collision_configuration(inner, (species, rng), (coulomb_log,))
        logs = coulomb_log isa Real ? T(coulomb_log) : Matrix{T}(coulomb_log)
        (h, ep, ns, uf, logs, typeof(rng), HybridPlasmaPIC._collisional_runtime())
    end
    return HybridPlasmaPIC._mpi_collective_stage(inner.g.fft.comm, context, "composition allocation") do
        h, ep, ns, uf, logs = parameters[1:5]
        populations = _PencilParticle{T}[ps for ps in species]
        owned_rng = deepcopy(rng)
        st = DistributedCollisionalEMPIC(inner, inner.g, populations,
            copy(inner._species_scalars), inner.c, inner.charge_rtol, h, logs, copy(logs), ep, ns, uf,
            owned_rng, deepcopy(owned_rng), HybridPlasmaPIC._coulomb_workspace(populations, T, _pencil_collision_cells(inner)),
            Base.Threads.Atomic{Bool}(false), :uninitialized, zero(T), 0)
        _pencil_cem_binding(st, species; counts = false)
        st
    end
end

function HybridPlasmaPIC.init!(st::DistributedCollisionalEMPIC{T}, species) where {T}
    context = "Distributed collisional EMPIC initialization"
    return _pencil_collision_operation(st, context) do
        _pencil_cem_agreement(st, species, context; ready = false)
        initial_time, initial_count = zero(T), 0
        HybridPlasmaPIC.init_empic!(st.inner, species)
        st.source_time = initial_time
        st.source_count = initial_count
        st.phase = :ready
        st
    end
end

function HybridPlasmaPIC.step!(st::DistributedCollisionalEMPIC{T}, species, dt) where {T}
    context = "Distributed collisional EMPIC step"
    return _pencil_collision_operation(st, context) do
        _pencil_cem_agreement(st, species, context)
        h, next_time, next_step = _agree_configuration(st.g.fft.comm, context) do
            dt isa Real || throw(ArgumentError("distributed collisional EMPIC dt must be real"))
            h = HybridPlasmaPIC._require_resolved_dt(context,
                HybridPlasmaPIC._validated_nonnegative_dt(T, dt; name = context))
            HybridPlasmaPIC._particle_stored_zero(h) && !iszero(dt) &&
                throw(ArgumentError("nonzero distributed collisional EMPIC dt must not underflow at grid precision"))
            HybridPlasmaPIC._particle_stored_zero(h) || isequal(h, st.dt) ||
                throw(ArgumentError("distributed collisional EMPIC uses its fixed dt"))
            if HybridPlasmaPIC._particle_stored_zero(h)
                h, st.source_time, st.source_count
            else
                time, count = HybridPlasmaPIC._validated_step_clock(st.source_time, st.source_count, h)
                time > st.source_time || throw(ArgumentError("distributed collisional EMPIC time must advance"))
                h, time, count
            end
        end
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "current ownership") do
            _pencil_owned_particles(species, st.g, st.inner._layout)
            _pencil_collision_cell_admission(st.inner, species)
        end
        HybridPlasmaPIC._particle_stored_zero(h) && return st
        active = HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "collision preparation") do
            active = HybridPlasmaPIC._coulomb_active(species, st.coulomb_log, T)
            active && copy!(st.candidate_rng, st.rng)
            active
        end
        _pencil_pic_charges(st.inner, species, context)
        prepared = _pencil_empic_advance_candidate!(st.inner, species, h, next_time, next_step, context)
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "local collision candidate") do
            HybridPlasmaPIC._validate_local_coulomb_species(prepared, Val(3), T)
            HybridPlasmaPIC._validate_coulomb_workspace(st.workspace, prepared, T,
                _pencil_collision_cells(st.inner), st.coulomb_log)
            _pencil_collision_cell_admission(st.inner, prepared)
            _pencil_cem_roles(st, prepared)
            if active
                HybridPlasmaPIC._collide_coulomb_prepared!(prepared, st.g.dx, h, st.epsilon0,
                    st.weight_to_number, st.u_floor, st.coulomb_log, st.candidate_rng,
                    st.workspace, (ps, p) -> _pencil_collision_cell(st.inner, ps, p))
                HybridPlasmaPIC._validate_collisional_rng(st.candidate_rng)
            end
        end
        binding = _pencil_empic_final(st.inner, species, prepared, context)
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "complete publication preflight") do
            _pencil_cem_binding(st, species)
            _pencil_cem_roles(st, prepared)
        end
        active && copy!(st.rng, st.candidate_rng)
        _pencil_empic_commit!(st.inner, species, prepared, binding; dt = h, clock = (next_time, next_step))
        st.source_time = next_time
        st.source_count = next_step
        st.phase = :ready
        st
    end
end

function HybridPlasmaPIC.electromagnetic_diagnostics(st::DistributedCollisionalEMPIC, species)
    context = "Distributed collisional EMPIC diagnostics"
    return _pencil_collision_operation(st, context) do
        _pencil_cem_agreement(st, species, context)
        merge(HybridPlasmaPIC.electromagnetic_diagnostics(st.inner, species),
            (split_policy = :field_then_coulomb, source_time = st.source_time, source_count = st.source_count))
    end
end
