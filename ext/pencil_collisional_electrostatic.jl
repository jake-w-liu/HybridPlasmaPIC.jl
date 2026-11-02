# Actual-pencil ES/Coulomb composition; collision physics stays in coulomb_local.jl.
mutable struct DistributedCollisionalElectrostaticPIC{T,S,P,R,L,G}
    const inner::S
    const g::G
    const species::Vector{P}
    const identities::Vector{Tuple{T,T}}
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

function _pencil_ces_cells(inner)
    owned = inner._layout.by_rank[inner._layout.rank]
    count = Base.checked_mul(Base.checked_mul(length(owned[1]), length(owned[2])), length(owned[3]))
    Base.checked_add(count, 1)
    return count
end

# Only this geometry seam differs from serial collisions. No global bins exist.
function _pencil_ces_cell(inner, ps, p)
    g, layout = inner.g, inner._layout
    index = HybridPlasmaPIC._particle_linear_cell(ps, g, p, HybridPlasmaPIC._particle_cell_strides(g.n))
    cell = CartesianIndices(g.n)[index]
    owned = layout.by_rank[layout.rank]
    _pencil_contains(owned, Tuple(cell)) || throw(
        ArgumentError("canonical collision cell is not owned by the migrated PIC rank"),
    )
    local_cell = ntuple(d -> cell[d] - first(owned[d]) + 1, 3)
    return LinearIndices(map(length, owned))[local_cell...]
end

function _pencil_ces_cell_admission(inner, species)
    for ps in species, p in eachindex(ps.weight)
        iszero(ps.weight[p]) || _pencil_ces_cell(inner, ps, p)
    end
    return nothing
end

# Register actual pencil/solver roots while retaining the canonical configuration
# traversal, native-resource rejection and serialization/alias semantics.
function _pencil_ces_configuration(inner, extras, configuration)
    return HybridPlasmaPIC._validate_checkpoint_configuration_registered(configuration) do remember, backing, owned
        HybridPlasmaPIC._visit_checkpoint_objects(remember, inner; backing)
        HybridPlasmaPIC._visit_checkpoint_objects(remember, extras; backing)
    end
end

function _pencil_ces_roles(st, extra = ())
    inner, g = st.inner, st.g
    writers = Any[parent(a) for a in (inner.E..., inner.rho, inner._candidate_E..., inner._candidate_rho,
        g.physical, g.spectral, g.spectral_work..., g.staging, g.staging_work...)]
    for populations in (st.species, inner._staged, extra), ps in populations
        append!(writers, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
    end
    for work in st.workspace.populations
        append!(writers, (work.velocity..., HybridPlasmaPIC._coulomb_index_arrays(work)...))
    end
    for rng in (st.rng, st.candidate_rng)
        append!(writers, HybridPlasmaPIC._collisional_rng_buffers(rng))
    end
    readonly = (g.midx..., g.kfull..., g.ik..., g.kvec..., inner._layout.by_rank,
        inner._layout.partitions..., inner._layout.starts..., inner._layout.lookup,
        inner._species_scalars, st.identities)
    for i in eachindex(writers)
        for j = 1:(i-1)
            HybridPlasmaPIC._particle_arrays_mightalias(writers[i], writers[j]) &&
                throw(ArgumentError("distributed collisional ES writable storage must be independent"))
        end
        for metadata in readonly
            HybridPlasmaPIC._particle_arrays_mightalias(writers[i], metadata) &&
                throw(ArgumentError("distributed collisional ES writes must not alias readonly metadata"))
        end
    end
    return nothing
end

function _pencil_ces_rngs(st; values = true)
    live, candidate = st.rng, st.candidate_rng
    live !== candidate || throw(ArgumentError("distributed collisional ES RNGs must be independent"))
    values && HybridPlasmaPIC._validate_collisional_rng_pair(live, candidate, "distributed collisional ES")
    for rng in (live, candidate)
        HybridPlasmaPIC._validate_collisional_rng_storage(rng)
        other = rng === live ? candidate : live
        values && _pencil_ces_configuration(st.inner, (st.species, st.identities, st.busy, st.workspace,
            st.coulomb_log, st.log_binding, HybridPlasmaPIC._collisional_rng_buffers(other)), (rng,))
        # Do not register own buffers against the whole RNG; inspect only its
        # read-only seed/history graph to reject self-cache and cross-cache links.
        values && _pencil_ces_configuration(st.inner, (st.species, st.identities, st.busy, st.workspace,
            st.coulomb_log, st.log_binding, HybridPlasmaPIC._collisional_rng_buffers(live),
            HybridPlasmaPIC._collisional_rng_buffers(candidate)), HybridPlasmaPIC._collisional_rng_metadata(rng))
    end
    if !values
        for a in HybridPlasmaPIC._collisional_rng_buffers(live), b in HybridPlasmaPIC._collisional_rng_buffers(candidate)
            HybridPlasmaPIC._particle_arrays_mightalias(a, b) &&
                throw(ArgumentError("distributed collisional ES RNG buffers must be independent"))
        end
    end
    return nothing
end

function _pencil_ces_binding(st::DistributedCollisionalElectrostaticPIC{T}, species; values = true) where {T}
    inner = st.inner
    _pencil_es_mode(T)
    st.g === inner.g || throw(ArgumentError("distributed collisional ES grid is fixed"))
    (!values || st.phase in (:uninitialized, :ready)) ||
        throw(ArgumentError("invalid distributed collisional ES phase"))
    _pencil_es_storage(inner)
    _pencil_es_species(inner, species; bound = true, values, factors = values)
    length(species) == length(st.species) == length(st.identities) ||
        throw(ArgumentError("distributed collisional ES ordered species are fixed"))
    for a in eachindex(species)
        species[a] === st.species[a] && isequal(inner._species_scalars[a], st.identities[a]) ||
            throw(ArgumentError("distributed collisional ES owners and q/m are fixed"))
    end
    HybridPlasmaPIC._particle_stored_zero(inner.background_charge) &&
        inner.neutrality === :require && inner.grid_dc === :project && isequal(inner.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("distributed collisional ES requires fixed zero-background neutral/projected ES"))
    HybridPlasmaPIC._require_resolved_dt("distributed collisional ES dt",
        HybridPlasmaPIC._require_finite_positive_real("dt", st.dt, T))
    HybridPlasmaPIC._require_finite_positive_real("epsilon0", st.epsilon0, T)
    HybridPlasmaPIC._require_finite_positive_real("weight_to_number", st.weight_to_number, T)
    HybridPlasmaPIC._require_finite_positive_real("u_floor", st.u_floor, T)
    isequal(st.epsilon0, st.weight_to_number) ||
        throw(ArgumentError("distributed collisional ES requires equal resolved epsilon0 and weight_to_number"))
    isequal(st.coulomb_log, st.log_binding) || throw(ArgumentError("distributed collisional ES logs are fixed"))
    HybridPlasmaPIC._validate_coulomb_logs(st.coulomb_log, length(species), T)
    for width in st.g.dx
        HybridPlasmaPIC._require_finite_positive_real("grid cell width", width, T)
    end
    values && HybridPlasmaPIC._validate_local_coulomb_species(species, Val(3), T)
    HybridPlasmaPIC._validate_coulomb_workspace(st.workspace, species, T, _pencil_ces_cells(inner), st.coulomb_log)
    rng_roots = values ? (st.rng, st.candidate_rng) :
        (HybridPlasmaPIC._collisional_rng_buffers(st.rng), HybridPlasmaPIC._collisional_rng_buffers(st.candidate_rng))
    _pencil_ces_configuration(inner, (st.species, st.identities, st.busy, st.workspace, rng_roots),
        (st.coulomb_log, st.log_binding))
    _pencil_ces_configuration(inner, (st.species, st.identities, st.busy, st.coulomb_log, st.log_binding,
        rng_roots), (st.workspace,))
    _pencil_ces_rngs(st; values)
    _pencil_ces_roles(st)
    return nothing
end

function _pencil_ces_boundary(st)
    inner = st.inner
    st.phase === :ready && inner.initialized[] ||
        throw(ArgumentError("initialize distributed collisional ES before use"))
    _pencil_es_phase(inner)
    isequal(st.source_time, inner.time[]) && st.source_count == inner.step[] ||
        throw(ArgumentError("distributed collisional ES source ledger differs from its inner clock"))
    (inner.step[] == 0 || isequal(inner.last_dt[], st.dt)) ||
        throw(ArgumentError("distributed collisional ES carried interval differs from fixed dt"))
    _finite_vector_magnitude(inner.E)
    _finite_magnitude(inner.rho)
    return nothing
end

function _pencil_ces_config(st)
    return (_pencil_es_bound_config(st.inner), st.dt, st.coulomb_log, st.epsilon0,
        st.weight_to_number, st.u_floor, typeof(st.rng), HybridPlasmaPIC._collisional_runtime())
end

function _pencil_ces_agreement(st, species, context; ready = true, values = true)
    return _agree_configuration(st.g.fft.comm, context) do
        _pencil_ces_binding(st, species; values)
        ready && _pencil_ces_boundary(st)
        (_pencil_ces_config(st), ready ? (st.phase, st.source_time, st.source_count) : ())
    end
end

# A failed collective acquisition releases only this call's locally acquired flag.
function _pencil_ces_operation(f, st, context)
    acquired = false
    try
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "exclusive acquisition") do
            Base.Threads.atomic_cas!(st.busy, false, true) &&
                throw(ArgumentError("concurrent or reentrant distributed collisional ES operation"))
            acquired = true
        end
        return f()
    finally
        acquired && (st.busy[] = false)
    end
end

"""
    distributed_collisional_electrostatic_pic(inner, species; dt, rng, coulomb_log,
        epsilon0, weight_to_number, u_floor=1e-3)

Collective fixed-step3D actual-pencil ES/Coulomb composition. Uses all-mobile
ordinary Vector species, zero background, physical neutrality and explicit grid
DC projection. Required positive resolved epsilon0 and weight_to_number must be
equal. Each rank owns a deeply detached Xoshiro or MersenneTwister; provide the
same algorithm, independently chosen states and an exclusively borrowed inner
model/grid/communicator. See docs/DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md.
"""
function HybridPlasmaPIC.distributed_collisional_electrostatic_pic(
    inner::DistributedElectrostaticPIC{T}, species; dt, rng, coulomb_log,
    epsilon0, weight_to_number, u_floor = 1e-3,
) where {T}
    context = "Distributed collisional ES construction"
    _pencil_es_agreement(inner, species, context; ready = false)
    parameters = _agree_configuration(inner.g.fft.comm, context) do
        _pencil_pic_phase_values(inner.initialized[], inner.primed[], inner.last_dt[], inner.time[], inner.step[])
        inner.step[] == 0 || throw(ArgumentError("construct collisional ES before inner stepping"))
        HybridPlasmaPIC._particle_stored_zero(inner.background_charge) && inner.neutrality === :require ||
            throw(ArgumentError("distributed collisional ES requires zero background and neutrality=:require"))
        h = HybridPlasmaPIC._require_resolved_dt(context, HybridPlasmaPIC._require_finite_positive_real("dt", dt, T))
        ep = HybridPlasmaPIC._require_finite_positive_real("epsilon0", epsilon0, T)
        ns = HybridPlasmaPIC._require_finite_positive_real("weight_to_number", weight_to_number, T)
        isequal(ep, ns) || throw(ArgumentError("epsilon0 and weight_to_number must resolve equally"))
        uf = HybridPlasmaPIC._require_finite_positive_real("u_floor", u_floor, T)
        HybridPlasmaPIC._validate_coulomb_logs(coulomb_log, length(species), T)
        rng isa Union{HybridPlasmaPIC.Random.Xoshiro,HybridPlasmaPIC.Random.MersenneTwister} ||
            throw(ArgumentError("distributed collisional ES requires stock Xoshiro or MersenneTwister"))
        HybridPlasmaPIC._validate_collisional_rng(rng)
        _pencil_ces_configuration(inner, (species, coulomb_log), (rng,))
        _pencil_ces_configuration(inner, (species, rng), (coulomb_log,))
        logs = coulomb_log isa Real ? T(coulomb_log) : Matrix{T}(coulomb_log)
        (h, ep, ns, uf, logs, typeof(rng), HybridPlasmaPIC._collisional_runtime())
    end
    return HybridPlasmaPIC._mpi_collective_stage(inner.g.fft.comm, context, "composition allocation") do
        h, ep, ns, uf, logs = parameters[1:5]
        populations = _PencilParticle{T}[ps for ps in species]
        owned_rng = deepcopy(rng)
        st = DistributedCollisionalElectrostaticPIC(inner, inner.g, populations,
            copy(inner._species_scalars), inner.charge_rtol, h, logs, copy(logs), ep, ns, uf,
            owned_rng, deepcopy(owned_rng), HybridPlasmaPIC._coulomb_workspace(populations, T, _pencil_ces_cells(inner)),
            Base.Threads.Atomic{Bool}(false), :uninitialized, zero(T), 0)
        _pencil_ces_binding(st, species)
        st
    end
end

function HybridPlasmaPIC.init!(st::DistributedCollisionalElectrostaticPIC{T}, species) where {T}
    context = "Distributed collisional ES initialization"
    return _pencil_ces_operation(st, context) do
        _pencil_ces_agreement(st, species, context; ready = false)
        initial_time, initial_count = zero(T), 0
        HybridPlasmaPIC.init_espic!(st.inner, species)
        st.source_time = initial_time
        st.source_count = initial_count
        st.phase = :ready
        st
    end
end

function HybridPlasmaPIC.step!(st::DistributedCollisionalElectrostaticPIC{T}, species, dt) where {T}
    context = "Distributed collisional ES step"
    return _pencil_ces_operation(st, context) do
        _pencil_ces_agreement(st, species, context)
        h, next_time, next_step = _agree_configuration(st.g.fft.comm, context) do
            dt isa Real || throw(ArgumentError("distributed collisional ES dt must be real"))
            h = HybridPlasmaPIC._require_resolved_dt(context,
                HybridPlasmaPIC._validated_nonnegative_dt(T, dt; name = context))
            HybridPlasmaPIC._particle_stored_zero(h) && !iszero(dt) &&
                throw(ArgumentError("nonzero distributed collisional ES dt must not underflow at grid precision"))
            HybridPlasmaPIC._particle_stored_zero(h) || isequal(h, st.dt) ||
                throw(ArgumentError("distributed collisional ES uses its fixed dt"))
            if HybridPlasmaPIC._particle_stored_zero(h)
                h, st.source_time, st.source_count
            else
                time, count = HybridPlasmaPIC._validated_step_clock(st.source_time, st.source_count, h)
                time > st.source_time || throw(ArgumentError("distributed collisional ES time must advance"))
                h, time, count
            end
        end
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "current ownership") do
            _pencil_owned_particles(species, st.g, st.inner._layout)
            _pencil_ces_cell_admission(st.inner, species)
        end
        HybridPlasmaPIC._particle_stored_zero(h) && return st
        active = HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "collision preparation") do
            active = HybridPlasmaPIC._coulomb_active(species, st.coulomb_log, T)
            active && copy!(st.candidate_rng, st.rng)
            active
        end
        prepared = _pencil_es_advance_candidate!(st.inner, species, h, context)
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "local collision candidate") do
            HybridPlasmaPIC._validate_local_coulomb_species(prepared, Val(3), T)
            HybridPlasmaPIC._validate_coulomb_workspace(st.workspace, prepared, T,
                _pencil_ces_cells(st.inner), st.coulomb_log)
            _pencil_ces_cell_admission(st.inner, prepared)
            _pencil_ces_roles(st, prepared)
            if active
                HybridPlasmaPIC._collide_coulomb_prepared!(prepared, st.g.dx, h, st.epsilon0,
                    st.weight_to_number, st.u_floor, st.coulomb_log, st.candidate_rng,
                    st.workspace, (ps, p) -> _pencil_ces_cell(st.inner, ps, p))
                HybridPlasmaPIC._validate_collisional_rng(st.candidate_rng)
            end
        end
        scalars = _pencil_es_final(st.inner, species, prepared, context)
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "complete publication preflight") do
            _pencil_ces_binding(st, species)
            _pencil_ces_roles(st, prepared)
        end
        active && copy!(st.rng, st.candidate_rng)
        _pencil_es_commit!(st.inner, species, prepared, scalars; dt = h, clock = (next_time, next_step))
        st.inner.primed[] = true
        st.source_time = next_time
        st.source_count = next_step
        st.phase = :ready
        st
    end
end

function HybridPlasmaPIC.electrostatic_diagnostics(st::DistributedCollisionalElectrostaticPIC, species)
    context = "Distributed collisional ES diagnostics"
    return _pencil_ces_operation(st, context) do
        _pencil_ces_agreement(st, species, context)
        merge(HybridPlasmaPIC.electrostatic_diagnostics(st.inner, species),
            (split_policy = :field_then_coulomb, source_time = st.source_time, source_count = st.source_count))
    end
end
