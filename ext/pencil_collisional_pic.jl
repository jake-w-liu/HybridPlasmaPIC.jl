# Bounded shared ownership for actual-pencil ES/EM Coulomb composition.

function _pencil_collision_cells(inner)
    owned = inner._layout.by_rank[inner._layout.rank]
    count = Base.checked_mul(Base.checked_mul(length(owned[1]), length(owned[2])), length(owned[3]))
    Base.checked_add(count, 1)
    return count
end

# Only this geometry seam differs from serial collisions. No global bins exist.
function _pencil_collision_cell(inner, ps, p)
    g, layout = inner.g, inner._layout
    index =
        HybridPlasmaPIC._particle_linear_cell(ps, g, p, HybridPlasmaPIC._particle_cell_strides(g.n))
    cell = CartesianIndices(g.n)[index]
    owned = layout.by_rank[layout.rank]
    _pencil_contains(owned, Tuple(cell)) ||
        throw(ArgumentError("canonical collision cell is not owned by the migrated PIC rank"))
    local_cell = ntuple(d -> cell[d] - first(owned[d]) + 1, 3)
    return LinearIndices(map(length, owned))[local_cell...]
end

function _pencil_collision_cell_admission(inner, species)
    for ps in species, p in eachindex(ps.weight)
        iszero(ps.weight[p]) || _pencil_collision_cell(inner, ps, p)
    end
    return nothing
end

# Register actual pencil/solver roots while retaining the canonical configuration
# traversal, native-resource rejection and serialization/alias semantics.
function _pencil_collision_configuration(inner, extras, configuration)
    return HybridPlasmaPIC._validate_checkpoint_configuration_registered(
        configuration,
    ) do remember, backing, owned
        HybridPlasmaPIC._visit_checkpoint_objects(remember, inner; backing)
        HybridPlasmaPIC._visit_checkpoint_objects(remember, extras; backing)
    end
end

function _pencil_collision_rngs(st; values = true)
    live, candidate = st.rng, st.candidate_rng
    live !== candidate ||
        throw(ArgumentError("distributed collisional PIC RNGs must be independent"))
    values && HybridPlasmaPIC._validate_collisional_rng_pair(
        live,
        candidate,
        "distributed collisional PIC",
    )
    for rng in (live, candidate)
        HybridPlasmaPIC._validate_collisional_rng_storage(rng)
        other = rng === live ? candidate : live
        values && _pencil_collision_configuration(
            st.inner,
            (
                st.species,
                st.identities,
                st.busy,
                st.workspace,
                st.coulomb_log,
                st.log_binding,
                HybridPlasmaPIC._collisional_rng_buffers(other),
            ),
            (rng,),
        )
        # Do not register own buffers against the whole RNG; inspect only its
        # read-only seed/history graph to reject self-cache and cross-cache links.
        values && _pencil_collision_configuration(
            st.inner,
            (
                st.species,
                st.identities,
                st.busy,
                st.workspace,
                st.coulomb_log,
                st.log_binding,
                HybridPlasmaPIC._collisional_rng_buffers(live),
                HybridPlasmaPIC._collisional_rng_buffers(candidate),
            ),
            HybridPlasmaPIC._collisional_rng_metadata(rng),
        )
    end
    if !values
        for a in HybridPlasmaPIC._collisional_rng_buffers(live),
            b in HybridPlasmaPIC._collisional_rng_buffers(candidate)

            HybridPlasmaPIC._particle_arrays_mightalias(a, b) &&
                throw(ArgumentError("distributed collisional PIC RNG buffers must be independent"))
        end
    end
    return nothing
end

# A failed collective acquisition releases only this call's locally acquired flag.
function _pencil_collision_operation(f, st, context)
    acquired = false
    try
        HybridPlasmaPIC._mpi_collective_stage(st.g.fft.comm, context, "exclusive acquisition") do
            Base.Threads.atomic_cas!(st.busy, false, true) && throw(
                ArgumentError("concurrent or reentrant distributed collisional PIC operation"),
            )
            acquired = true
        end
        return f()
    finally
        acquired && (st.busy[] = false)
    end
end

function _pencil_collision_roles(st, writers, readonly, extra = ())
    for populations in (st.species, st.inner._staged, extra), ps in populations
        append!(writers, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
    end
    for work in st.workspace.populations
        append!(writers, (work.velocity..., HybridPlasmaPIC._coulomb_index_arrays(work)...))
    end
    for rng in (st.rng, st.candidate_rng)
        append!(writers, HybridPlasmaPIC._collisional_rng_buffers(rng))
    end
    for i in eachindex(writers)
        for j = 1:(i-1)
            HybridPlasmaPIC._particle_arrays_mightalias(writers[i], writers[j]) && throw(
                ArgumentError("distributed collisional PIC writable storage must be independent"),
            )
        end
        for metadata in readonly
            HybridPlasmaPIC._particle_arrays_mightalias(writers[i], metadata) && throw(
                ArgumentError(
                    "distributed collisional PIC writes must not alias readonly metadata",
                ),
            )
        end
    end
    return nothing
end

# Full incoming values, including read-only seed/history graphs, precede detach.
function _pencil_collision_incoming_configuration(st, inner, logs, rng)
    typeof(rng) === typeof(st.rng) ||
        throw(ArgumentError("Pencil collisional PIC RNG algorithm differs"))
    HybridPlasmaPIC._validate_collisional_rng(rng)
    HybridPlasmaPIC._validate_collisional_rng_storage(rng)
    # Check original incoming links before detachment can hide them. Old stream
    # values may be stale, but incoming RNG must be detached from both streams.
    _pencil_collision_configuration(
        st.inner,
        (
            inner,
            logs,
            st.species,
            st.identities,
            st.busy,
            st.workspace,
            st.coulomb_log,
            st.log_binding,
            st.rng,
            st.candidate_rng,
        ),
        (rng,),
    )
    _pencil_collision_configuration(
        st.inner,
        (
            inner,
            logs,
            st.species,
            st.identities,
            st.busy,
            st.workspace,
            st.coulomb_log,
            st.log_binding,
            st.rng,
            st.candidate_rng,
            HybridPlasmaPIC._collisional_rng_buffers(rng),
        ),
        HybridPlasmaPIC._collisional_rng_metadata(rng),
    )
    _pencil_collision_configuration(
        st.inner,
        (inner, rng, st.species, st.identities, st.busy, st.workspace, st.rng, st.candidate_rng),
        (logs,),
    )
    return nothing
end

# Common fixed-collision payload metadata; family tags/inner schemas stay explicit.
function _pencil_collision_checkpoint_config(st, inner_config, model)
    return (
        schema = 1,
        model = model,
        split_policy = :field_then_coulomb_v1,
        collision_position_policy = :new_integer,
        collision_velocity_policy = :carried_half,
        collision_cell_policy = :canonical_owned_floor_v1,
        inner_config = inner_config,
        dt = st.dt,
        coulomb_log = copy(st.coulomb_log),
        epsilon0 = st.epsilon0,
        weight_to_number = st.weight_to_number,
        u_floor = st.u_floor,
        rng_algorithm = typeof(st.rng),
        runtime = HybridPlasmaPIC._collisional_runtime(),
    )
end

_pencil_collision_checkpoint_phase(st) = (
    initialized = st.inner.initialized[],
    primed = st.inner.primed[],
    last_dt = st.inner.last_dt[],
    time = st.inner.time[],
    step = st.inner.step[],
    phase = st.phase,
    source_time = st.source_time,
    source_count = st.source_count,
)

function _pencil_collision_checkpoint_common(record, st, expected, ::Type{T}) where {T}
    for name in keys(expected)
        value, reference = getproperty(record, name), getproperty(expected, name)
        typeof(value) === typeof(reference) && isequal(value, reference) ||
            throw(ArgumentError("Pencil collisional PIC checkpoint $name does not match"))
    end
    record.initialized === true &&
        record.primed isa Bool &&
        record.last_dt isa T &&
        record.time isa T &&
        record.step isa Int &&
        record.phase === :ready &&
        record.source_time isa T &&
        record.source_count isa Int || throw(
        ArgumentError("Pencil collisional PIC checkpoint requires a typed completed phase/clock"),
    )
    _pencil_pic_phase_values(
        record.initialized,
        record.primed,
        record.last_dt,
        record.time,
        record.step,
    )
    isequal(record.source_time, record.time) &&
        record.source_count == record.step &&
        (record.step == 0 || isequal(record.last_dt, st.dt)) || throw(
        ArgumentError("Pencil collisional PIC checkpoint source ledger or fixed interval differs"),
    )
    return nothing
end
