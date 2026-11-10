# Dynamic topology composes canonical static-composite epochs. Every numerical
# step/regrid is prepared privately; no field equation or particle mover is copied.
struct DynamicCompositeESRegridReport{T<:Union{Float32,Float64}}
    time::T
    step::Int
    regrid::Int
    before_refinements::Tuple
    after_refinements::Tuple
    before_leaves::Int
    after_leaves::Int
    field_before::_CompositeQ
    field_after::_CompositeQ
    kinetic_before::_CompositeQ
    kinetic_after::_CompositeQ
    mesh_work::_CompositeQ
    particle_charge_change::_CompositeQ
    background_charge_change::_CompositeQ
    left_plate_change::_CompositeQ
    right_plate_change::_CompositeQ
    max_gauss_residual::_CompositeQ
    max_constitutive_residual::_CompositeQ
    selection::Symbol
    threshold::Union{Nothing,T}
    buffer::Int
end

mutable struct DynamicCompositeElectrostaticPIC{
    T<:Union{Float32,Float64},
    S<:YeeTriangularCloud,
    B,
    P,
    F<:Tuple,
}
    mesh::StaticAMRHierarchy1D{T}
    E::Vector{T}
    phi::Vector{T}
    time::T
    step::Int
    last_dt::T
    initialized::Bool
    epoch_start_time::T
    epoch_start_step::Int
    regrid_count::Int
    initial_energy::_CompositeQ
    physical_defect::_CompositeQ
    mesh_work::_CompositeQ
    last_regrid_report::Union{Nothing,DynamicCompositeESRegridReport{T}}
    shape::S
    background_charge::B
    left_plate_charge::P
    _factors::F
    _configuration::NamedTuple
    _state::StaticCompositeElectrostaticPIC{T}
    _history::Union{Nothing,NamedTuple}
    _owners::Vector{_ElectrostaticHostParticle{1,T}}
    _components::Vector{Tuple}
    _counts::Vector{Int}
    _busy::Threads.Atomic{Bool}

    # Explicit field types keep the implicit constructor's array/range
    # conversions out of inference.
    function DynamicCompositeElectrostaticPIC{T,S,B,P,F}(
        mesh::StaticAMRHierarchy1D{T},
        E::Vector{T},
        phi::Vector{T},
        time::T,
        step::Int,
        last_dt::T,
        initialized::Bool,
        epoch_start_time::T,
        epoch_start_step::Int,
        regrid_count::Int,
        initial_energy::_CompositeQ,
        physical_defect::_CompositeQ,
        mesh_work::_CompositeQ,
        last_regrid_report::Union{Nothing,DynamicCompositeESRegridReport{T}},
        shape::S,
        background_charge::B,
        left_plate_charge::P,
        _factors::F,
        _configuration::NamedTuple,
        _state::StaticCompositeElectrostaticPIC{T},
        _history::Union{Nothing,NamedTuple},
        _owners::Vector{_ElectrostaticHostParticle{1,T}},
        _components::Vector{Tuple},
        _counts::Vector{Int},
        _busy::Threads.Atomic{Bool},
    ) where {T,S,B,P,F}
        return new{T,S,B,P,F}(
            mesh,
            E,
            phi,
            time,
            step,
            last_dt,
            initialized,
            epoch_start_time,
            epoch_start_step,
            regrid_count,
            initial_energy,
            physical_defect,
            mesh_work,
            last_regrid_report,
            shape,
            background_charge,
            left_plate_charge,
            _factors,
            _configuration,
            _state,
            _history,
            _owners,
            _components,
            _counts,
            _busy,
        )
    end
end

function _dce_configuration(st::StaticCompositeElectrostaticPIC)
    # Retain the canonical static policy and original scalar bits. Only topology
    # becomes variable; all other constructor semantics stay in that owner.
    canonical = merge(_sce_configuration(st), (refinements = (),))
    return (
        canonical = canonical,
        remap = :original_cloud_constraint_mesh_exchange_v1,
        tagging = :base_electric_curvature_v1,
        phase = :global_fixed_duration_static_epochs_v1,
    )
end

function _dce_copy_species(species, ::Type{T}) where {T}
    copies = _ElectrostaticHostParticle{1,T}[]
    _pic_prepare_species_stages!(copies, species)
    return copies
end

_dce_private_species(st::StaticCompositeElectrostaticPIC) = st._binding[].owners
_dce_energy(st::StaticCompositeElectrostaticPIC) =
    _sce_field_energy(st.mesh, st.E) + _kinetic_energy_rational(_dce_private_species(st))
_dce_particle_charge(species) = sum(
    (_exact_rational(ps.q) * _exact_float_sum(ps.weight) for ps in species);
    init = zero(_CompositeQ),
)

function _dce_species_equal(left, right)
    length(left) == length(right) || return false
    for (a, b) in zip(left, right)
        isequal((a.q, a.m), (b.q, b.m)) || return false
        all(
            pair -> isequal(pair[1], pair[2]),
            zip(_sce_particle_arrays(a), _sce_particle_arrays(b)),
        ) || return false
    end
    return true
end

function _dce_acquire(st::DynamicCompositeElectrostaticPIC)
    Threads.atomic_cas!(st._busy, false, true) == false ||
        throw(ArgumentError("dynamic composite ES already has an active operation"))
    return nothing
end
_dce_release(st::DynamicCompositeElectrostaticPIC) =
    (Threads.atomic_xchg!(st._busy, false); nothing)

function _dce_independent(arrays)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("dynamic composite physical/private/history storage roles overlap"))
    end
    return nothing
end

function _dce_owned_arrays(st::DynamicCompositeElectrostaticPIC, species)
    private = _dce_private_species(st._state)
    writable, readonly = _sce_arrays(st._state, private)
    arrays = vcat(writable, readonly)
    for ps in species
        append!(arrays, _sce_particle_arrays(ps))
    end
    push!(arrays, st._counts)
    if st._history !== nothing
        append!(arrays, _dce_history_arrays(st._history, eltype(st.E)))
    end
    return arrays
end

function _dce_storage(
    st::DynamicCompositeElectrostaticPIC{T},
    species;
    counts = true,
    values = true,
) where {T}
    _sce_mode()
    _sce_population_storage(species, T)
    length(species) ==
    length(st._owners) ==
    length(st._components) ==
    length(st._counts) ==
    length(st._factors) ||
        throw(ArgumentError("dynamic composite ordered population count changed; call init!"))
    all(component -> length(component) == 7, st._components) ||
        throw(ArgumentError("dynamic composite component binding layout changed"))
    st.mesh === st._state.mesh && st.E === st._state.E && st.phi === st._state.phi ||
        throw(ArgumentError("dynamic composite physical views no longer match their epoch"))
    isequal(st._configuration, _dce_configuration(st._state)) &&
        isequal(st._factors, st._state._factors) &&
        isequal(st.shape, st._state.shape) &&
        isequal(st.background_charge, st._state.background_charge) &&
        isequal(st.left_plate_charge, st._state.left_plate_charge) ||
        throw(ArgumentError("dynamic composite immutable configuration changed"))
    private = _dce_private_species(st._state)
    _sce_restore_storage(st._state, private)
    for (k, ps) in enumerate(species)
        ps === st._owners[k] ||
            throw(ArgumentError("dynamic composite population object/order changed; call init!"))
        all(pair -> pair[1] === pair[2], zip(_sce_particle_arrays(ps), st._components[k])) || throw(
            ArgumentError("dynamic composite particle component identities changed; call init!"),
        )
        counts &&
            length(ps.weight) != st._counts[k] &&
            throw(ArgumentError("dynamic composite particle count changed; call init!"))
    end
    values &&
        !_dce_species_equal(species, private) &&
        throw(ArgumentError("dynamic composite physical/source values changed; call init!"))
    _dce_independent(_dce_owned_arrays(st, species))
    return nothing
end

function _dce_epoch_clock(
    time::T,
    step::Int,
    dt::T,
    epoch_time::T,
    epoch_step::Int,
    local_time::T,
    local_step::Int,
    local_dt::T,
) where {T}
    isfinite(time) &&
        _stored_nonnegative(time) &&
        step >= 0 &&
        isfinite(dt) &&
        _stored_nonnegative(dt) &&
        isfinite(epoch_time) &&
        _stored_nonnegative(epoch_time) &&
        0 <= epoch_step <= step &&
        epoch_time <= time &&
        0 <= local_step <= step &&
        local_step == step - epoch_step &&
        isfinite(local_time) &&
        _stored_nonnegative(local_time) &&
        local_time <= time ||
        throw(ArgumentError("dynamic composite global/epoch clock is inconsistent"))
    if iszero(step)
        _particle_stored_zero(time) && _particle_stored_zero(dt) ||
            throw(ArgumentError("unadvanced dynamic composite global clock/history must be zero"))
    else
        _stored_positive(time) && _stored_positive(dt) && dt <= time ||
            throw(ArgumentError("advanced dynamic composite duration must be positive and <= time"))
    end
    if iszero(epoch_step)
        _particle_stored_zero(epoch_time) ||
            throw(ArgumentError("unadvanced dynamic epoch origin must have stored-zero time"))
    else
        _stored_positive(epoch_time) && dt <= epoch_time || throw(
            ArgumentError("advanced dynamic epoch origin must retain a completed fixed duration"),
        )
    end
    if iszero(local_step)
        isequal(time, epoch_time) &&
            step == epoch_step &&
            _particle_stored_zero(local_time) &&
            _particle_stored_zero(local_dt) ||
            throw(ArgumentError("empty dynamic epoch must retain its exact global origin"))
    else
        time > epoch_time &&
            _stored_positive(local_time) &&
            isequal(local_dt, dt) &&
            local_dt <= local_time ||
            throw(ArgumentError("advanced dynamic epoch has inconsistent fixed-duration history"))
    end
    return nothing
end

function _dce_regrid_report(
    st::DynamicCompositeElectrostaticPIC{T},
    report,
    epoch_time,
    epoch_step,
    count,
    mesh,
) where {T}
    report isa DynamicCompositeESRegridReport{T} &&
        report.regrid == count &&
        report.step == epoch_step &&
        isequal(report.time, epoch_time) &&
        isequal(report.after_refinements, mesh.refinements) &&
        report.after_leaves == length(mesh.widths) ||
        throw(ArgumentError("dynamic composite regrid report does not match the current epoch"))
    before = StaticAMRHierarchy1D(mesh.nbase, mesh.L, report.before_refinements)
    report.before_leaves == length(before.widths) &&
        !isequal(report.before_refinements, report.after_refinements) ||
        throw(ArgumentError("dynamic composite regrid report has invalid topology provenance"))
    for name in (
        :field_before,
        :field_after,
        :kinetic_before,
        :kinetic_after,
        :mesh_work,
        :particle_charge_change,
        :background_charge_change,
        :left_plate_change,
        :right_plate_change,
        :max_gauss_residual,
        :max_constitutive_residual,
    )
        _sce_q_finite(getfield(report, name)) ||
            throw(ArgumentError("dynamic composite regrid report requires finite exact scalars"))
    end
    min(
            report.field_before,
            report.field_after,
            report.kinetic_before,
            report.kinetic_after,
            report.max_gauss_residual,
            report.max_constitutive_residual,
        ) >= 0 &&
        report.mesh_work == report.field_after - report.field_before &&
        report.kinetic_before == report.kinetic_after &&
        all(
            iszero,
            (
                report.particle_charge_change,
                report.background_charge_change,
                report.left_plate_change,
                report.right_plate_change,
            ),
        ) || throw(ArgumentError("dynamic composite mesh work/inventory report is inconsistent"))
    if report.selection === :manual
        report.threshold === nothing && report.buffer == 0 ||
            throw(ArgumentError("manual dynamic regrid has adaptive selection metadata"))
    elseif report.selection === :base_electric_curvature
        report.threshold isa T &&
            isfinite(report.threshold) &&
            _stored_nonnegative(report.threshold) &&
            report.buffer >= 0 &&
            all(spec -> spec[1] == 1, report.after_refinements) ||
            throw(ArgumentError("dynamic adaptive selection metadata is invalid"))
    else
        throw(ArgumentError("unknown dynamic composite regrid selection policy"))
    end
    return nothing
end

function _dce_phase(
    st::DynamicCompositeElectrostaticPIC{T},
    current,
    history,
    time,
    step,
    dt,
    epoch_time,
    epoch_step,
    regrids,
    initial,
    defect,
    work,
    report,
) where {T}
    time isa T &&
        step isa Int &&
        dt isa T &&
        epoch_time isa T &&
        epoch_step isa Int &&
        regrids isa Int &&
        regrids >= 0 || throw(ArgumentError("dynamic composite phase scalar types are invalid"))
    _dce_epoch_clock(
        time,
        step,
        dt,
        epoch_time,
        epoch_step,
        current.time[],
        current.step[],
        current.last_dt[],
    )
    all(_sce_q_finite, (initial, defect, work)) &&
        initial >= 0 &&
        _dce_energy(current) - initial == defect + work ||
        throw(ArgumentError("dynamic composite cumulative energy accounting is inconsistent"))
    iszero(step) &&
        !iszero(defect) &&
        throw(ArgumentError("unadvanced dynamic state cannot have physical-step energy defect"))
    if regrids == 0
        report === nothing &&
            iszero(work) &&
            _particle_stored_zero(epoch_time) &&
            epoch_step == 0 ||
            throw(ArgumentError("unregridded dynamic state has mesh work or an epoch reset"))
    else
        _dce_regrid_report(st, report, epoch_time, epoch_step, regrids, current.mesh)
        if current.step[] == 0
            report.field_after == _sce_field_energy(current.mesh, current.E) &&
                report.kinetic_after == _kinetic_energy_rational(_dce_private_species(current)) ||
                throw(
                    ArgumentError(
                        "dynamic regrid report does not match its current physical state",
                    ),
                )
            rhs = fill(zero(_CompositeQ), length(current.phi))
            _sce_source_rhs!(
                rhs,
                _dce_private_species(current),
                current.mesh,
                current.shape,
                current.background_charge,
                current.left_plate_charge,
                current.right_plate_charge[],
            )
            residuals = _sce_residuals(current.mesh, current.E, current.phi, rhs)
            report.max_gauss_residual == residuals.max_gauss &&
                report.max_constitutive_residual == residuals.max_constitutive ||
                throw(ArgumentError("dynamic regrid residual report does not match current fields"))
        end
    end
    if current.step[] == 0 && step > 0
        _dce_validate_history(history, st, current, time, step, dt)
    else
        history === nothing || throw(
            ArgumentError("dynamic prior-step snapshot is present outside a regridded empty epoch"),
        )
    end
    return nothing
end

function _dce_ready(st::DynamicCompositeElectrostaticPIC, species)
    _dce_storage(st, species)
    st.initialized || throw(ArgumentError("call init! before dynamic composite evolution or save"))
    _sce_ready(st._state, _dce_private_species(st._state))
    _dce_phase(
        st,
        st._state,
        st._history,
        st.time,
        st.step,
        st.last_dt,
        st.epoch_start_time,
        st.epoch_start_step,
        st.regrid_count,
        st.initial_energy,
        st.physical_defect,
        st.mesh_work,
        st.last_regrid_report,
    )
    return nothing
end

function _dce_publish!(
    st::DynamicCompositeElectrostaticPIC{T},
    current,
    owners::Vector{_ElectrostaticHostParticle{1,T}},
    components::Vector{Tuple},
    counts::Vector{Int},
    history,
    time,
    step,
    dt,
    epoch_time,
    epoch_step,
    regrids,
    initial,
    defect,
    work,
    report,
) where {T}
    # All numerical/graph/ledger work is complete. Only same-typed host stores follow.
    st._state = current
    st.mesh, st.E, st.phi = current.mesh, current.E, current.phi
    st._owners, st._components, st._counts = owners, components, counts
    st._history = history
    st.time, st.step, st.last_dt = time, step, dt
    st.epoch_start_time, st.epoch_start_step = epoch_time, epoch_step
    st.regrid_count = regrids
    st.initial_energy, st.physical_defect, st.mesh_work = initial, defect, work
    st.last_regrid_report = report
    st.initialized = true
    return st
end

"""
    DynamicCompositeElectrostaticPIC(mesh, species; shape, background_charge=0,
                                     left_plate_charge=0)

Completed-boundary adaptive CPU 1D physical-cloud electrostatic PIC. Canonical
static-composite KDK/Neumann epochs provide every numerical step; `regrid!` and
`adapt!` conserve the particle/source inventory and report actual mesh energy
exchange without resetting the global physical clock. See the family guide.
"""
function DynamicCompositeElectrostaticPIC(
    mesh::StaticAMRHierarchy1D{T},
    species::AbstractVector;
    shape::YeeTriangularCloud,
    background_charge = zero(T),
    left_plate_charge = zero(T),
) where {T}
    _static_amr_validate(mesh)
    owned_mesh = StaticAMRHierarchy1D(mesh.nbase, mesh.L, mesh.refinements)
    # Validate the ORIGINAL graph before copying can conceal particle aliases.
    template = StaticCompositeElectrostaticPIC(
        owned_mesh,
        species;
        shape,
        background_charge,
        left_plate_charge,
    )
    private = _dce_copy_species(species, T)
    current = StaticCompositeElectrostaticPIC(
        owned_mesh,
        private;
        shape,
        background_charge = template.background_charge,
        left_plate_charge = template.left_plate_charge,
    )
    owners = _ElectrostaticHostParticle{1,T}[ps for ps in species]
    components = Tuple[_sce_particle_arrays(ps) for ps in species]
    z = zero(_CompositeQ)
    st = DynamicCompositeElectrostaticPIC{
        T,
        typeof(shape),
        typeof(current.background_charge),
        typeof(current.left_plate_charge),
        typeof(current._factors),
    }(
        owned_mesh,
        current.E,
        current.phi,
        zero(T),
        0,
        zero(T),
        false,
        zero(T),
        0,
        0,
        z,
        z,
        z,
        nothing,
        shape,
        current.background_charge,
        current.left_plate_charge,
        current._factors,
        _dce_configuration(current),
        current,
        nothing,
        owners,
        components,
        Int[length(ps.weight) for ps in species],
        Threads.Atomic{Bool}(false),
    )
    _dce_storage(st, species)
    return st
end

function init!(st::DynamicCompositeElectrostaticPIC{T}, species::AbstractVector) where {T}
    _dce_acquire(st)
    try
        _sce_mode()
        _sce_population_storage(species, T)
        _sce_values(st.mesh, species, st._factors)
        # This disposable canonical constructor admits original source/storage
        # links before the private copies and all later publication preparation.
        StaticCompositeElectrostaticPIC(
            st.mesh,
            species;
            shape = st.shape,
            background_charge = st.background_charge,
            left_plate_charge = st.left_plate_charge,
        )
        private = _dce_copy_species(species, T)
        current = StaticCompositeElectrostaticPIC(
            st.mesh,
            private;
            shape = st.shape,
            background_charge = st.background_charge,
            left_plate_charge = st.left_plate_charge,
        )
        init!(current, private)
        owners = _ElectrostaticHostParticle{1,T}[ps for ps in species]
        components = Tuple[_sce_particle_arrays(ps) for ps in species]
        counts = Int[length(ps.weight) for ps in species]
        initial, z = _dce_energy(current), zero(_CompositeQ)
        _dce_phase(st, current, nothing, zero(T), 0, zero(T), zero(T), 0, 0, initial, z, z, nothing)
        _dce_prepared_roles(st, current, species, nothing)
        return _dce_publish!(
            st,
            current,
            owners,
            components,
            counts,
            nothing,
            zero(T),
            0,
            zero(T),
            zero(T),
            0,
            0,
            initial,
            z,
            z,
            nothing,
        )
    finally
        _dce_release(st)
    end
end

function step!(st::DynamicCompositeElectrostaticPIC{T}, species::AbstractVector, dt::Real) where {T}
    _dce_acquire(st)
    try
        _sce_mode()
        tau = _validated_nonnegative_dt(T, dt; name = "dynamic composite ES")
        _particle_stored_zero(tau) &&
            !iszero(dt) &&
            throw(ArgumentError("nonzero dynamic composite duration must not round to zero"))
        _dce_ready(st, species)
        _particle_stored_zero(tau) && return st
        iszero(st.step) ||
            isequal(tau, st.last_dt) ||
            throw(ArgumentError("dynamic composite duration is fixed after first success"))
        time, count = _validated_step_clock(st.time, st.step, tau)
        time > st.time || throw(ArgumentError("dynamic composite global time must advance"))
        current = _dce_clone_epoch(st)
        private = _dce_private_species(current)
        step!(current, private, tau)
        defect = st.physical_defect + current.last_report[].total_energy_defect
        _dce_phase(
            st,
            current,
            nothing,
            time,
            count,
            tau,
            st.epoch_start_time,
            st.epoch_start_step,
            st.regrid_count,
            st.initial_energy,
            defect,
            st.mesh_work,
            st.last_regrid_report,
        )
        _dce_prepared_roles(st, current, species, nothing)
        # Public particle identities/counts and all allocation/range checks passed.
        for (destination, prepared) in zip(species, private)
            copyto!(destination.x[1], prepared.x[1])
            foreach(copyto!, destination.v, prepared.v)
        end
        return _dce_publish!(
            st,
            current,
            st._owners,
            st._components,
            st._counts,
            nothing,
            time,
            count,
            tau,
            st.epoch_start_time,
            st.epoch_start_step,
            st.regrid_count,
            st.initial_energy,
            defect,
            st.mesh_work,
            st.last_regrid_report,
        )
    finally
        _dce_release(st)
    end
end

function _dce_regrid!(
    st::DynamicCompositeElectrostaticPIC{T},
    species,
    mesh,
    selection,
    threshold,
    buffer,
) where {T}
    # Both public callers hold the exclusive guard and have admitted the ready
    # state. Avoid repeating global/history/source validation in this private row.
    mesh isa StaticAMRHierarchy1D{T} ||
        throw(ArgumentError("dynamic regrid requires a matching-precision interval hierarchy"))
    _static_amr_validate(mesh)
    mesh.nbase == st.mesh.nbase && isequal(mesh.L, st.mesh.L) || throw(
        ArgumentError("dynamic regrid must preserve its base grid and original physical domain"),
    )
    isequal(mesh.refinements, st.mesh.refinements) && return st
    count = Base.checked_add(st.regrid_count, 1)
    owned = StaticAMRHierarchy1D(mesh.nbase, mesh.L, mesh.refinements)
    private = _dce_copy_species(_dce_private_species(st._state), T)
    current = StaticCompositeElectrostaticPIC(
        owned,
        private;
        shape = st.shape,
        background_charge = st.background_charge,
        left_plate_charge = st.left_plate_charge,
    )
    init!(current, private)
    oldprivate = _dce_private_species(st._state)
    oldfield, newfield = _sce_field_energy(st.mesh, st.E), _sce_field_energy(owned, current.E)
    oldkinetic, newkinetic = _kinetic_energy_rational(oldprivate), _kinetic_energy_rational(private)
    solution = _sce_ready(current, private)
    residuals = _sce_residuals(owned, current.E, current.phi, solution.rhs)
    exchange = newfield - oldfield
    report = DynamicCompositeESRegridReport{T}(
        st.time,
        st.step,
        count,
        st.mesh.refinements,
        owned.refinements,
        length(st.mesh.widths),
        length(owned.widths),
        oldfield,
        newfield,
        oldkinetic,
        newkinetic,
        exchange,
        _dce_particle_charge(private) - _dce_particle_charge(oldprivate),
        _exact_rational(st.background_charge) *
        (_exact_rational(owned.L) - _exact_rational(st.mesh.L)),
        _exact_rational(current.left_plate_charge) - _exact_rational(st.left_plate_charge),
        current.right_plate_charge[] - st._state.right_plate_charge[],
        residuals.max_gauss,
        residuals.max_constitutive,
        selection,
        threshold,
        buffer,
    )
    history = if st._state.step[] > 0
        (
            epoch_start_time = st.epoch_start_time,
            epoch_start_step = st.epoch_start_step,
            global_time = st.time,
            global_step = st.step,
            state = _checkpoint_state(st._state, oldprivate),
        )
    else
        st._history
    end
    work = st.mesh_work + exchange
    _dce_phase(
        st,
        current,
        history,
        st.time,
        st.step,
        st.last_dt,
        st.time,
        st.step,
        count,
        st.initial_energy,
        st.physical_defect,
        work,
        report,
    )
    _dce_prepared_roles(st, current, species, history)
    return _dce_publish!(
        st,
        current,
        st._owners,
        st._components,
        st._counts,
        history,
        st.time,
        st.step,
        st.last_dt,
        st.time,
        st.step,
        count,
        st.initial_energy,
        st.physical_defect,
        work,
        report,
    )
end

"Change the completed composite hierarchy without moving particles or global time."
function regrid!(st::DynamicCompositeElectrostaticPIC, species::AbstractVector, mesh)
    _dce_acquire(st)
    try
        _dce_ready(st, species)
        return _dce_regrid!(st, species, mesh, :manual, nothing, 0)
    finally
        _dce_release(st)
    end
end

"Adapt base-cell refinement using the canonical undivided base-electric-field curvature."
function adapt!(
    st::DynamicCompositeElectrostaticPIC{T},
    species::AbstractVector;
    threshold::T,
    buffer::Integer = 1,
) where {T}
    _dce_acquire(st)
    try
        _dce_ready(st, species)
        field = patch_electric_field(st._state, 1)
        specs = _dce_tagged_refinements(field, threshold, buffer)
        mesh = StaticAMRHierarchy1D(st.mesh.nbase, st.mesh.L, specs)
        return _dce_regrid!(st, species, mesh, :base_electric_curvature, threshold, Int(buffer))
    finally
        _dce_release(st)
    end
end
