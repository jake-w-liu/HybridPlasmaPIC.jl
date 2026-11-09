"""
    composite_electrostatic_diagnostics(state::DynamicCompositeElectrostaticPIC, species)

Canonical current-mesh charge/field diagnostics plus global/epoch phase, exact
physical-step defect and mesh-energy exchange. The preceding physical report
retains its actual old mesh after regridding; returned data are detached.
"""
function composite_electrostatic_diagnostics(
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector,
)
    _dce_acquire(st)
    try
        _dce_ready(st, species)
        current = composite_electrostatic_diagnostics(st._state, _dce_private_species(st._state))
        history = if st._state.step[] > 0
            (
                epoch_start_time = st.epoch_start_time,
                epoch_start_step = st.epoch_start_step,
                global_time = st.time,
                global_step = st.step,
                state = _checkpoint_state(st._state, _dce_private_species(st._state)),
            )
        else
            deepcopy(st._history)
        end
        return merge(
            current,
            (
                time = st.time,
                step = st.step,
                last_dt = st.last_dt,
                epoch_start_time = st.epoch_start_time,
                epoch_start_step = st.epoch_start_step,
                epoch_time = st._state.time[],
                epoch_step = st._state.step[],
                mesh = deepcopy(st.mesh),
                E = copy(st.E),
                phi = copy(st.phi),
                regrid_count = st.regrid_count,
                initial_energy = deepcopy(st.initial_energy),
                cumulative_physical_defect = deepcopy(st.physical_defect),
                cumulative_mesh_work = deepcopy(st.mesh_work),
                accounting_defect = current.total_energy - st.initial_energy - st.physical_defect -
                                    st.mesh_work,
                last_regrid_report = deepcopy(st.last_regrid_report),
                last_physical_step = history,
            ),
        )
    finally
        _dce_release(st)
    end
end

"Return the canonical width-restricted covered/leaf electric field for one current patch."
function patch_electric_field(st::DynamicCompositeElectrostaticPIC, patch_id::Integer)
    _dce_acquire(st)
    try
        _dce_ready(st, st._owners)
        return patch_electric_field(st._state, patch_id)
    finally
        _dce_release(st)
    end
end
