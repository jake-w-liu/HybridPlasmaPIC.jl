using HybridPlasmaPIC

function sourced_cavity_example()
    # Prescribed current in a conducting 1D cavity, epsilon0=1 and unit transverse area.
    # This example has no particles: current is supplied by the caller.
    grid = YeeGrid((64,), (1.0,))
    model = MaxwellVacuum(grid)
    dt = 0.25 * maxwell_cfl(model)
    current_impulse = map(a -> zero(a), model.fields.E)
    x = only(yee_coordinates(grid, :electric, 2))
    current_impulse[2] .= dt * 0.2 .* sinpi.(x)
    current_impulse[2][1] = current_impulse[2][end] = 0.0
    old_E = map(copy, model.fields.E)
    initial_modified = yee_modified_energy(model, dt)
    work = 0.0
    for _ = 1:128
        foreach(copyto!, old_E, model.fields.E)
        step_maxwell!(model, dt, current_impulse)
        work += yee_source_work(old_E, model.fields.E, current_impulse, grid)
    end
    charge = zeros(grid.n[1] + 1)
    yee_electric_charge!(charge, model.fields.E, grid)
    return (
        time = model.time[],
        physical_energy = em_field_energy(model),
        modified_energy = yee_modified_energy(model, dt),
        source_work = work,
        modified_balance_residual = yee_modified_energy(model, dt) - initial_modified + work,
        all_node_charge = sum(charge),
    )
end

println(sourced_cavity_example())
