using HybridPlasmaPIC

# Normalized c = 1, homogeneous PEC walls at x = 0 and x = 1. Both E and B
# are stored at the same time; no manual half-step initialization is needed.
grid = YeeGrid((128,), (1.0,))
model = MaxwellVacuum(grid)
x = only(yee_coordinates(grid, :electric, 2))
model.fields.E[2] .= sinpi.(x)
apply_yee_pec!(model.fields, grid)
initial_energy = em_field_energy(model)
dt = 0.5 * maxwell_cfl(model)
for _ = 1:256
    step_maxwell!(model, dt)
end
divB = zeros(grid.n)
yee_divergence!(divB, model.fields.B, grid)
println((
    time = model.time[],
    initial_energy,
    energy = em_field_energy(model),
    max_divB = maximum(abs, divB),
    outward_flux = yee_poynting_flux(model),
))
