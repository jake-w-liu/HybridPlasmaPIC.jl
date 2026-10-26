# Prescribed currents on a 1D Yee grid

`step_maxwell!(model, dt, impulse)` advances the existing synchronized Yee field
state with a caller-supplied current impulse, `impulse = integral(J dt)`, in units
with epsilon0=1. It supports CPU Float32/Float64 in one spatial dimension, with
periodic or homogeneous perfectly conducting boundaries. The two-argument call
continues to advance vacuum fields in one to three dimensions.

`impulse` is a tuple of three ordinary vectors with the shapes and precision of
`model.fields.E`. Arrays must be finite and independent of one another and all
physical/staging model arrays. Conducting boundaries require zero tangential
current impulse at both endpoints; incompatible data is rejected rather than
projected, before any staging write. Inputs remain unchanged. Use the electric
coordinates returned by `yee_coordinates` when sampling a source, and account for
the time integral in the supplied values. The solver does not choose a source
quadrature for the user.

The magnetic half-step, sourced electric step and final magnetic half-step use
the same native stages as vacuum stepping. The electric stage is
`(Eold + dt*c²*curl(Bhalf)) - impulse`, rounded once (in Float64 when a Float32
multiplier is outside the normal range); an entry whose native value overflows
is evaluated exactly with the impulse included. Each stage meets the
[normwise contract](BOUNDED_MAXWELL_IMPLEMENTATION.md#numerical-contract-normwise)
with `|impulse|` added to its scale; a zero impulse entry, including a negative
zero, reproduces the vacuum value bitwise. Nonfinite complete fields, invalid
inputs, CFL or clock failures reject the step without changing physical fields,
their array identities, time or count. Positive dt must be strictly below
`maxwell_cfl(model)`; arbitrary variable-step sequences have no stability guarantee.

The diagnostics expose three different quantities:

| Operation | Meaning |
| --- | --- |
| `yee_electric_charge!(out, E, grid)` | Extensive charge on every 1D node: PEC surface charges included, periodic incidence wrapped; no division by dx. |
| `yee_source_work(Eold, Enew, impulse, grid)` | `dx*sum(impulse*(Eold+Enew))/2`, positive for energy transferred out of fields. Caller supplies both electric states. |
| `yee_modified_energy(model, dt)` | Fixed-step quadratic `dx/2*(sum(E²)+c²*sum(B²)-dt²*c²*sum(curl(E)²)/4)`. |

All-node charge is distinct from the existing interior-only electric
`yee_divergence!`. Its PEC entries are `[Ex[1], diff(Ex)..., -Ex[end]]` with unit
transverse area; each entry is one correctly rounded difference and their exact
sum, including both surfaces, is zero. The exact update is
`Qnew-Qold = -incidence(impulse)`. No initial Gauss projection, particle charge
deposition or automatic wall charge model is supplied.

At fixed dt, exact arithmetic gives `Hnew-Hold = -source_work`; rounded stage
results add numerical residuals. This follows by writing the two magnetic states
around their common half-stage B: the curl/transpose work terms cancel the
electric update, leaving `-dx*sum(impulse*(Eold+Enew))/2`. The negative curl term
has magnitude at most `(dt*c/dx)²` times the electric norm in 1D, so the quadratic
is positive definite strictly inside the CFL bound. Changing dt changes H itself.
Physical `em_field_energy` and physical Poynting diagnostics retain their meanings;
they are not this modified-energy invariant. Work and modified energy use
pairwise summation (Float64 for Float32 data), agree with their exact values to a
few `eps` per summed term times the sum of absolute contributions, complete
overflow or non-negligible underflow exactly and return signed infinity if the
exact result is unrepresentable.

See [`examples/sourced_maxwell_cavity.jl`](../examples/sourced_maxwell_cavity.jl)
for a driven cavity and explicit work balance; verification is recorded in the
[implementation record](BOUNDED_MAXWELL_IMPLEMENTATION.md#prescribed-current-slice).

This capability is prescribed-source field evolution. Coupled particle motion,
deposition/gather, particle absorption, wall-current ledgers and nonlinear event
solves remain incomplete. The retained fixed-macrocharge self-force and coupled
impact-velocity convergence failures are still open; this prerequisite does not
resolve them. Higher-dimensional impressed sources, open/PML boundaries, circuits,
particle injection, reflecting particles, GPU/MPI and restart remain later slices.
