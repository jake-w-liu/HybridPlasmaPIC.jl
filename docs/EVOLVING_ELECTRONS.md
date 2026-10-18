# Evolving electron pressure

`EvolvingElectrons(γ = 5/3)` gives periodic `HybridStepper` runs an independent
scalar electron pressure. It supports one proton population, ordered ion
collections, and replicated MPI. Set `stepper.fields.pe` before `init!`;
initialization preserves it. The default zero pressure is a valid cold state.
For a single species with a different charge, use a one-element ion collection
to obtain charge-weighted moments; the scalar interface retains proton normalization.

The closure requires finite `γ > 1` and `η = ηH = de2 = 0`, including after
conversion to the grid precision. CAM-CL and Hall-MHD reject this closure.
Heating, heat flux, and anisotropic pressure evolution are separate model
extensions.

## Initialize and advance

```julia
using HybridPlasmaPIC

g = FourierGrid((32,), (2pi,))
ions = ParticleSet{1,Float64}(128)
load_lattice_1d!(ions, 0.0, g.L[1])
set_density_weight!(ions, 1.0, g)

closure = EvolvingElectrons(5 / 3)
st = HybridStepper(g, HybridModel(closure), CIC(), nparticles(ions))
x = (0:g.n[1]-1) .* g.dx[1]
st.fields.pe .= 0.7 .+ 0.04 .* cos.(x)
init!(st, ions)
step!(st, ions, 0.001; NB = 2)

budget = energy_budget(ions, st.fields.B, st.fields.n, closure, g;
                       pe = st.fields.pe)
```

The explicit timestep in this example is for this initial condition. Choose
and refine the timestep for the actual pressure, flow, and magnetic scales.
`recommended_dt` estimates magnetic and gyro limits; it does not bound
electron-pressure advection, compression, or ion-acoustic accuracy.

The equation is

```
∂p_e/∂t = −u_e·∇p_e − γ p_e ∇·u_e,
u_e = u_i − J/max(n, nfloor).
```

Each magnetic subcycle advances B and pressure together through RK4. Every
stage uses its own B, current, and pressure, with charge-weighted ion moments
held at the particle-step midpoint. The outer particle/field coupling is
second order; RK4's fourth order applies to the subsystem with fixed moments.
Changing `dt` or `NB` preserves the pressure state.

After a step, B and pressure are at integer time; particle velocities retain
their half-time level. Recomputing carried E reads the advanced pressure
without evolving it again. Direct `ohms_law!` also reads pressure without
changing it. `electron_pressure!(pe, n, closure)` rejects this closure because
density alone cannot reconstruct the stored pressure.

## Rejected states and retries

Pressure must be finite and nonnegative on entry, at every RK4 stage, and in
each completed subcycle. Zero is accepted. Any negative value is rejected,
including a small negative roundoff value. There is no clipping or pressure
floor. The explicit spectral method does not preserve positivity for all
data and timesteps; a smaller timestep or better-resolved positive profile
may be necessary.

A rejected serial `init!` or `step!` restores public field arrays, particle
positions and velocities, the density-floor counter, and simulation clocks.
Array identities are preserved, so existing views remain valid. Private
workspaces remain scratch. A smaller-step retry can use the restored state.
A zero timestep returns without validating or changing pressure.

Pressure must have independent storage. Public Ohm-law, stepping, and
persistence entry points reject pressure that overlaps their other state or
work arrays before mutation.

MPI ranks must start with identical replicated fields, including pressure,
and matching grid/model configuration. Pressure validation failures are
coordinated across ranks before further moment reductions or migration; the
numerical transaction restores each rank's local public state. MPI transport,
process, and migration failures are outside this guarantee.

## Restart and energy diagnostics

`save_checkpoint`, `load_checkpoint!`, `archive_run`, and distributed
checkpoints preserve pressure along with the existing field and particle
state. A saved evolving-pressure state must contain valid pressure; restore
does not reconstruct it from density or call `init!`. Use the same numerical
environment and stepping inputs for bitwise continuation.

Pass `pe = st.fields.pe` to `electron_internal_energy` or `energy_budget`.
The electron reservoir is `∫ p_e/(γ−1) dV`. Omitting `pe` for this closure is
an error; algebraic closures reject an alternate `pe` argument. Collections
share one electron reservoir.
The stored-pressure integral preserves the ordinary floating-point evaluation
order and uses exact arithmetic when intermediate range loss or possible final
overflow requires it. This protection applies to the electron-pressure term;
[other diagnostic range defects](HYBRID_DIAGNOSTIC_RANGE_PLAN.md) remain tracked.

Energy diagnostics use the supplied particle velocities without a hidden
time-centering kick. Compare reservoirs at a common time when studying energy
convergence. Ideal, periodic, floor-inactive continuum energy conservation
does not imply exact conservation by the explicit integrator. See the
[energy-model balance](HYBRID_ENERGY_MODEL_PLAN.md) for the scope of that
identity and the separate dissipative and inertial model limitations.

## Workspace cost

Only this closure allocates pressure stages and transaction buffers. For D
dimensions, Ncell grid cells and Nparticle local particles, their numeric
storage is `(27 + D)*Ncell + (D + 3)*Nparticle` grid-precision values, plus
array headers and bookkeeping. Particle buffers resize with local particle
counts. Successful steps copy public state into the saved buffers; rejected
steps also copy it back. Algebraic closures do not allocate these buffers.
