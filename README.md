# HybridPlasmaPIC.jl

A research-grade **hybrid particle-in-cell (PIC) plasma solver** in Julia:
kinetic ions (Boris-pushed macroparticles), massless fluid electrons
(generalized Ohm's law), and spectral fields on periodic grids — plus a family
of companion solvers (Hall-MHD, electrostatic PIC, full electromagnetic PIC),
collision operators, shock drivers, geometric-optics ray tracing, and an
extensive analytic-oracle verification suite.

The periodic hybrid core supports 1D3V, 2D3V, and 3D3V: positions carry `D`
coordinates and velocities carry 3. Particle, field and integrator state lives
in the `ParticleSet`, `HybridFields` and `HybridStepper` objects you own.
Companion models have the dimensional and physical limits described below.

**Contents:**
[Install](#installation) ·
[Quickstart](#quickstart) ·
[Physics & units](#the-physics-model) ·
[Time integration](#time-integration) ·
[Other models](#beyond-the-hybrid-model) ·
[Collisions](#collisions-and-extra-particle-physics) ·
[Shocks](#shocks) ·
[Diagnostics](#diagnostics) ·
[Ray tracing](#ray-tracing) ·
[Boundaries, parallelism, I/O](#boundaries-parallelism-and-io) ·
[Verification](#verification) ·
[Examples](#examples) ·
[Layout](#package-layout) ·
[References](#references)

## Installation

Requires **Julia 1.12 or later**. `SpectralOperators.jl` (the sibling package
providing the Fourier/SBP operators) is in the General registry, so it resolves
automatically as a dependency:

```julia
using Pkg
Pkg.add(url = "https://github.com/jake-w-liu/HybridPlasmaPIC.jl.git")
```

## Quickstart

An ion-acoustic wave in the hybrid model — load a quiet lattice of ions, seed a
small velocity perturbation, and step (runs in ~20–25 s including compilation;
`examples/ion_acoustic.jl` is the full version that also measures ω against
k·c_s):

```julia
using HybridPlasmaPIC, Random

g  = FourierGrid((64,), (2π,))                 # 64 cells, box length 2π d_i
ps = ParticleSet{1,Float64}(64 * 400)          # 400 particles per cell
load_lattice_1d!(ps, 0.0, 2π)                  # quiet spatial load
set_density_weight!(ps, 1.0, g)                # weights so that n = n0 = 1
load_quiet_velocities!(ps, MersenneTwister(1), (0, 0, 0), (0, 0, 0))  # cold
for p in eachindex(ps.weight)                  # seed a small acoustic wave
    ps.v[1][p] += 0.01 * sin(ps.x[1][p])
end

st = HybridStepper(g, HybridModel(IsothermalElectrons(1.0)), CIC(), nparticles(ps))
init!(st, ps)
for _ = 1:1000
    step!(st, ps, 0.02)
end

b = energy_budget(ps, st.fields.B, st.fields.n, IsothermalElectrons(1.0), g)
println(b)   # (kinetic = ..., magnetic = ..., electron_internal = ..., electron_inertia = ..., total = ...)
```

Two things every setup needs:

- `set_density_weight!(ps, 1.0, g)`, so a uniform load deposits `n = n0 = 1`
  — otherwise every `1/n` term in Ohm's law is silently rescaled;
- a timestep from `recommended_dt(g; NB, integrator = :rk4)` (see
  [Time integration](#time-integration)) unless you know your whistler CFL.

For several positive ion populations, pass a vector of `ParticleSet`s to the
same initialization and stepping functions: `init!(st, species)` and
`step!(st, species, dt; NB)`, or the CAM-CL equivalents. Each population has its
own charge, mass, weights and particle count. They share one field update with
`n_e = sum(q_s*n_s)` and charge-weighted ion velocity; the density floor is
applied after summation. The single-population interface retains its proton
normalization. See [the proton/alpha example](examples/two_ion_hybrid.jl).

Populations need independent storage and a common velocity time level. Empty
populations and changing particle counts are supported; reinitialize when
adding particles loaded at a new physical velocity time level. Collection
stepping currently uses the serial periodic solvers. Its explicit predictor
has second-order momentum error in the uniform two-ion convergence test;
discrete momentum is not conserved exactly. Use `kinetic_energy(species)`,
`total_momentum(species)`, and the collection forms of `energy_budget` and
`momentum_budget` for totals; shared field energy is counted once.

## The physics model

The core solver integrates the standard massless-electron hybrid model in
**Ω_ci-normalized units** (proton `q/m = 1`, `μ0 = 1`):

```
ions:      dx/dt = v,   dv/dt = E + v×B          (Boris macroparticles)
current:   J   = ∇×B
electrons: u_e = u_i − J/n
Ohm's law: E   = −u_i×B + (J×B)/n − ∇p_e/n + ηJ − ηH∇²J
Faraday:   ∂B/∂t = −∇×E,    ∇·B = 0  (spectrally enforced)
```

| quantity | unit |
|---|---|
| length | ion inertial length `d_i = c/ω_pi` |
| time | inverse ion gyrofrequency `Ω_ci⁻¹` |
| velocity | Alfvén speed `v_A` |
| B, E | `B0`, `v_A B0` |
| density | `n0` |
| current | `B0/(μ0 d_i)` |
| pressure | `B0²/μ0` |

`PlasmaUnits(; n0, B0, mi)` plus `to_SI` / `to_normalized` convert the table's
quantities between this system and SI; `alfven_speed`, `gyrofrequency`, and
`inertial_length` give the reference scales. All six stored references, including
physical constants supplied to a positional constructor, must be finite and
positive. Ordinary conversions keep their arithmetic order. Range-sensitive
conversions evaluate the complete expression before rounding, recovering
representable results even when the scale itself cannot be represented.
`BigFloat` arithmetic uses the current working precision and rounding mode;
typed construction retains the precision of existing `BigFloat` references.

CUDA/Metal array conversion inspects ordinary results on the device and
synchronizes. Exceptional cases use host arithmetic, then return the original
backend and shape; zero-dimensional arrays return a scalar. Native Metal Float16/32
arrays and reference scales have been tested; native NVIDIA execution remains
unverified. See
[NORM-001](docs/NORMALIZATION_RANGE_PLAN.md) for coverage, costs and open gates.

### Electron closures

- `IsothermalElectrons(Te)` — `p_e = Te·n`
- `PolytropicElectrons(pe0, n0, γ)` — `p_e = pe0·(n/n0)^γ`
- `EvolvingElectrons(γ = 5/3)` — independent adiabatic scalar pressure,
  advanced with B by `HybridStepper`. Initialize `fields.pe` before `init!`;
  requires `γ > 1` and `η = ηH = de2 = 0`. Supports ion collections and MPI.
- `CGLElectrons(p⊥0, p∥0, n0, B0)` — anisotropic double-adiabatic
  (Chew–Goldberger–Low): `p_⊥ ∝ n·B`, `p_∥ ∝ n³/B²`, entering Ohm's law
  through the full gyrotropic pressure-tensor force `−(∇·P_e)/n`. (Note the
  fluid mirror mode has no short-wavelength cutoff, so strongly
  mirror-unstable CGL states blow up — a property of the CGL model itself.)

Evolving pressure must remain finite and nonnegative at every RK4 stage;
rejected steps restore public state without clipping pressure. See the
[initialization, timestep, restart, and energy guide](docs/EVOLVING_ELECTRONS.md).
The standalone `advance_electron_pressure!` remains a forward-Euler building
block; it does not provide the coupled stepper's validation or rollback.

### Optional physics (off by default; `nfloor` is an always-on guard)

`HybridModel(closure; η = 0.0, ηH = 0.0, de2 = 0.0, nfloor = 1e-6)`:

- `η` — resistivity (`+ηJ` in Ohm's law),
- `ηH` — hyper-resistivity (`−ηH∇²J`),
- `de2` — electron inertia `d_e²`: keeps the leading finite-electron-mass term
  by filtering the **transverse** part of `E` with `1/(1 + d_e²k²)`, capping
  the whistler frequency at the `d_e` scale (the longitudinal/electrostatic
  field is untouched, as the physics requires); `de2 = 0` recovers the
  massless model exactly,
- `nfloor` — density floor protecting the `1/n` divisions.

## Time integration

Two production steppers advance the hybrid model, both explicit and 2nd-order,
both subcycling `B` through the stiff whistler branch with frozen midpoint ion
moments:

- **`HybridStepper`** (`init!`, `step!(st, ps, dt; NB)`) — the default.
  Boris leapfrog for particles; `B` advanced `n → n+1` by `NB` RK4 substeps;
  the carried `E` recomputed each step. Supports every closure, including CGL;
  `EvolvingElectrons` advances pressure through the same RK4 stages as B.
- **`CAMCLStepper`** (`init_camcl!`, `step_camcl!(st, ps, dt; NB ≥ 2)`) —
  Matthews' CAM-CL (J. Comput. Phys. 112, 102, 1994): current advance method
  plus a cyclic-leapfrog `B` subcycle whose two staggered field copies persist
  across particle steps (neutrally stable, non-dissipative, for whistler
  `ω·h < 2` per substep; two Ohm evaluations per substep instead of RK4's
  four). Algebraic scalar closures only.

**Choosing dt.** `recommended_dt(g; NB = 1, integrator = :rk4, safety = 0.8,
η = 0.0, ηH = 0.0)` returns a conservative particle step from the whistler CFL
evaluated at the spectral-corner wavenumber `kmax = √Σ_d (π/Δx_d)²` (the
stiffest representable mode in 2D/3D — a single-axis Nyquist estimate is up to
`D×` too optimistic), additionally capped by Boris gyro-accuracy
(`Ω_ci·Δt ≤ 0.3`) and, when `η`/`ηH` are supplied, by the resistive real-axis
stability of the chosen integrator. In 1D, in the fine-grid limit
`kmax·d_i ≫ 1`, it reduces to the familiar `Ω_ci Δt ≲ (1/π)(Δx/d_i)²`.
These limits do not bound evolving-pressure advection/compression or
ion-acoustic accuracy; choose and refine `dt` to resolve those dynamics too.

`run_whistler` is a separate linear Fourier reference for Crank–Nicolson versus
forward Euler. Its `cn_multiplier` preserves finite phase components across
extreme products, including near quadrature; the driver combines the complete
frequency before exceptional rounding. Frequencies beyond finite Float64 range
raise `ArgumentError`. Rational-only multiplier calls remain exact and widen
the integer type when needed. See the [linear CN contract](docs/LINEAR_CN_RANGE_PLAN.md)
for compatibility and validation. Hardware floating timesteps with fixed-width
step counts reuse each mode's multiplier while preserving its ordered history;
the [packet allocation plan](docs/WHISTLER_PACKET_ALLOCATION_PLAN.md) records the
optimization and acceptance gates. Nonlinear coupled implicit evolution remains
in the [implementation plan](docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md).

## Beyond the hybrid model

| model | what it is | entry points |
|---|---|---|
| **Hall-MHD** | periodic single-fluid Hall-MHD with algebraic electron closures; ion pressure `p_i = Ti·n`; the electron pressure force acts in both Ohm's law and the bulk momentum equation; explicit RK4 with divergence-free projection | `HallMHDModel`, `HallMHDState`, `step_hall_mhd!` |
| **Electrostatic PIC** | arbitrary kinetic species on periodic CPU grids, explicit charge/background policy, shared spectral Poisson field and synchronized leapfrog in 1D/2D/3D; scalar electron compatibility | `ElectrostaticPIC` (or `Electrostatic1D`), `init_espic!`, `step_espic!` |
| **Electromagnetic PIC** | full EM PIC with kinetic electrons, optional mobile ions (`mobile = true`), optional relativistic Boris push (`relativistic = true`), spectral current correction for charge continuity, Gauss-law initialization, electron/Maxwell subcycling (`n_sub`) in 1D–3D | `EMPIC`, `EMPIC1D`, `init_empic!`, `step_empic!`, `charge_conservation_residual` |
| **Bounded Maxwell vacuum** | staggered Yee fields in 1D–3D with periodic or homogeneous conducting axes, synchronized second-order stepping, CFL rejection, compatible divergence and energy/Poynting diagnostics; native stencils under the normwise accuracy contract; no particle coupling or absorbing layer | `YeeGrid`, `MaxwellVacuum`, `step_maxwell!`, `yee_divergence!`, `yee_poynting_flux` |
| **Prescribed-current Yee fields** | 1D current impulses on electric lattices, all-node charge including conducting surfaces, source work and fixed-step modified energy; [scope and example](docs/SOURCED_YEE_FIELDS.md); particle coupling remains open | `step_maxwell!(model, dt, impulse)`, `yee_electric_charge!`, `yee_source_work`, `yee_modified_energy` |
| **Fusion geometry / gyrokinetics / AMR** | integrated tokamak-oriented building blocks: toroidal curvilinear metric operators, guiding-centre gyrokinetics, and finite-volume AMR primitives alongside the existing mesh/model subsystems | `ToroidalGrid`, `metric_gradient`, `GuidingCentre`, `push_guiding_centre!`, `AMRGrid`, `refine` |
| **KdV reference** | integrating-factor pseudo-spectral KdV solver + analytic soliton, used as a nonlinear-wave verification target | `kdv_solve`, `kdv_soliton` |

The bounded vacuum model uses a constant-step Maxwell stability bound.
Individually sub-CFL steps do not guarantee stability for an arbitrary variable
step sequence; see its [centering and stability contract](docs/BOUNDED_MAXWELL_IMPLEMENTATION.md#discrete-contract-and-derivation).

The PIC models share the grid, particle, and diagnostic infrastructure with
the hybrid stepper. Their species and time-level contracts differ;
[the electrostatic collection example](examples/electrostatic_species.jl) shows a
neutral pair plasma. Hall-MHD shares the grid and field machinery, and the KdV
solver is a standalone verification tool.

For an initialized ES collection, call `poisson_E!(state)` after manually changing
positions, weights or charges, or adding particles, to refresh the self-field
without resetting the velocity phase. New velocities must use the existing
shared half level. Replacing or reordering populations instead requires
`init_espic!` and physical full-time velocities.

General ES/EM longitudinal projections preserve finite results when intermediate
wavenumber squares or products lose range. The EM continuity diagnostic also
retains range-sensitive density differences and timestep-weighted currents through
normalization. Its existing spectral input restrictions and FFT accuracy limits
still apply, particularly when large terms nearly cancel. See the
[diagnostic scope](docs/CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md) and remaining
[current-correction and ES/Gauss range limits](docs/EMPIC_RANGE_PLAN.md).

For finite inputs and charge snapshots, EM PIC with built-in NGP/CIC/TSC shapes
applies a uniform correction based on actual signed particle weights, background
charge and the complete product of stored cell widths. The final total differs
from that target only by final cell-rounding errors; nonzero net charge is
retained. Local stencil and accumulation errors retain their separate arithmetic
limits. Custom shapes retain their represented coefficients. See the
[charge-capture contract and measured costs](docs/CHARGE_DEPOSITION_CONSERVATION_PLAN.md).

Relativistic stepping and restart require subluminal stored velocities. Boundary
validation compares the original components exactly when rounding could hide
`|v| ≥ c`; inputs whose rounded speed ratio is already at least one remain
unsupported. [Validation scope and evidence](docs/RELATIVISTIC_VALIDATION_PLAN.md).

## Collisions and extra particle physics

- `collide_coulomb!(ps, gcoeff, dt)` — Takizuka–Abe (1977) binary Coulomb
  scattering, with the Higginson–Holod–Link (JCP 2020) correction for unequal
  macro-particle weights (relaxation rates independent of the weight
  distribution; conservation exact for equal weights, in expectation
  otherwise), an isotropic large-angle fallback when the sampled scattering
  variance exceeds unity, and the T&A triplet for odd particle counts.
  Exact identity rotations preserve original velocity bits and RNG draw order.
  Binary-scaled rates and direct velocity increments retain small physical
  kicks; updates follow the normwise standard — widened arithmetic replays
  only range loss, while near-cancelling components keep fast rounding within
  a few ulps of the pair scale.
  Finite Gaussian draws use the scaled tangent when boundary roundoff would
  overflow its intermediate value; nonfinite draws still reject the call.
  See the [arithmetic validation scope](docs/COULOMB_SMALL_ANGLE_PLAN.md).
  The [population mean](docs/COULOMB_WEIGHT_MEAN_PLAN.md) retains
  large counts and sparse weights for hardware-float populations; BigFloat
  accumulation uses local guard bits without changing global precision.
- `collide_coulomb!(species, grid, dt; coulomb_log, epsilon0, weight_to_number)`
  — cell-local self and interspecies Coulomb scattering for ordered CPU
  populations on a periodic `FourierGrid` (Float32/Float64, 1D–3D). Pairs form
  within each cell from the physical charges, masses, macro-weights and cell
  volume; unequal cross counts correct each reused particle's rate, and odd
  counts use the half-strength triplet. Reuse a `CoulombCollisionWorkspace` to
  stay allocation-free. `coulomb_normalization(units)` supplies `epsilon0` and
  `weight_to_number` for the package's normalized units. Velocities stage until
  the whole call succeeds. [User guide](docs/LOCAL_COULOMB.md) and
  [two-species example](examples/local_coulomb.jl).
- `collide_bgk!(ps, ν, dt)` — BGK relaxation toward the set's own drifting
  Maxwellian, conserving momentum and energy exactly.
- `collide_neutral_mcc!` / `ionize_mcc!` — Monte-Carlo elastic scattering off
  a neutral Maxwellian reservoir, and electron-impact ionization with exact
  energy-threshold bookkeeping and unique newborn particle ids.
- `apply_radiation_reaction!(ps, E, B, dt; K, c)` — the dominant `γ²` term of
  the reduced Landau–Lifshitz drag, using an exponential substep with the initial
  drag rate held fixed. Float16/32/64 particles are scaled by a factor in [0, 1],
  so no stored speed increases, and stay normwise accurate across extreme finite
  light speeds, fields and coefficient–timestep products through guarded wider
  evaluation. Inputs must remain finite in particle precision;
  [arithmetic and accuracy record](docs/RADIATION_REACTION_PLAN.md).

## Shocks

A self-contained perpendicular collisionless-shock laboratory:

- **1D**: `PerpShock` — reflecting-wall piston shock on an SBP grid
  (non-periodic), with polytropic or Leroy-1982 electron-energy closures;
  `LeroyBoundary` + `step_leroy_shock!` run the wall-less shock-rest-frame
  variant with flux-weighted particle recycling at both ends.
  `run_perp_shock`, `run_perp_shock_rh`, `run_perp_shock_leroy` are complete
  drivers returning compression, reflected-ion fractions, and the independent
  Rankine–Hugoniot prediction (the sustained rh/leroy drivers additionally
  return the magnetic overshoot).
- **2D/3D**: `PerpShock2D` / `PerpShock3D` (mixed SBP-x × Fourier-transverse
  grids), rippled-front surface diagnostics, initial upstream
  Alfvénic-fluctuation seeding (3D, `db_turb`), and a campaign layer
  (`shock_campaign_3d`, `production_3d_case`) with bitwise
  checkpoint/restart. For `PerpShock3D`, `field_method=:cn` solves the frozen-moment
  induction equation with a checked true residual. Reuse
  `ShockCNWorkspace(sh, ps)` through `cn_workspace=`; configure `cn_rtol`,
  `cn_atol` and `cn_maxiter` for the solve. CN supports CPU Float32/Float64 with
  ordinary `Vector` particle storage. A failed solve leaves the entire caller
  state unchanged and reports `ShockCNConvergenceError`; nonfinite arithmetic
  also rejects the step. The default remains RK4. This magnetic solve retains
  the reduced induction model and explicit particle time-resolution requirements.
  See the [CN contract and measured costs](docs/SHOCK_CN_SOLVE_PLAN.md).
- **Jump conditions**: `rankine_hugoniot` / `rh_branches` solve the coplanar
  MHD Rankine–Hugoniot system with residuals of all six conservation laws
  reported; the field-aligned (`Bt = 0`) bifurcation is handled explicitly —
  inside the switch-on window the solver constructs the switch-on branch in
  closed form and selects the evolutionary solution.
- **Boundary machinery**: flux-correct open-boundary injection
  (`ShockInjector`, `flux_speed`, `inject_face_1d!`) samples inflow velocities
  from the drifting-Maxwellian flux distribution `p(v) ∝ v·f_M`.

## Diagnostics

- **Energy**: `energy_budget(ps, B, n, closure, g; de2 = 0, pe = nothing)` returns kinetic,
  magnetic, electron-internal/free-energy, and optional quadratic
  electron-inertia reservoirs. Periodic ideal massless scalar-closure
  evolution conserves the continuum total when density floors are inactive;
  the explicit integrator is checked through energy convergence. CGL has an
  additional anisotropic battery term. The reduced inertia reservoir closes
  the uniform-density fixed-ion field balance; moving-ion and dissipative energy accounting
  have further terms described in the [model balance](docs/HYBRID_ENERGY_MODEL_PLAN.md).
  For `EvolvingElectrons`, supply `pe = st.fields.pe` explicitly. Budgets use
  the supplied particle velocities, which are staggered in hybrid runs;
  energy-convergence comparisons need reservoirs at a common time.
  `kinetic_energy_relativistic(ps, c)` covers relativistic EMPIC runs.
  `resistive_dissipation` reports the Ohmic quadratic form;
  `jdotE_density`, `electric_work`, and `momentum_budget` report work and momentum.
- **Velocity space**: `pressure_tensor!` + `temperatures_par_perp` (T∥/T⊥
  relative to the local field), `pressure_strain` (the −P:∇u interaction),
  `velocity_histogram`, `phase_space_histogram`.
- **Waves & turbulence**: `power_spectrum`, `mode_amplitude`.
- **Shock analysis**: shock-front tracking, `boundary_energy_flux`,
  `CrossingLogger` (per-particle crossing counts and energy gain),
  reflected-ion classification (`classify_reflected`) plus the de
  Hoffmann–Teller frame transform (`dehoffmann_teller_velocity`), and
  synthetic spacecraft (`SyntheticProbe`, four-spacecraft timing) for
  comparing runs against heliospheric data conventions.
- **Housekeeping**: particle sorting (`sort_particles!`), occupancy and
  load-balance planning (`particle_load_balance`), provenance logging.

## Ray tracing

Two geometric-optics subsystems complement the PIC solvers:

- **Hybrid-branch WKB rays** (`trace_ray`, `AnalyticRayMedium`,
  `GridRayMedium`): Hamiltonian ray tracing of the warm Hall-MHD branches
  (slow / ion-cyclotron / fast-whistler) through analytic media or snapshots
  of the simulation's own `n` and `B` fields, in the code's normalized units.
  Dispersion roots are verified against the independent HYB-006 eigenvalue
  oracle (no shared code path); dispersion-function derivatives are
  complex-step exact, with spectral medium gradients for grid media; a
  recorded dispersion residual monitors accuracy along each ray.
- **`HybridPlasmaPIC.Raycon`**: a full Julia port of the RAYCON MATLAB
  package (Tracy, Kaufman & Jaun) for RF ray tracing **with linear mode
  conversion** in tokamak plasmas — Solovev equilibrium, cold-plasma
  eigenvalue dispersion, adaptive ray integration with conversion-monitor
  events, saddle-point analysis with transmission `τ = e^{−πη²}` and
  conversion coefficient `β`, ray splitting, and WKB amplitude transport
  (focusing, Maslov caustic switching, per-species damping and power
  deposition). The user-facing interface takes a `PlasmaUnits` first argument
  and works entirely in `d_i`/`Ω_ci` units; the SI layer underneath is
  regression-pinned against the original MATLAB code, and the two layers
  upstream shipped broken or disabled (cld3x3 conversion coefficients,
  amplitude transport) are completed with corrected math. See
  `docs/superpowers/specs/2026-07-05-raycon-port-notes.md` and
  `examples/raycon_cmod.jl` (Alcator C-Mod ICRF).

## Boundaries, parallelism, and I/O

**Particle boundaries.** `apply_periodic!` (half-open `[lo, hi)` contract,
float-rounding-safe), `apply_reflecting!` (specular), `apply_absorbing!`
(removal with count); all three reject non-finite positions loudly. Field
damping toward upstream in the shock steppers is done by SAT boundary terms.

**Threads.** `deposit_scalar_threaded!` / `density_threaded!` parallelize
deposition with per-thread private grids (no atomics).

**MPI.** Real MPI (via MPI.jl) covers: Cartesian communicators matched to a
`LogicalRankLayout`, slab halo exchange for fields and ghost moments,
destination-routed particle migration (`mpi_migrate_particles!`), reduced
diagnostics, distributed per-rank checkpoint/restart, and a replicated-field
reference stepper (`mpi_init!`/`mpi_step!`) — the building blocks are tested
at 2/4/8 ranks in CI, while production-scale domain-decomposed stepping and
cluster scaling data remain future work (see the honest list below).
`benchmark/mpi_scaling.jl` is the measurement harness.
`mpi_init!`, `mpi_step!`, and `mpi_compute_moments!` also accept ordered ion
collections. Charge moments are summed across populations and ranks before the
density floor is applied; migration preserves each population separately. All
ranks must supply the same ordered charge/mass parameters and stepping settings.
Use `migrate_particles=false` with particle views; migration requires Vector
storage. `save_mpi_checkpoint` and `load_mpi_checkpoint!` accept scalar populations
and ordered collections in a shared directory. Recoverable save failures before
publication preserve the previous checkpoint. Keep the same grid/layout and
exclusive directory access during checkpoint I/O; schema 1 checkpoints remain
readable. A cleanup error after publication leaves the new checkpoint usable.

**I/O.** `save_checkpoint` / `load_checkpoint!` restart `HybridStepper`,
`CAMCLStepper`, `Electrostatic1D`, `ElectrostaticPIC`, `EMPIC1D`, `EMPIC`, and `HallMHDState` runs bitwise-identically
for the same stepping inputs and numerical environment. Checkpoints preserve particle provenance and
timestep history, including CAM-CL's persistent magnetic copies and subcycle phase.
Evolving-electron runs also preserve stored pressure through serial and MPI
checkpoints and archives; missing or invalid saved pressure is rejected.
Hybrid restores replace model/shape parameters with the saved values;
electrostatic destinations require matching background density `n0`; ES
collections also require the same neutrality policy and charge tolerance.
Model/shape parameters must be independent of evolving solver and particle
storage; capture and restore reject retained links to that storage. Custom
configuration uses Julia Serialization. Preflight rejects native FFT plans,
non-null raw pointers and sharing that serialization copies independently.
Independent views and shared Array/Ref parameters remain supported. Custom
serializers own their reconstruction semantics; already-detached saved graphs
cannot recover their original links. See the
[configuration contract and limits](docs/CHECKPOINT_CONFIGURATION_PLAN.md).
Hybrid and electrostatic collections use `save_checkpoint(path, st, species)` and
`load_checkpoint!(st, species, path)`, preserving population order, charge, mass
and particle data. Destination storage is retained; resizable arrays may change
length. Every destination is checked before any population changes;
the population count must match, and scalar/collection layouts stay distinct.
`archive_run(path, st, species; rng_seed=...)` captures the same collection state.
EM checkpoints preserve both species, field/current diagnostics, configuration,
and electron/ion timestep history. For mobile ions, use
`save_checkpoint(path, es, electrons, ions)` and
`load_checkpoint!(es, electrons, ions, path)`; the destination must have matching
ion mobility and independent species storage. `archive_run` accepts the same species.
Hall-MHD uses `save_checkpoint(path, hall)`, `load_checkpoint!(hall, path)`, and
`archive_run(path, hall; rng_seed=...)`. It preserves the model, fluid and derived
fields, retained hyperresistive stage values, and clock without particle arguments.
`save_run`/`load_run`/`archive_run`
add schema versioning, SHA-256 checksums, and full reproducibility metadata
(git commit, Julia version, project/manifest hashes, RNG seed, hardware).
Serial checkpoints, archives and `async_save` publish a completed file by atomic
rename; failures before publication retain the previous file. Existing rwx modes
are retained and new files use 0600. Symlinks retain their link and update the
referent; other hardlinks keep the old contents. Directory write access and no
concurrent directory mutation are required. This provides no crash-durability
guarantee. [Publication contract](docs/SERIAL_PUBLICATION_PLAN.md).
Field dumps use a self-describing binary format with zero heavy dependencies;
HDF5 output activates as a package extension when HDF5.jl is loaded, and
distributed FFTs (PencilFFTs.jl) as another for fully periodic 3D domains.
`distributed_fft_plan` is collective: all communicator ranks must provide equal
grid sizes, transforms and precision. Select `real_type=Float32` for single
precision or retain the default `Float64`. It checks configuration agreement, rejects empty
pencil partitions, and reports recoverable local validation or planning failures
on every rank. Transform storage and transposes remain owned by PencilFFTs.
[Collective contract](docs/DISTRIBUTED_FFT_COLLECTIVE_PLAN.md) and
[precision verification](docs/DISTRIBUTED_FFT_PRECISION_PLAN.md).
`distributed_fft_roundtrip_error` supports both in-place and out-of-place plans
and returns the maximum absolute roundtrip error across the communicator. It
allocates a snapshot of the input; in-place plans reuse that input as output.
[Diagnostic contract](docs/DISTRIBUTED_FFT_DIAGNOSTIC_PLAN.md).
`distributed_fourier_grid(n, L; comm)` collectively builds a periodic 3D Fourier
workspace whose Float32/Float64 precision follows the domain lengths, and
`distributed_field(grid)` allocates a rank-owned real field on its pencil.
`deriv!(out, field, grid, axis)` and `deriv(field, grid, axis)` differentiate
without gathering the global field; input and output may alias, and local
validation failures are reported on every rank before any output is written.
[Derivative contract](docs/DISTRIBUTED_DERIVATIVE_PLAN.md).
The same grid supports `gradient!`, `divergence!`, `curl!`, `laplacian!`,
`project_divfree!`, `exp_filter!` and `dealias_two_thirds!`, with vector fields
as three-tuples of rank-owned fields and the serial alias rules (curl outputs
distinct from inputs; vector outputs and projection components pairwise
distinct). All outputs are published together after the last collective check.
[Operator plan](docs/DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md).
GPU backend storage is exposed through `particle_storage_backend`. The Metal
extension includes particle gather, deposit and push kernels; the CUDA extension
currently supplies storage and telemetry. Hybrid and EM field storage remains
CPU-based, so complete GPU-resident coupled stepping is still planned.
`backend_memory_status` returns independent telemetry counters as nonnegative
integers fitting `Int`, or `nothing` when unavailable. Every `BackendMemoryStatus`
constructor enforces that counter contract; [validation details](docs/BACKEND_MEMORY_CONSTRUCTOR_PLAN.md).

## Verification

The suite is oracle-driven: analytic dispersion relations, an independent
two-fluid eigenvalue solver, conservation laws, convergence orders, and
external reference comparisons. CI runs `Pkg.test()` on Julia 1.12 and the latest
stable Julia release. The [audit tracker](docs/DEEP_DEBUG_TRACKER.md) records tested
revisions, version-specific results, and outstanding coverage. The
[implementation plan](docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) tracks
missing end-to-end capabilities, owners, dependencies and acceptance gates.
The expanded simulator goal remains incomplete; current solver support and
standalone building blocks are distinguished in that plan. Fixes and their
regression oracles are described in `CHANGELOG.md`.

| area | oracles (selection) |
|---|---|
| spectral operators | derivative exactness, `∇·(∇×A) = 0`, projection idempotence, integration by parts, Nyquist conventions, zero steady-state allocations |
| particles & coupling | Boris gyration/E×B/2nd-order convergence, shape-function partition of unity, deposit↔gather adjointness, quiet-start current-noise spectrum |
| hybrid dispersion | ion-acoustic ω = k·c_s (0.05%), Alfvén/whistler/ion-cyclotron (~1–2%), oblique magnetosonic speeds, all against the independent two-fluid eigenvalue oracle |
| integrators | energy convergence, subcycling invariance (NB = 1..8), 2D corner-mode and CAM-CL stability at `recommended_dt` |
| Hall-MHD | Ohm/RHS oracles, mass conservation, warm-electron ion-acoustic dispersion ω = k√(Ti+Te) |
| ES/EM PIC | Langmuir & two-stream, spectral Poisson in 2D/3D, charge-conservation residual at roundoff, EMPIC ↔ EMPIC1D equivalence, relativistic beam energetics |
| collisions | conservation (exact / in-expectation), weight-independent relaxation rates, monotonicity in collisionality, odd-N handling |
| shocks | Rankine–Hugoniot residuals < 1e-10 incl. the switch-on branch, Leroy 1982 reflected-ion fraction and overshoot, published Preisser 2020 reference summary (DOI 10.5281/zenodo.3697360), 1D ≡ 3D cross-checks |
| budgets & boundaries | energy convergence and fixed-ion inertia balance, boundary seam/edge contracts, bitwise checkpoint/restart |
| ray tracing | WKB dispersion roots vs the two-fluid oracle (~1e-11), stratified-medium invariants; Raycon pinned layer-by-layer against MATLAB reference data, amplitude transport vs symplectic tangent-map oracles |

The `validation/` suite includes analytic benchmarks, published reference
comparisons and instability/control checks. Its [case catalog](validation/README.md)
describes the actual comparison in each case:

```bash
julia --project=validation validation/run_validation.jl --quick   # or --all
```

The optional [NHDS dispersion](validation/27_nhds_dispersion_comparison/README.md)
and [Hybrid-VPIC shock](validation/28_hybridvpic_perp_shock/README.md) cases build
pinned external source revisions with explicitly installed toolchains. They
verify source, input, tool and output identities before using cached results,
retain failed build attempts, and copy a provenance manifest into their artifacts.
The [reference implementation tracker](docs/EXTERNAL_REFERENCE_PROVENANCE_PLAN.md)
records acceptance evidence and remaining production-dispersion/profile work.

[Production wave-history validation](validation/34_production_parallel_dispersion/README.md)
measures signed frequency and growth/damping from evolved cold parallel hybrid
and Hall-MHD modes. It retains raw histories, residuals, half-window sensitivity,
independent grid/time/particle refinements and source/environment hashes.
Conditional fit errors do not establish statistical confidence for correlated
PIC noise. Broader oblique, thermal, kinetic and EM comparisons remain in the
[production dispersion plan](docs/PRODUCTION_DISPERSION_PLAN.md).

### Honest limits (not claimed)

- **Complete GPU-resident coupled stepping and multi-GPU scaling.** Metal
  uniform-field particle pushing has passed actual-device checks on supported
  Julia versions; that does not establish full GPU solver support. See the
  [implementation and evidence](docs/METAL_UNIFORM_PUSH_PLAN.md).
- **MPI cluster scaling data**: all building blocks pass real-MPI tests at
  2–8 ranks locally/CI; production runs need a cluster.
- **Full-profile replay** of the Preisser 2020 dataset (the compact scalar
  summary with checksums is bundled; the full HDF5 requires downloading the
  upstream data).
- Non-periodic EM-PIC particle coupling (`EMPIC` and `EMPIC1D` use periodic geometry).
  The separate [bounded Maxwell vacuum solver](docs/BOUNDED_MAXWELL_IMPLEMENTATION.md)
  supplies conducting cavities and mixed periodic/conducting fields; charge-conserving
  particle coupling, absorbing layers and restart integration remain open.
- Cell-local Coulomb collisions are a standalone CPU operator. Automatic
  collision/field splitting with saved RNG state, MPI cell ownership and device
  kernels are separate open stages ([CAP-010](docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)).
- The [normalization repair](docs/NORMALIZATION_RANGE_PLAN.md)
  does not establish range correctness of every compound RayCon expression.
  The shared BigFloat exact-conversion repair has passed its
  [local acceptance and exact-revision CI checks](docs/EXACT_BIGFLOAT_CONVERSION_PLAN.md).

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Focused real-MPI checks and the scaling harness:

```bash
julia --project=. test/test_mpi.jl
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. test/test_mpi_multirank.jl`)'
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. benchmark/mpi_scaling.jl --reps 3`)'
```

## Examples

Each runs standalone in a few seconds to ~30 s (including compilation):

```bash
julia --project=. examples/ion_acoustic.jl   # hybrid run; measures ω against k·c_s (0.05%)
julia --project=. examples/two_ion_hybrid.jl # proton/alpha run against a uniform-field analytic orbit
julia --project=. examples/ray_tracing.jl    # WKB fast-wave ray through a density wave
julia --project=. examples/raycon_cmod.jl    # RAYCON: C-Mod ICRF ray + mode conversion + power deposition
```

## Package layout

```
src/
  boundaries/             particle boundaries, flux-correct injection reservoir
  coupling/               NGP/CIC/TSC shapes, deposition, gather, moments
  diagnostics/            budgets, spectra, phase space, shock & spacecraft diagnostics, WKB rays
  electrons/              closures, generalized Ohm's law, pressure evolution
  fields/                 HybridFields state container
  initial_conditions/     Rankine–Hugoniot, shocks, ramps, KdV
  integrators/            HybridStepper, CAM-CL, recommended_dt, CN whistler demo
  io/                     checkpoint/restart, metadata, field dumps
  meshes/                 Cartesian, local FD, toroidal curvilinear, and AMR mesh utilities
  models/                 hybrid PIC, Hall-MHD, electrostatic/electromagnetic PIC, gyrokinetics
  parallel/               rank layouts, MPI, threads, GPU/extension hooks
  particles/              ParticleSet{D}, loaders, Boris push, collisions, sorting
  raycon/                 RAYCON port (tokamak RF mode conversion)
  utils/                  argument-validation helpers
  verification/           analytic oracles, shock/instability campaigns, PlasmaUnits
```

Reusable Fourier, filtering, projection, and mixed SBP/Fourier operators live
in [`SpectralOperators.jl`](https://github.com/jake-w-liu/SpectralOperators.jl)
and are re-exported here.

## References

- A.P. Matthews, *Current advance method and cyclic leapfrog for 2D multispecies
  hybrid plasma simulations*, J. Comput. Phys. **112**, 102 (1994).
- T. Takizuka & H. Abe, *A binary collision model for plasma simulation with a
  particle code*, J. Comput. Phys. **25**, 205 (1977).
- D.P. Higginson, I. Holod & A. Link, *A corrected method for Coulomb scattering
  in arbitrarily weighted particle-in-cell plasma simulations*, J. Comput.
  Phys. **413**, 109450 (2020).
- M.M. Leroy et al., *The structure of perpendicular bow shocks*, J. Geophys.
  Res. **87**, 5081 (1982).
- E.R. Tracy, A.N. Kaufman & A. Jaun, Phys. Lett. A **290**, 309 (2001);
  A. Jaun, E.R. Tracy & A.N. Kaufman, Plasma Phys. Control. Fusion **49**, 43
  (2007); E.R. Tracy, A.N. Kaufman & A. Jaun, Phys. Plasmas **14**, 082102
  (2007) — the RAYCON mode-conversion papers.
- L. Preisser et al. 2020 dataset, Zenodo, DOI
  [10.5281/zenodo.3697360](https://doi.org/10.5281/zenodo.3697360) (CC-BY-4.0).
