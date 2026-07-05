# HybridPlasmaPIC.jl

A research-grade **hybrid particle-in-cell (PIC) plasma solver** in Julia:
kinetic ions (Boris-pushed macroparticles), massless fluid electrons
(generalized Ohm's law), and spectral fields on periodic grids — plus a family
of companion solvers (Hall-MHD, electrostatic PIC, full electromagnetic PIC),
collision operators, shock drivers, geometric-optics ray tracing, and an
extensive analytic-oracle verification suite (108,600 passing tests).

Everything is dimension-parametric: the same code runs 1D3V, 2D3V, and 3D3V
(positions carry exactly `D` coordinates, velocities always carry 3), with no
global mutable state — simulations are explicit
`HybridStepper` / `ParticleSet` / `HybridFields` objects you own.

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

Requires **Julia 1.12 or later**. Until `SpectralOperators.jl` (the sibling
package providing the Fourier/SBP operators) is registered, add it first:

```julia
using Pkg
Pkg.add(url = "https://github.com/jake-w-liu/SpectralOperators.jl.git")
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
`inertial_length` give the reference scales.

### Electron closures

- `IsothermalElectrons(Te)` — `p_e = Te·n`
- `PolytropicElectrons(pe0, n0, γ)` — `p_e = pe0·(n/n0)^γ`
- `CGLElectrons(p⊥0, p∥0, n0, B0)` — anisotropic double-adiabatic
  (Chew–Goldberger–Low): `p_⊥ ∝ n·B`, `p_∥ ∝ n³/B²`, entering Ohm's law
  through the full gyrotropic pressure-tensor force `−(∇·P_e)/n`. (Note the
  fluid mirror mode has no short-wavelength cutoff, so strongly
  mirror-unstable CGL states blow up — a property of the CGL model itself.)

A standalone adiabatic electron pressure-evolution step
(`advance_electron_pressure!`) is exported as a building block; the built-in
steppers use the algebraic closures above.

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
  the carried `E` recomputed each step. Supports every closure, including CGL.
- **`CAMCLStepper`** (`init_camcl!`, `step_camcl!(st, ps, dt; NB ≥ 2)`) —
  Matthews' CAM-CL (J. Comput. Phys. 112, 102, 1994): current advance method
  plus a cyclic-leapfrog `B` subcycle whose two staggered field copies persist
  across particle steps (neutrally stable, non-dissipative, for whistler
  `ω·h < 2` per substep; two Ohm evaluations per substep instead of RK4's
  four). Scalar closures only.

**Choosing dt.** `recommended_dt(g; NB = 1, integrator = :rk4, safety = 0.8,
η = 0.0, ηH = 0.0)` returns a conservative particle step from the whistler CFL
evaluated at the spectral-corner wavenumber `kmax = √Σ_d (π/Δx_d)²` (the
stiffest representable mode in 2D/3D — a single-axis Nyquist estimate is up to
`D×` too optimistic), additionally capped by Boris gyro-accuracy
(`Ω_ci·Δt ≤ 0.3`) and, when `η`/`ηH` are supplied, by the resistive real-axis
stability of the chosen integrator. In 1D, in the fine-grid limit
`kmax·d_i ≫ 1`, it reduces to the familiar `Ω_ci Δt ≲ (1/π)(Δx/d_i)²`.

## Beyond the hybrid model

| model | what it is | entry points |
|---|---|---|
| **Hall-MHD** | periodic single-fluid Hall-MHD with the same electron closures; ion pressure `p_i = Ti·n`; the electron pressure force acts in both Ohm's law and the bulk momentum equation; explicit RK4 with divergence-free projection | `HallMHDModel`, `HallMHDState`, `step_hall_mhd!` |
| **Electrostatic PIC** | full-kinetic electrons over a neutralizing background; spectral Poisson solve; leapfrog with one-time velocity priming; 1D/2D/3D | `ElectrostaticPIC` (or `Electrostatic1D`), `init_espic!`, `step_espic!` |
| **Electromagnetic PIC** | full EM PIC with kinetic electrons, optional mobile ions (`mobile = true`), optional relativistic Boris push (`relativistic = true`), spectral charge-conserving current (continuity holds to roundoff on every represented mode), Gauss-law initialization; `EMPIC1D` additionally offers electron subcycling `n_sub` | `EMPIC`, `EMPIC1D`, `init_empic!`, `step_empic!`, `charge_conservation_residual` |
| **KdV reference** | integrating-factor pseudo-spectral KdV solver + analytic soliton, used as a nonlinear-wave verification target | `kdv_solve`, `kdv_soliton` |

The PIC models share the grid, particle, and diagnostic infrastructure with
the hybrid stepper, so moving a setup between them is mostly a constructor
swap; Hall-MHD shares the grid and field machinery, and the KdV solver is a
standalone verification tool.

## Collisions and extra particle physics

- `collide_coulomb!(ps, gcoeff, dt)` — Takizuka–Abe (1977) binary Coulomb
  scattering, with the Higginson–Holod–Link (JCP 2020) correction for unequal
  macro-particle weights (relaxation rates independent of the weight
  distribution; conservation exact for equal weights, in expectation
  otherwise), an isotropic large-angle fallback when the sampled scattering
  variance exceeds unity, and the T&A triplet for odd particle counts.
- `collide_bgk!(ps, ν, dt)` — BGK relaxation toward the set's own drifting
  Maxwellian, conserving momentum and energy exactly.
- `collide_neutral_mcc!` / `ionize_mcc!` — Monte-Carlo elastic scattering off
  a neutral Maxwellian reservoir, and electron-impact ionization with exact
  energy-threshold bookkeeping and unique newborn particle ids.
- `apply_radiation_reaction!(ps, E, B, dt; K, c)` — the dominant `γ²` term of
  the reduced Landau–Lifshitz drag, integrated as an exact exponential decay
  (unconditionally stable).

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
  checkpoint/restart.
- **Jump conditions**: `rankine_hugoniot` / `rh_branches` solve the coplanar
  MHD Rankine–Hugoniot system with residuals of all six conservation laws
  reported; the field-aligned (`Bt = 0`) bifurcation is handled explicitly —
  inside the switch-on window the solver constructs the switch-on branch in
  closed form and selects the evolutionary solution.
- **Boundary machinery**: flux-correct open-boundary injection
  (`ShockInjector`, `flux_speed`, `inject_face_1d!`) samples inflow velocities
  from the drifting-Maxwellian flux distribution `p(v) ∝ v·f_M`.

## Diagnostics

- **Energy**: `energy_budget(ps, B, n, closure, g; de2 = 0)` returns kinetic,
  magnetic, electron-internal, and (for `de2 > 0`) electron-inertia
  reservoirs, and closes (total conserved at `η = 0`) for every wired
  closure — exactly for the scalar closures (polytropic `∫p_e/(γ−1)`,
  isothermal via the free energy `T_e∫n ln n dV`, its exact `γ→1` limit) and
  modulo the small anisotropic battery term for CGL (gyrotropic
  `∫(p_⊥ + p_∥/2) dV`). `kinetic_energy_relativistic(ps, c)` covers
  relativistic EMPIC runs. `resistive_dissipation`, `jdotE_density`,
  `electric_work`, and `momentum_budget` complete the ledgers (with honest
  docstrings about which channel receives what).
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

**I/O.** `save_checkpoint` / `load_checkpoint!` restart runs
bitwise-identically (stepping uses no RNG); `save_run`/`load_run`/`archive_run`
add schema versioning, SHA-256 checksums, and full reproducibility metadata
(git commit, Julia version, project/manifest hashes, RNG seed, hardware).
Field dumps use a self-describing binary format with zero heavy dependencies;
HDF5 output activates as a package extension when HDF5.jl is loaded, and
distributed FFTs (PencilFFTs.jl) as another for fully periodic 3D domains.
GPU array-backend plumbing (CUDA/Metal extensions) exists behind
`particle_storage_backend`; production GPU kernels are not claimed.

## Verification

The suite is oracle-driven: analytic dispersion relations, an independent
two-fluid eigenvalue solver, conservation laws, convergence orders, and pinned
external references — not just smoke tests. `Pkg.test()` runs 62 test files,
**108,600 tests**, all passing on Julia 1.12. A full physics/math audit
(2026-07) re-derived every subsystem from first principles, adversarially
verified 37 findings, and fixed all of them with new discriminating regression
oracles (see `CHANGELOG.md` for the complete list).

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
| budgets & boundaries | closed energy budgets for every closure, boundary seam/edge contracts, bitwise checkpoint/restart |
| ray tracing | WKB dispersion roots vs the two-fluid oracle (~1e-11), stratified-medium invariants; Raycon pinned layer-by-layer against MATLAB reference data, amplitude transport vs symplectic tangent-map oracles |

There is also a 32-study `validation/` suite (analytic benchmarks through
published-physics reproductions — Leroy 1982, Hellinger 2002, NHDS dispersion,
Hybrid-VPIC shock, firehose/ion-cyclotron/Weibel/reconnection growth rates):

```bash
julia --project=validation validation/run_validation.jl --quick   # or --all
```

### Honest limits (not claimed)

- **GPU** production kernels and multi-GPU scaling (no hardware here; the
  array-backend plumbing and extensions exist).
- **MPI cluster scaling data**: all building blocks pass real-MPI tests at
  2–8 ranks locally/CI; production runs need a cluster.
- **Full-profile replay** of the Preisser 2020 dataset (the compact scalar
  summary with checksums is bundled; the full HDF5 requires downloading the
  upstream data).
- Multidimensional electron subcycling and non-periodic EM-PIC boundaries
  (`EMPIC` is periodic; `n_sub` belongs to `EMPIC1D`).

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
  meshes/                 pointer notes (concrete operators live in SpectralOperators.jl)
  models/                 hybrid PIC, Hall-MHD, electrostatic & electromagnetic PIC
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
