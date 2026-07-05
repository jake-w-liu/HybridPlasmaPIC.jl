# HybridPlasmaPIC.jl

A dimension-parametric (1D3V / 2D3V / 3D3V) **hybrid particle-in-cell** plasma
solver: kinetic ions (Boris-pushed macroparticles), massless fluid electrons
(generalized Ohm's law), and spectral fields. The package is the standalone
hybrid-PIC solver split out of the plasma-solver implementation checklist.

Design rules followed: **design for 3D, verify in 1D**; no global mutable state
(explicit `HybridStepper`/`ParticleSet`/`HybridFields`); the spatial dimension
`D` is a type parameter; particle positions carry exactly `D` coordinates while
velocities always carry 3.

## Governing equations (massless-electron hybrid, Ω_ci-normalized)

```
ions:   dx/dt = v,   dv/dt = E + v×B            (Boris, q/m = 1 for protons)
        J = ∇×B
        u_e = u_i − J/n
        E = −u_i×B + (J×B)/n − ∇p_e/n + η J      (η optional, off by default)
        ∂B/∂t = −∇×E,    ∇·B = 0
```

Electron closures: isothermal `p_e = T_e·n` or polytropic `p_e = p_e0 (n/n0)^γ`.

### Normalization (§7 of the checklist)

| quantity | unit |
|---|---|
| length | ion inertial length `d_i = c/ω_pi` |
| time   | `Ω_ci^{-1}` |
| velocity | Alfvén speed `v_A` |
| B | `B0` |
| E | `v_A B0` |
| density | `n0` |
| J | `B0/(μ0 d_i)` |
| pressure | `B0²/μ0` |

Use `set_density_weight!(ps, 1.0, g)` so a uniform load deposits `n = n0 = 1` —
otherwise every `1/n` term in Ohm's law is silently rescaled.

## Time integration

Explicit hybrid leapfrog with documented time levels
(`src/integrators/extrapolated_leapfrog.jl`):
the Boris push is centered at step `n`; `B` advances `n→n+1` via an RK4-subcycled
Faraday step with the `n+1/2` ion moments frozen (CAM/CL structure); the carried
`E` is recomputed each step. Magnetic subcycling (`NB`) handles the whistler CFL
(`Ω_ci Δt ≲ (Δx/d_i)²`).

## Ray tracing (supporting methods)

Two geometric-optics subsystems complement the PIC solvers:

- **Hybrid-branch WKB rays** (`trace_ray`, `AnalyticRayMedium`,
  `GridRayMedium`): Hamiltonian ray tracing of the warm Hall-MHD branches
  (slow / Alfvén-ion-cyclotron / fast-whistler) through analytic media or
  snapshots of the simulation's own `n` and `B` fields, in the code's
  Ω_ci-normalized units. The scalar dispersion relation is verified against
  the HYB-006 two-fluid eigenvalue oracle; Hamiltonian derivatives are exact
  (complex-step); the recorded `|D|` residual monitors integration accuracy.
- **`Raycon`** (`HybridPlasmaPIC.Raycon`): a full Julia port of the RAYCON
  MATLAB package (Tracy, Kaufman & Jaun — PLA 290 (2001) 309; PPCF 49 (2007)
  43; PoP 14 (2007) 082102) for RF ray tracing **with linear mode conversion**
  in tokamak plasmas: Solovev equilibrium, cold-plasma eigenvalue dispersion,
  adaptive ray integration with conversion-monitor events, saddle-point
  analysis with transmission τ = e^{−πη²} and conversion coefficient β, and
  ray splitting. **The user-facing interface follows this package's Ω_ci
  normalization**: pass a `PlasmaUnits` as the first argument
  (`RayconProblem(units; …)`, `launch_ray(units, prob; …)`,
  `trace_rays(units, prob; …)`) and work in d_i / Ω_ci / v_A / B0 / n0
  throughout; the same functions without `PlasmaUnits` are the raw SI engine
  that is regression-pinned against the original MATLAB code
  (`cmod_units()` + `cmod_parameters(units)` give the Alcator C-Mod ICRF
  reference case). See
  `docs/superpowers/specs/2026-07-05-raycon-port-notes.md` for the file-by-file
  mapping, preserved upstream quirks, and the upstream bugs found & fixed.

## Verification status

Benchmarks below are checked against **independent analytic oracles**.
Post-extraction verification on Julia 1.12.6: `Pkg.test()` passed 100,630 tests.
Tolerances are the checklist's initial engineering targets.

| Benchmark | What | Status |
|---|---|---|
| OP-001/002/003/005/006 | derivative, div(curl)=0, projection, IBP, Nyquist | ✅ verified |
| OP-004 / OP-005-SBP | mixed SBP-x/Fourier-y operator, SBP identity, convergence | ✅ verified |
| zero-alloc operators | steady-state allocations = 0 | ✅ verified |
| PUSH-001..004 | gyration, E×B, ∥E, 2nd-order convergence | ✅ verified |
| LOAD-001/002 | Maxwellian moments, quiet start | ✅ verified |
| DEP-001..005 | partition of unity, conservation, linear gather, translation, N⁻¹ᐟ² noise | ✅ verified |
| HYB-001 | uniform equilibrium stationary | ✅ verified |
| HYB-002 | ion-acoustic ω = k·c_s (0.05%) | ✅ verified |
| HYB-003/004 | Alfvén / whistler / ion-cyclotron branches (~1–2%) | ✅ verified |
| HYB-005 | fast/slow magnetosonic c_{f,s} (oblique) + PIC fast mode | ✅ verified |
| HYB-006 | independent two-fluid eigenvalue oracle (`test/oracles/`) | ✅ verified |
| HYB-007/008 | adiabatic energy convergence, subcycling (NB=1,2,4,8) | ✅ verified |
| HALL-MHD | periodic Hall-MHD Ohm/RHS/RK4 mass conservation oracles | ✅ verified |
| ESPIC | electrostatic PIC: 1D Langmuir/two-stream + 2D/3D Poisson oracles | ✅ verified |
| EMPIC | EM PIC: 1D Esirkepov + 2D/3D spectral charge-conserving current | ✅ verified |
| KDV-001 | KdV soliton + 2/3 dealiasing | ✅ verified |
| RAY-001..009 | WKB hybrid-branch ray tracing: dispersion roots vs HYB-006 oracle (~1e-11), group velocity, stratified k(x)=ω/c(x) invariant, 1D≡2D≡3D media | ✅ verified |
| RCN-001..011 | `Raycon` (RAYCON port, tokamak RF + mode conversion): FD oracles for every derivative layer, exact identities, C-Mod conversion run (τ, β) | ✅ verified |
| SHK-001 | Rankine–Hugoniot solver, residuals < 1e-10 | ✅ verified |
| SHK-005 | Published external hybrid-code reference metadata + scalar comparison target | ✅ verified |
| DIM-001 | 1D ≡ y-invariant 2D (operators + integrator) | ✅ verified |
| DFFT | optional PencilFFTs/PencilArrays distributed FFT extension, COMM_SELF parity | ✅ verified |
| IO-001 | checkpoint/restart bitwise | ✅ verified |

### Remaining Gated Work

These require hardware or external references not available here and are **not**
claimed as done:

- **GPU** (CUDA/Metal) production plasma kernels — no NVIDIA hardware here.
- **MPI** cluster scaling data. Focused real-MPI Cartesian mapping, diagnostic
  Allreduce, destination-routed particle migration, slab field/moment halo
  exchange, time-advanced particle budget invariance, distributed
  checkpoint/restart bitmatch, field-coupled serial-vs-MPI agreement, optional
  PencilFFTs distributed FFT one-rank parity, and the `benchmark/mpi_scaling.jl`
  harness pass locally. Production scaling still requires cluster runs.
- **External hybrid-code full-profile replay**: SHK-005 now bundles DOI,
  checksum, and scalar summaries from the published Preisser et al. 2020 Zenodo
  hybrid-simulation dataset (`10.5281/zenodo.3697360`). Replaying the full HDF5
  profiles remains optional and requires downloading the upstream data.
- **Advanced multidimensional full-PIC features**: `EMPIC` covers periodic
  1D/2D/3D electromagnetic PIC with spectral charge-conserving current
  correction. Multidimensional electron subcycling and non-periodic EM PIC
  boundary conditions are not claimed.
- **Multi-GPU** restart and scaling.

## Package layout

```
src/
  boundaries/             particle and field boundary helpers
  coupling/               NGP/CIC/TSC shapes, deposition, gather, moments
  diagnostics/            conservation, spectra, phase-space, shock diagnostics
  electrons/              closures, generalized Ohm's law, pressure evolution
  fields/                 HybridFields state container
  initial_conditions/     waves, shocks, ramps, and uniform plasma setup
  integrators/            HybridStepper, CAM-CL, semi-implicit helpers
  io/                     checkpoint/restart and reproducibility metadata
  meshes/                 Cartesian and local finite-difference mesh adapters
  models/                 HybridPIC, electrostatic/full-PIC, and Hall-MHD models
  parallel/               threaded CPU/MPI helpers plus GPU and I/O extension hooks
  particles/              ParticleSet{D}, loaders, Boris mover, sorting
  verification/           analytic oracles and campaign helpers
```

Reusable Fourier, filtering, projection, and mixed SBP/Fourier operators live in
`SpectralOperators.jl` and are re-exported by HybridPlasmaPIC. `Project.toml`
resolves that package from `https://github.com/jake-w-liu/SpectralOperators.jl.git`,
so a sibling checkout is not required for development from this repository.

## Installation

Until `SpectralOperators.jl` is registered, add it explicitly before adding
HybridPlasmaPIC:

```julia
using Pkg
Pkg.add(url = "https://github.com/jake-w-liu/SpectralOperators.jl.git")
Pkg.add(url = "https://github.com/jake-w-liu/HybridPlasmaPIC.jl.git")
```

Requires Julia 1.12 or later.

## Run the tests

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Focused MPI checks:

```bash
julia --project=. test/test_mpi.jl
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project=. test/test_mpi_multirank.jl`)'
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. test/test_mpi_multirank.jl`)'
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 8 $(Base.julia_cmd()) --project=. test/test_mpi_multirank.jl`)'
```

MPI scaling smoke/benchmark harness:

```bash
julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. benchmark/mpi_scaling.jl --reps 3`)'
```

## Minimal example

See `examples/ion_acoustic.jl`. Sketch:

```julia
using HybridPlasmaPIC
g  = FourierGrid((64,), (2π,))
ps = ParticleSet{1,Float64}(64*400)
load_lattice_1d!(ps, 0.0, 2π); set_density_weight!(ps, 1.0, g)
load_quiet_velocities!(ps, MersenneTwister(1), (0,0,0), (0,0,0))   # cold
for p in eachindex(ps.weight); ps.v[1][p] += 0.01*sin(ps.x[1][p]); end
st = HybridStepper(g, HybridModel(IsothermalElectrons(1.0)), CIC(), nparticles(ps))
init!(st, ps)
for _ in 1:1000; step!(st, ps, 0.02); end
```
