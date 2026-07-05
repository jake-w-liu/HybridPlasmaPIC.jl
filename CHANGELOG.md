# Changelog

All notable changes to HybridPlasmaPIC.jl are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **WKB ray tracing for the hybrid wave branches**
  (`src/diagnostics/ray_tracing.jl`): warm Hall-MHD scalar dispersion relation
  (verified against the HYB-006 eigenvalue oracle to ~1e-11), branch
  frequencies/wavenumbers/group velocities, complex-step Hamiltonian
  derivatives, and RK4 `trace_ray` through analytic media or snapshots of the
  simulation's own grid fields (`AnalyticRayMedium`, `GridRayMedium`, 1D/2D/3D).
  Tests RAY-001..009.
- **`Raycon` submodule** (`src/raycon/`): full Julia port of the RAYCON
  MATLAB package (Jaun–Kaufman–Tracy) for RF ray tracing with linear mode
  conversion in tokamak plasmas — Solovev equilibrium, magnetic geometry with
  basis-vector curvature derivatives, cold-plasma eigenvalue dispersion
  (cld2x2/cld3x3), adaptive DP4(5) ray integration with conversion-monitor
  events, saddle-point mode-conversion analysis with transmission
  `τ = e^{−πη²}` and conversion coefficient `β` (complex Γ), ray splitting,
  and the Alcator C-Mod ICRF reference case. The user-facing interface is
  unified with the package's Ω_ci normalization via `PlasmaUnits`-first
  methods (`RayconProblem(units; …)`, `launch_ray(units, …)`,
  `trace_rays(units, …)`, `integrate_ray(units, …)`, `cmod_units()`); the
  SI layer remains as the MATLAB-pinned engine. Beyond the port, two layers
  upstream ships broken or disabled are completed: cld3x3 mode-conversion
  coefficients (exact determinant-sampling derivatives + Tracy–Kaufman
  near-null-subspace coupling, reducing to the pinned 2×2 when the
  electrostatic branch decouples) and WKB amplitude transport
  (`integrate_ray_amplitude`, `trace_rays(amplitude=true)`: Riccati focusing
  tensor, lnE²/eikonal phase, Maslov caustic switching, cyclotron/Landau/TTMP
  damping with per-species deposition, conversion amplitude bookkeeping),
  verified against symplectic tangent-map identities, stationary-phase
  representation checks, and exact energy bookkeeping. Port notes, preserved
  upstream quirks, and the upstream bugs found & fixed are documented in
  `docs/superpowers/specs/2026-07-05-raycon-port-notes.md`. Tests RCN-001..015
  plus regression against reference data dumped from the original MATLAB code
  (`tools/raycon_reference.m`).
- **Package quality and CI hardening.** JET static-analysis smoke test
  (report-only), a `.JuliaFormatter.toml` (default style, 100-column margin), a
  formatting-check CI job, a scheduled (cron) CI run, a Windows runner in the
  test matrix for the non-MPI core, and upload of test logs as CI artifacts.
- Project metadata: `CITATION.cff`, `CONTRIBUTING.md`, and this changelog.
- SHK-005 now includes a compact published external hybrid-code reference from
  the Preisser et al. 2020 Zenodo dataset (`10.5281/zenodo.3697360`): source DOI,
  license, HDF5 checksum, derived scalar `Bavg_y` summaries, and comparison tests.
- Implemented periodic Hall-MHD with `HallMHDModel`, `HallMHDState`,
  generalized Ohm's law, continuity/momentum/Faraday RHS, RK4 stepping, and
  analytic verification tests.
- Added dimension-parametric `ElectrostaticPIC` for periodic 1D/2D/3D
  electrostatic full-PIC, with 2D/3D spectral Poisson oracles and 2D validation
  and equilibrium tests while preserving the existing `Electrostatic1D` API.
- Added dimension-parametric `EMPIC` for periodic 1D/2D/3D electromagnetic
  full-PIC with mobile electrons, optional mobile ions, relativistic Boris
  support, spectral Gauss-law initialization, and spectral current correction
  that enforces charge conservation on represented Fourier modes.
- Added optional `HybridPlasmaPICPencilFFTSExt` integration with
  `PencilFFTs.jl`/`PencilArrays.jl` for fully periodic 3D distributed FFT plans,
  input/output allocation, forward/inverse transforms, and round-trip checks.

### Fixed

- **Full physics/math audit of the simulator core** (37 adversarially verified
  findings; every fix carries a new discriminating regression oracle):
  - Hall-MHD momentum equation now includes the electron-pressure force
    `−∇p_e/n` (scalar closures) / `−(∇·P_e)/n` (CGL). Warm-electron runs
    previously had no electron acoustic physics at all; ion-acoustic dispersion
    now verified against `ω = k√(Ti+Te)` to 8 digits, and `Te = 0` trajectories
    are bit-identical to before.
  - Electron-inertia filter `1/(1+d_e²k²)` now acts only on the k-transverse
    projection of `E`; the longitudinal/electrostatic field used for the ion
    push is no longer spuriously damped (was up to ~53% wrong at `k·d_e ~ 1`).
  - Quiet-start velocity pairing switched from `i ↔ i+N/2` to adjacent lattice
    indices: deposited thermal-current noise at odd harmonics drops ~900× below
    the random-load level (the old pairing *amplified* it by √2). Odd particle
    counts are now accepted (last particle carries the drift exactly).
  - `recommended_dt` uses the spectral-corner `|k|max = √Σ(π/dx_d)²` (the old
    single-axis Nyquist under-resolved multi-D whistlers by up to D×, blowing
    up 2D diagonal-B₀ runs at the recommended step) and accepts `η`/`ηH`
    keywords to cap the step by the resistive/hyper-resistive real-axis rates.
  - CAM-CL cyclic-leapfrog B-subcycle now keeps its two staggered copies
    persistent across particle steps (Matthews-style output-only averaging,
    drift-triggered re-sync). The previous restart-and-average composition was
    weakly unstable at every `ω·h > 0`; the scheme is now neutrally stable to
    its true measured limit `ω·h ≤ 2` and remains 2nd-order accurate.
  - EMPIC's spectral continuity correction now zeroes pure-Nyquist deposited
    current (all-`k`-zeroed mode combinations other than DC); EMPIC1D likewise
    zeroes the transverse Nyquist current. This removes a random-walking
    grid-Nyquist sawtooth in `E` and restores EMPIC ↔ EMPIC1D equivalence to
    machine precision.
  - Rankine-Hugoniot solver handles the field-aligned (`Bt1 = 0`) bifurcation:
    `rh_branches` now returns the switch-on branch (closed-form `X = M_An²`
    solution, flux residuals ~1e-16) and `rankine_hugoniot` prefers the
    evolutionary branch, instead of silently returning the non-evolutionary
    gasdynamic compression (1.74× wrong in the switch-on window).
  - Leroy shock-boundary reinsertion now draws flux-weighted velocities
    (`flux_speed`) at both ends instead of density-weighted half-Maxwellians
    (removes a ~30% spurious boundary-layer density pile-up); `PerpShock` with
    `closure = :energy` now advances `p_e` in `step_shock!` (it was frozen at
    the initial profile); `flux_speed`'s Rayleigh branch uses `log1p(−u)` so a
    `rand() == 0` draw can no longer inject an infinite-speed particle.
  - Takizuka–Abe Coulomb collisions: unequal macro-weight pairs now scatter
    symmetrically about the true velocity midpoint with the
    Higginson–Holod–Link (JCP 2020) rejection/variance-scaling correction
    (relaxation rates are now weight-independent; weights previously acted as
    masses, biasing rates by O(1)); a large-angle isotropic fallback replaces
    the Gaussian `tan(Θ/2)` sample when `⟨δ²⟩ > 1` (relaxation now saturates
    monotonically with collisionality instead of stalling); odd particle
    counts use the T&A triplet at half variance (no particle skipped).
  - Energy budgets close for every wired closure: isothermal
    `electron_internal_energy` returns the free energy `T_e∫n ln n dV` (the
    exact `γ→1` invariant — the documented "no closed invariant" claim was
    false), CGL runs get the gyrotropic `∫(p_⊥+p_∥/2) dV` via a new
    `electron_internal_energy(n, B, closure, g)` method, and
    `energy_budget(...; de2)` adds the electron-inertia reservoir
    `∫d_e²|J|²/2 dV`. New exported `kinetic_energy_relativistic(ps, c)`
    completes energy bookkeeping for relativistic EMPIC runs. Docstrings for
    `resistive_dissipation`, `jdotE_density`, and `momentum_budget` no longer
    misattribute energy/momentum channels.
  - `weibel_growth` classifies stability with the finite-box criterion
    `A > (c·k_min/ω_pe)²` (the previous `A > 1` threshold does not exist in
    the linear theory); `reproduce_established_shock` gates the measured
    compression against the run's actual `X_RH` instead of the `M→∞` ceiling.
  - Boundary/indexing edge cases: `apply_periodic!` can no longer return
    `x == hi` (half-open contract restored), `apply_absorbing!` rejects
    non-finite positions instead of counting blow-ups as absorption, and
    `cell_index` wraps periodically exactly like the deposition mesh instead
    of clamping.
  - `semi_implicit` whistler demo gives the Nyquist bin its physical
    `ω = c·k²` (the field is complex; no reality constraint), and the
    radiation-reaction docstring states the `c`-dependent cooling law
    `d(1/γ)/dt = K v⁴B²/c²`.
  - Replaced ~10 tautological test oracles (tests that encoded the same wrong
    assumption as the implementation) with discriminating physics oracles:
    quiet-start current spectrum, longitudinal-inertia invariance, 2D
    diagonal-B₀ and CAM-CL stability at `recommended_dt`, switch-on RH window,
    warm-electron ion-acoustic dispersion, unequal-weight collision rates and
    `gcoeff` monotonicity, isothermal/CGL/inertia budget conservation, EMPIC
    Nyquist exclusion, and boundary-seam contracts (`test_boundary_edges.jl`).
- Hall-MHD CGL electron closures now use the anisotropic pressure-force Ohm-law path.
- Hybrid and CAM-CL steppers now resize particle workspaces after particle injection/removal.
- Logical rank layouts now reject decompositions whose total rank count overflows `Int`.

## [0.1.0]

First tagged release. HybridPlasmaPIC.jl is a dimension-parametric (1D3V / 2D3V /
3D3V) hybrid particle-in-cell plasma solver with kinetic ions, massless fluid
electrons, and spectral fields, in Ω_ci-normalized units. There is no global
mutable state; all configuration is passed explicitly.

### Spectral operators (`SpectralOperators.jl`, `src/meshes/local_finite_difference.jl`)

- `FourierGrid(n, L)` periodic spectral grid carrying angular wavenumbers
  (`kvec`, Nyquist mode zeroed), cached forward/inverse FFT plans, and scratch
  buffers.
- `deriv!`/`deriv`, `gradient!`, `divergence!`, `curl!`, `laplacian!`, and
  `project_divfree!` for divergence cleaning.
- `SBP1D` summation-by-parts finite-difference operator and `sbp_deriv!` /
  `fourier_deriv_y!` for mixed FD-x / spectral-y discretizations.

### Particles (`src/particles/`, `src/boundaries/`)

- `ParticleSet{D,T}` structure-of-arrays macroparticle container (positions
  carry exactly `D` coordinates; velocities always carry 3).
- Loaders: `load_uniform!`, `load_lattice!`, `load_lattice_1d!`,
  `load_maxwellian!`, `load_quiet_velocities!`, `set_density_weight!`.
- `boris_kick` relativistic-free Boris velocity rotation; `push_uniform!`.
- Boundary conditions: `apply_periodic!`, `apply_reflecting!`,
  `apply_absorbing!`.

### Deposition and moments (`src/coupling/`)

- Shape functions `NGP`, `CIC`, `TSC`.
- `deposit_scalar!`, `gather_scalar!`, `gather_vector!`, and the moment kernels
  `density!`, `momentum!`, `current!`, `pressure_tensor!`,
  `temperature_components`.

### Hybrid model and integrator (`src/models/hybrid_pic.jl`, `src/integrators/`)

- Electron closures `IsothermalElectrons`, `PolytropicElectrons`.
- `HybridModel` / `HybridFields`, generalized Ohm's law (`ohms_law!`), Faraday
  RHS (`faraday_rhs!`), divergence cleaning (`project_b!`), and moment
  assembly (`compute_moments!`, `compute_moments_multi!`).
- `HybridStepper` with `init!` / `step!`, plus energy diagnostics
  (`kinetic_energy`, `magnetic_energy`, `electron_internal_energy`) and
  `mode_amplitude`.

### Dispersion benchmarks (`test/test_dispersion.jl`)

- Whistler / Alfvén-branch dispersion verified against the analytic
  cold-plasma hybrid dispersion relation.

### Rankine-Hugoniot and shocks (`src/initial_conditions/`, `src/verification/`, `src/boundaries/`)

- `MHDState`, `rankine_hugoniot`, and `rh_branches` MHD jump-condition solver.
- `PerpShock` perpendicular-shock setup with `init_shock!`, `step_shock!`,
  `deposit_moments!`, `compute_E!`, and `shock_density_weight`.
- Boundary injection: `flux_speed`, `flux_per_density`, `inject_face_1d!`.

### Electrostatic PIC (`src/models/electrostatic.jl`)

- `Electrostatic1D` 1D electrostatic PIC: `init_espic!`, `step_espic!`,
  `poisson_E!`, `field_energy`.

### KdV reference (`src/initial_conditions/linear_waves.jl`)

- `kdv_soliton`, `kdv_solve` analytic/numerical Korteweg-de Vries reference.

### Diagnostics (`src/diagnostics/`)

- `total_momentum`, `electric_work`, `temperatures_par_perp`,
  `velocity_histogram`, `phase_space_histogram`, `power_spectrum`,
  `pressure_strain`, `shock_front`.

### Particle sorting, smoothing, spacecraft (`src/particles/`, `src/diagnostics/`)

- `cell_index`, `sort_particles!`, `particles_per_cell`, `memory_bytes`.
- `binomial_smooth!`, `smoothing_transfer`.
- `gather_at`, `SyntheticProbe`, `sample!`, `advance!`, `shock_frame`,
  `dehoffmann_teller_velocity`, `classify_reflected`.

### I/O, metadata, normalization (`src/io/`, `src/verification/normalization.jl`)

- `save_checkpoint` / `load_checkpoint!`, `save_run` / `load_run`,
  `RunMetadata`, `capture_metadata`, `CHECKPOINT_SCHEMA_VERSION`.
- `PlasmaUnits`, `alfven_speed`, `gyrofrequency`, `inertial_length`, `to_SI`,
  `to_normalized`.

[Unreleased]: https://github.com/jake-w-liu/HybridPlasmaPIC.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/jake-w-liu/HybridPlasmaPIC.jl/releases/tag/v0.1.0
