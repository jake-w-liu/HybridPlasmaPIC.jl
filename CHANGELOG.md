# Changelog

All notable changes to HybridPlasmaPIC.jl are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
