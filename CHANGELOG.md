# Changelog

All notable changes to HybridPlasmaPIC.jl are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `EvolvingElectrons(γ = 5/3)` advances independent scalar pressure with the
  magnetic field in periodic `HybridStepper` runs, including ion collections
  and replicated MPI. Set `fields.pe` before initialization. The ideal massless
  closure requires `γ > 1` and `η = ηH = de2 = 0`; CAM-CL and Hall-MHD reject it.
  Invalid pressure stages roll back public state without clipping. Checkpoints
  and archives preserve pressure, and energy diagnostics accept it through
  `pe=`. The guide covers timestep limits, positivity, time levels and retries.

- MPI checkpoint/restart supports ordered ion collections, unequal local counts,
  empty populations, writable views and continued particle migration. Reads,
  destination validation and species/phase agreement finish on all ranks before
  restoration. Schema 2 uses generation directories and still reads schema 1.

- MPI hybrid initialization, stepping and moment deposition accept ordered ion
  collections. Shared serial kernels retain each population's charge, mass and
  identity; raw moments are reduced across ranks before normalization. Collection
  input validation is collective, and views support stepping without migration.
  The MPI benchmark now measures single- and two-ion evolution and allocations.

- Multispecies kinetic energy, total momentum and hybrid energy/momentum
  budgets sum each population's mass-weighted contribution and count shared
  field/electron reservoirs once. Typed empty collections return zero particle
  contributions; heterogeneous nonempty momentum collections remain supported.

- Hybrid multispecies checkpoint/restart and archives preserve ordered ion
  populations and their shared integrator state. Restore preflights every
  destination before mutation, supports independent writable views and empty
  populations, and retains the scalar checkpoint representation unchanged.

- Multispecies `HybridStepper` and `CAMCLStepper` initialization and evolution
  share charge-weighted moments and fields while retaining each positive ion
  population's charge, mass and storage. Supports different counts, writable
  views, empty populations, changing timesteps and magnetic subcycles. Invalid
  populations and shared storage are rejected before initialization or stepping
  mutates state. Includes independent two-ion orbit, discrete-force and wave
  convergence checks, an executable example and a runtime/allocation benchmark.

- Hall-MHD checkpoint/restart and archives preserve model parameters, fluid and
  derived fields, simulation clock, and inactive hyperresistive stage history.
  Restoration validates the complete payload before changing existing field
  storage and continues exactly with the same stepping inputs.

- Electromagnetic PIC checkpoint/restart and archives for `EMPIC1D` and general
  `EMPIC` preserve both species, physics settings, fields/current diagnostics,
  and fast/full timestep history. Both particle destinations are checked before
  restoration; independent writable views and empty populations are supported.

- General `EMPIC` supports electron/Maxwell subcycling (`n_sub`) in 1D–3D.
  Mobile ions receive one kick per full step while both species drift and deposit
  current at every fast substep. Changing the timestep or subcycle count preserves
  the separate electron, ion, and magnetic-field phases.

- Electrostatic PIC checkpoint/restart and archives for the legacy 1D and general
  1D–3D solvers preserve electric fields, density, particles, and leapfrog phase.
  Restores reject incompatible background density, grid, precision, and integrator
  before modifying state, and support writable particle views and empty populations.

- CAM-CL checkpoint/restart and run archives preserve persistent magnetic copies,
  subcycle phase, and particle timestep history. Restores validate the integrator
  and all saved field components before modifying the destination.

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
  with logged findings and a bounded regression gate, a `.JuliaFormatter.toml`
  (default style, 100-column margin), a
  formatting-check CI job, Ubuntu package jobs for the minimum supported Julia
  minor and latest stable release, real 2-/4-/8-rank MPI jobs, and uploaded test
  logs. JET is required and permits the 0.12 series for current Julia releases.
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

### Changed

- Shock boundary totals reuse exact nodal particle sums already computed in
  the same call. The checked 256-/4,096-particle exceptional cases allocate
  50.4–50.8% of the previous implementation, while ordinary outputs and
  allocations remain unchanged. Uncached boundaries retain complete evaluation.

- Exceptional electromagnetic continuity transforms extract dyadic digits with
  integer shifts and masks, preserving signed digits and transform order while
  reducing measured complete-call allocations. Ordinary outputs and allocation
  totals remain unchanged in the checked diagnostic configurations.

- The legacy EM continuity diagnostic reuses the grid's FFT plans. Checked
  results remain identical and the measured allocation falls by 544 bytes per
  call. Existing grids keep the FFT settings selected at construction.

- Hall-MHD reuses one forward FFT per velocity component when computing
  advection gradients in 2D and 3D. The reproducible benchmark reports runtime
  and the additional 384 bytes allocated per RK4 step; 1D keeps its prior path.

- Vector field gathering reuses each particle stencil across all three components,
  preserving scalar-gather results and zero steady-state allocations. The CPU
  benchmark now compares fused and scalar gathers in 1D, 2D, and 3D.

- `SpectralOperators.jl` is now resolved from the General registry instead of a
  pinned GitHub URL; the `[sources]` override was removed from `Project.toml`
  and `validation/Project.toml`. Installing HybridPlasmaPIC no longer requires
  adding the operator package by URL first.
- Raised the `SpectralOperators` compat floor to `0.1.1`. The re-exported
  `fft_friendly_size` in 0.1.0 searched upward one integer at a time, which does
  not terminate in practical time above the largest 7-smooth `Int`; 0.1.1
  enumerates `2^a·3^b·5^c·7^d` instead. 0.1.1 also rejects grids whose spacings
  or Fourier wavenumbers are not finite and representable in the element type.

### Fixed

- Shock boundary kinetic and electron enthalpy fluxes retain finite complete
  values when intermediate products overflow or underflow. Particle flux also
  retains CIC weights, nodal cancellation and the final SBP norm. Certified
  ordinary arithmetic and the singular gamma-one enthalpy convention are
  preserved. The total retains complete component expressions when their
  separate rounding or overflowing cancellation would lose its range.

- CGL electron-energy diagnostics retain complete finite density, field,
  invariant-coefficient and cell-spacing expressions through range-sensitive
  arithmetic. Certified ordinary cases keep their existing arithmetic and
  allocations; the zero-field parallel-pressure singularity remains explicit.

- Nonrelativistic particle-energy diagnostics retain complete finite scalar and
  multispecies expressions when intermediate arithmetic loses range, including
  signed cancellation between overflowing populations and hidden final overflow.
  Range checks retain the original arithmetic when they certify it as safe.
- Hybrid magnetic energy shares the electromagnetic complete-integral evaluator.
  Both electromagnetic PIC families also detect finite-looking sums that hide
  true output overflow. Finite fields and cell spacings are retained through
  final rounding when their intermediate products lose range.

- Clarify hybrid energy accounting: the Ohmic quadratic diagnostic is not
  generally the net tracked-energy loss, and the reduced inertia reservoir
  does not establish moving-ion total-energy conservation. Numerical outputs
  and forces are unchanged. The README and diagnostic docstrings use the same
  qualified conservation statements.

- Legacy electromagnetic PIC removes the transverse Nyquist mode without
  overflowing its alternating sum or prematurely rounding a subnormal mean.
  Selected finite inputs use the complete cell expression, preserving finite
  currents and electric fields in extreme public-step cases. Ordinary arithmetic,
  odd-grid behavior and nonfinite propagation retain their existing behavior.

- Both electromagnetic PIC families reconcile finite built-in-shape charge
  snapshots to actual particle and background charge using complete stored cell
  geometry. Each corrected cell is rounded once, preserving nonzero net charge
  and bounding the total discrepancy by final cell rounding. Charge snapshots
  intentionally change, so bitwise trajectories across versions are not promised.
  Raw density deposition, custom-shape coefficients and saved checkpoint arrays
  retain their existing behavior; validation tolerances are unchanged.

- Both electromagnetic PIC families combine complete species currents when
  intermediate arithmetic overflows or underflows, including the legacy
  longitudinal mean. Midpoint deposition restores borrowed particle positions
  when deposition fails and propagates the original error.

- CI groups runs by workflow, ref and commit so a new pushed revision can finish
  its checks without cancelling an earlier revision's acceptance run.

- Both electromagnetic PIC families preserve complete electron/ion charge
  cancellation and background subtraction when intermediate density arithmetic
  overflows or underflows. Ordinary arithmetic and electron-density scratch
  remain unchanged; initialization and subcycled stepping are covered.

- Per-species density, momentum and current preserve complete finite results
  across intermediate overflow/underflow, including threaded deposition and
  both EMPIC current paths. Serial and threaded callers share one stencil
  kernel; ordinary arithmetic, generic raw deposition and scratch contents
  retain their existing behavior.

- Relativistic kinetic energy retains finite complete weighted results when
  intermediate squares, products or sums exceed floating-point range, and
  resolves uncertain speed limits from the exact stored velocity components.
  Ordinary results retain their arithmetic order; signed weights, nonfinite
  factors, BigFloat precision contexts and writable views remain supported.

- Threaded scalar and density deposition share serial input validation and
  reject output aliases of values or particle storage before clearing arrays.
- KdV dealiasing projects both operands and the quadratic product onto the
  strict integer-mode band `3abs(m) < n`, preventing excluded-input and cutoff
  aliases while preserving linear evolution of the complete state.
- RayCon uses analytic Solovev magnetic-field derivatives for geometry and
  amplitude transport, reducing sensitivity to last-bit launch perturbations.
- MPI checkpoint-publication tests synchronize all readers before renaming
  fixture files. The eight-rank CI step has a bounded 40-minute window to
  finish the full suite; two- and four-rank steps retain 20 minutes.

- Run archives retain the validated byte-vector payload type through checksum
  verification and decoding, eliminating unrelated loader inference reports
  when HDF5 and Pencil extensions are loaded. Archive format and error handling
  are unchanged; concrete loader compiler checks cover String and SubString paths.

- Checkpoint and archive preflight rejects configuration aliases that Julia
  reconstructs independently, including copied views/indices, reconstructed
  container buffers and repeated BigInt/Regex objects. Independent parameters,
  shared containers and cycles remain supported. Native FFT plans and non-null
  pointers reject before publication or restore mutation, including MPI saves.

- Both EM continuity diagnostics retain finite complete residual ratios when
  density differences, transforms or timestep-weighted currents exceed the
  element type's intermediate range. Ordinary values and allocations are
  preserved. Spectral input restrictions and FFT cancellation accuracy remain
  part of the supported numerical contract.

- Electromagnetic field energy retains finite complete integrals when squared
  fields, light speed or cell-volume intermediates exceed floating-point range.
  Both EMPIC interfaces preserve ordinary summation order and expose genuinely
  unrepresentable final energies as infinity.

- Serial checkpoints, archives and asynchronous saves publish only after a
  same-directory staging file is serialized, flushed and closed. Earlier errors
  preserve the previous file; cleanup failures remain visible. Existing rwx
  modes and symlink referents are retained; new files use 0600 and other hardlinks
  retain their old contents after replacement.

- General ES/EM Poisson, Gauss and current projections distinguish nonzero
  wavevectors even when their squares lose range. Range-sensitive modes retain
  the complete longitudinal expression and original timestep before rounding;
  true DC and Nyquist behavior is preserved. Density/FFT range boundaries remain
  tracked separately in `docs/EMPIC_RANGE_PLAN.md`.

- Checkpoint and archive validation rejects model/shape links to evolving
  particle/solver storage before serialization or restore mutation, including
  owned FFT plans, clocks, packed buffers and dictionary backing allocations.
  Independent mutable parameters, tested cycles and disjoint logical views
  remain supported. Internal configuration views/native resources have separate
  serialization limitations recorded in the configuration plan.

- Relativistic velocity preflight rejects stored superluminal components that
  previously rounded below the light-speed boundary. Radiation, both EM PIC
  APIs, and electromagnetic checkpoint validation share the repair and reject
  before mutation. Ordinary validation remains allocation-free; near-boundary
  inputs use an exact comparison of the original components.

- Radiation reaction validates fields, coefficient and timestep after conversion
  to particle precision, including no-op calls, before changing particles.
  Shared uniform-field validation constructs error labels only on failure,
  retaining allocation-free ordinary calls. Extreme-range radiation arithmetic
  remains under repair; the documentation now states its limitation.

- MPI checkpoint saves validate all ranks and publish a new manifest only after
  all new rank files close. Recoverable failures before publication preserve the
  previous checkpoint. Cleanup removes only owned files, rejects generation
  symlinks and reports failures without hiding the original error.

- Gyro averaging retains finite ring coordinates and sample means across large
  sample counts, extreme field geometry and cancellation. Callbacks run once per
  point in phase order; invalid coordinates are rejected before their callback.
  Cached lossless alignment reduces measured default 16-point costs by 18–37%
  relative to the correctness repair, retaining zero ordinary allocations.
  Small-ring and fallback timings remain mixed; compilation costs increase.

- Both electromagnetic PIC solvers preserve finite electric-field updates when
  squared light speed, magnetic force or timestep multiplication loses range.
  An exact fallback evaluates the complete field expression, including the
  specialized longitudinal update. Ordinary arithmetic and allocations remain
  unchanged. Range-sensitive spectral curls use scaled derivative bands through
  the complete electric or magnetic update; magnetic priming, centering and
  advancement also retain timesteps whose half interval would underflow.

- Moment deposition validates every output, scratch buffer and particle layout
  before writing. Momentum, current, pressure and single/multispecies hybrid
  moments reject overlapping storage without corrupting inputs or outputs. MPI
  moment validation completes collectively before deposition; independent views
  and sequential scratch reuse remain supported.

- Public E×B, grad-B and curvature drifts retain raw cross-product terms
  through their physical coefficients, preserving finite results across range
  loss and near-parallel cancellation. Exact fallbacks share the existing
  conversion logic and retain stored BigFloat inputs and rounding modes.

- Shared checkpoint validation rejects nonfinite or negative simulation time
  before saving or restoring hybrid, CAM-CL, electromagnetic PIC and Hall-MHD
  state. Valid time origins, signed zero and terminal step counters remain valid.

- General EMPIC charge-continuity diagnostics retain the true zero-frequency
  mode, so changes in total charge contribute to the residual. Pure Nyquist
  combinations outside the spectral derivative range remain excluded.

- Checkpoint destinations reject overlapping particle components before any
  resize or write, including empty and reinterpreted storage. Independent
  writable views remain supported. Hybrid collection validation shares the
  same allocation-free component check.

- MPI checkpoint restore completes reads and destination validation on every
  rank before changing any destination. Recoverable file and validation errors
  retain their original exception on the failing rank and notify peer ranks.

- Collision rollback tests initialize their position fixture, avoiding
  nondeterministic NaN comparisons in the full test suite.

- Guiding-centre position updates retain the complete streaming and drift
  displacement through the final coordinate addition. Finite updates survive
  intermediate range loss and cancellation; exact fallbacks preserve stored
  BigFloat inputs and round hardware outputs independently of BigFloat settings.
  Rejected updates leave the guiding-centre state unchanged.

- AMR exact fallbacks round hardware floats directly from the complete rational
  result, avoiding a second rounding through ambient BigFloat precision. Shared
  exact conversion also retains stored BigFloat inputs and the caller's BigFloat
  rounding mode. Tests use independent rounding-cell and dyadic oracles.

- Particle storage checks recognize shared empty vectors, reinterpret wrappers,
  and views with shared mutable indices. Checkpoint restoration, hybrid
  collection evolution and ionization reject these coupled buffers before
  mutation, while independent empty storage and disjoint nonempty views remain
  supported.

- AMR restriction preserves cancellation when BigFloat child values retain
  more precision than the caller's current setting. The complete average is
  rounded once at the working precision; hardware-float paths are unchanged.

- AMR prolongation rounds each complete reconstructed child value, preserving
  small results when a large slope cancels the cell mean. Exact fallback also
  handles the smallest-normal boundary and retained BigFloat input precision;
  ordinary measured hardware-float profiles remain allocation-free.

- AMR cell centers retain representable coordinates across intermediate overflow
  and large half-cell indices, including Float16 rounding boundaries and stored
  BigFloat precision.
- AMR refinement uses power-of-two scaling and certified rounding for its
  second difference, preserving zero curvature on affine profiles and tiny
  representable curvature near cancellation. Uncertain results use exact
  arithmetic; BigFloat inputs retain their stored precision and respect the
  caller's rounding mode. Ordinary Float32/64 reductions remain allocation-free.

- Guiding-centre parallel velocity evaluates the combined electric, mirror,
  timestep, mass and initial-velocity expression without avoidable intermediate
  overflow or underflow. Exact arithmetic resolves uncertain cancellation;
  certified exact equilibria keep an allocation-free path.

- Hall-MHD rejects a nonfinite final electric field before committing the
  candidate step. This includes overflow created by the final magnetic
  projection; rejected steps preserve public fields and clock.

- Hall-MHD accepts density exactly equal to `nfloor`, consistently with the
  shared Ohm-law floor. Densities below the floor still reject the step before
  changing public fields or the simulation clock.

- Checkpoint capture and restoration retain concrete field types for static
  analysis, including CAM-CL magnetic copies and EM species presence. The JET
  regression gate passes on Julia 1.12 and 1.13 without raising its limit;
  serialized checkpoint bytes and restoration behavior are unchanged.

- Hybrid, CAM-CL, Hall-MHD, electromagnetic PIC, and MPI reference steps reject
  time and step-counter overflow before changing state. EM clock validation
  preserves the existing sequence of rounded fast-substep additions.

- Hall-MHD completes final magnetic projection and Ohm-law evaluation on its
  candidate state before committing. A failed final evaluation preserves public
  fields, time, and step count; successful steps retain public array identities.

- Guiding-centre grad-B drift uses the stored magnetic moment directly, preserving
  finite drifts when reconstructing the perpendicular speed squared would overflow
  or underflow.

- Timestep recommendations retain representable whistler, resistive, and gyro
  limits when intermediate powers or products exceed floating-point range.

- CI tests both Julia 1.12 and the latest stable Julia release. Test failures
  propagate through the log-capture pipeline, and each version retains its own log.

- SHA compatibility accepts both the 0.7 and 1.x standard-library versions,
  including SHA 1.0 shipped with Julia 1.13.

- Electromagnetic PIC tracks particle and magnetic-field half-time levels when
  timesteps change, including separate electron/ion intervals in subcycled
  `EMPIC1D`. Invalid subcycle counts and unrepresentable substeps fail before mutation.

- Hybrid, CAM-CL, and replicated MPI stepping retain particle half-time levels
  across timestep changes. Hybrid checkpoints and archives now preserve timestep
  history and model/shape parameters, restoring the saved physics. Writable
  particle views restore in place; unsupported resizing is rejected before mutation.

- HDF5 field writes preserve logical values for ranges, transposes, and strided
  views. Contiguous arrays and views retain the direct write path.

- Electrostatic PIC retains the previous timestep so variable-step particle kicks
  span the correct half-time levels; fixed-step results remain unchanged.

- Field-file reads validate header extents and available payload before allocating
  arrays, rejecting truncated and overflowing dimensions with bounded memory.

- Run archives hash the exact serialized state bytes, so dictionary/set storage
  reconstruction no longer causes false checksum failures. Schema 4 verifies
  integrity before decoding state; older schemas are rejected explicitly.

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
