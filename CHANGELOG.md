# Changelog

All notable changes to HybridPlasmaPIC.jl are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Entries describe source changes. The [audit tracker](docs/DEEP_DEBUG_TRACKER.md)
records scoped verification and remaining exact-revision acceptance gates.

### Added

- `CartesianAMRHierarchy` and `AdaptiveCompositeElectrostaticPIC` add coupled
  CPU 2D3V/3D3V grounded-box evolution: balanced dyadic patches, exact hanging
  Q1 constraints, full tensor physical clouds, cold matrix-free field solve,
  KDK and signed charge/field/kick/trajectory work. Real multilevel adaptation
  records mesh-energy exchange and preserves prior physical-step history on its
  own mesh. Complete topology/count-aware raw/archive persistence and compact
  independent references are registered. Source review is complete; all new
  numerical/compiler/format/scientific/resource acceptance remains pending.

- `FieldAlignedConductingElectrons` adds varying-B tensor heat transport to
  canonical scalar-pressure stages, with diagonal Nyquist complements, two
  bounded gradient/flux preparations and complete signed-source publication.
  Instantaneous diagnostics privately redeposit density and measure the channel
  adjoint defect. Policy 3 preserves full carry and incoming direction admission;
  plain and collisional archive readers delegate to canonical restore. Source
  and essential independent controls are reviewed and registered; all new
  runtime, scientific and resource acceptance remains pending.

- Canonical `ImplicitHybridStepper` now jointly solves ion velocity, B and
  evolving scalar electron pressure, with fixed third residual scale, exact final
  pressure gate, nonnegative admission, synchronous endpoint pe and complete
  six-term stored-operand work receipt. Isothermal 19-field report/schema1 stay;
  evolving schema2 preserves original coefficient identity and full phase. The
  existing raw/archive owners now also provide the missing implicit-hybrid
  archive reader. Source and compact essential controls are reviewed; runtime,
  historical-schema continuation, compiler/format and full CI remain pending.

- `ToroidalAnnulusHybridPIC` supplies coupled CPU 3D Float32/64 ions, magnetic
  flux and adiabatic electron pressure on a positive-volume toroidal annulus.
  Every RK4 stage refreshes real moments and fields; electrical contacts and
  incoming thermal reservoirs retain actual signed work. Canonical CIC/Lorentz/
  pressure/RK owners, detached diagnostics and complete guarded raw/archive
  continuation are integrated. Source review is complete; essential remote
  controls, full CI and scientific/resource acceptance remain pending. Solid-axis,
  material-wall, equilibrium/import and backend stages remain open.

- The optional resident `CUDAHybridStepper` family implements a single-device
  D1/3V Float64 isothermal hybrid lifecycle with shared CIC/mover/Ohm/RK/Fourier
  owners, magnetic subcycles, complete predictor/phase diagnostics and distinct
  full physical checkpoint/archive persistence. CUDA 6.4.1 within 6.x is required;
  the source is reviewed, while actual NVIDIA numerical/resource acceptance is open.

- Actual physical-pencil ES/EM/isothermal-hybrid and vacuum-Maxwell lifecycle
  checks cover charged nonuniform evolution, rank ownership, one-rank rejection,
  complete matching-layout restart and variable-step continuation. Independent
  CIC density and Maxwell KDK references supplement full carried-state checks.
  Minimum/latest Julia two-rank GitHub acceptance remains pending.

- `DynamicCompositeElectrostaticPIC` adds completed-boundary 1D CPU refine/coarsen
  and base-field curvature adaptation over canonical physical-cloud Neumann/KDK
  epochs. Regrids preserve particles/global phase and record actual mesh-energy
  exchange. Diagnostics and distinct complete topology-aware checkpoint/archive
  persistence retain the prior physical step on its own mesh. Source is reviewed;
  independent regrid/Verlet/history/continuation controls pass on minimum/latest
  Julia in scoped run 42553061121; full integrated CI remains required.

- Independent compatible-Yee transverse references compare Float32/64 vacuum
  Maxwell updates, including the true Yee Nyquist mode, and magnetized unequal-q/m
  cold-plasma midpoint updates against exact rational Cayley maps. They are
  registered in the full suite and pass focused minimum/latest GitHub checks
  in run 42552199520; full integrated CI remains required.

- `ImplicitYeeEMPICStepper` adds a separate synchronous CPU 1D3V all-mobile,
  exact-neutral periodic Yee cloud midpoint Vlasov–Maxwell family on exactly
  tiled periodic grids. Every trial rebuilds exact endpoint transport, charge
  and the same-coefficient orbit current/gather from one velocity/E/By/Bz
  unknown, publishes a Gauss-consistent endpoint Ex from the exact endpoint
  charge and solved mean, and shares the bounded Newton/Jv/GMRES/Armijo driver.
  Constant Bx, exact work/defect ledgers, 256eps constraint gates, diagnostics
  and full checkpoint/archive persistence are included. See the
  [model](docs/IMPLICIT_YEE_EM.md).

- Canonical particle, coupling, moment, Ohm and RK4 scalar rows now accept an
  explicit device arithmetic policy. Default CPU entries retain their existing
  operation order, validation and exact repairs. The CUDA arithmetic policy
  remains opt-in and has separate pending device acceptance. Source review
  preserves mixed-precision custom shape weights. Ten whole-file inverses and
  all 48 Boris node expansions
  were checked. Runtime/compiler/allocation and device acceptance remain pending;
  this dependency does not itself provide a resident solver.

- `StaticAMRHierarchy1D` and `StaticCompositeElectrostaticPIC` add a bounded CPU
  1D nested 2:1 composite electrostatic family: unique-leaf/shared-node ownership,
  one mixed Neumann solve, physical clouds, synchronous Verlet and coarse/fine
  crossings. Exact charge/current/work/rounding diagnostics and full physical IO
  are included. The shared scalar-source helper moves byte-identically; existing
  cloud/force/plate owners remain canonical. Review corrected phase admission
  and pre-publication preparation. Essential remote regressions are prepared;
  runtime, self-force/resource and dynamic AMR acceptance remain open. See the
  [model](docs/STATIC_COMPOSITE_ELECTROSTATIC.md).

- `ImplicitEMPICStepper` adds a separate synchronous CPU D1–3 all-mobile,
  exact-neutral CIC Vlasov–Maxwell family on odd periodic grids. Every joint
  velocity/E/B residual and probe re-deposits endpoint charge and midpoint
  current, applies canonical continuity correction and recomputes both curls.
  Hybrid and EM share one bounded Newton/Jv/GMRES/Armijo driver; three canonical
  owner extractions retain whole-file literal inverses. Exact signed work,
  actual Gauss/divB/continuity gates and full physical checkpoint/archive are
  included. Review added rejection of impossible last_dt>time history.
  Essential remote regressions are prepared; runtime/scientific/resource
  acceptance remains open. See the [model](docs/IMPLICIT_EM.md).

- `PeriodicShearGrid` and `ShearedPeriodicHybridPIC` add a bounded CPU 3D
  periodic nonorthogonal full-f hybrid family: logical positions, Cartesian
  velocities, contravariant magnetic flux, canonical TSC/Ohm/Fourier owners,
  genuine continuous Lorentz derivatives, four collective RK4 stages and
  endpoint refresh. Source-risk rejects before mapped Faraday differentiation.
  Quiet loading, physical work/energy/flux diagnostics and complete synchronous
  raw/checksummed archive restore are included. Source review added the accepted
  receipt invariant divergence<=tolerance. Essential remote checks are prepared;
  runtime/scientific/resource and broader toroidal geometry remain open. See
  the [model and limits](docs/SHEARED_PERIODIC_HYBRID.md).

- `OpenHallHybrid` adds a separate synchronous CPU D2/3 isothermal positive-ion
  numerical reservoir family with full Hall induction, actual SBP-volume CIC,
  constrained projected RK4 trials, original-factor Boris velocities and two-face
  births/exits. Bounded stock RNG requests, stable IDs, last/cumulative exchanges,
  signed field/SAT/projection/energy ledgers and complete physical/RNG
  checkpoint/archive prepare before publication. Six canonical seam extractions
  retain literal legacy source inverses. Source review is complete and essentials
  are prepared for GitHub; runtime, physical-boundary, compiler and resource
  acceptance remain open. See the [model](docs/OPEN_HALL_HYBRID.md) and
  [persistence contract](docs/OPEN_HALL_HYBRID_CHECKPOINT.md).

- `SlabGyrokineticPIC` and true `GyrocentreSet` populations add a bounded CPU
  full-f electrostatic slab family: reciprocal ring-TSC coupling, full-k
  polarization and modified adiabatic electrons, collective four-stage RK4
  with endpoint field refresh, physical energy diagnostics and complete
  synchronous checkpoint/archive restoration. Publication stages privately;
  periodic wrapping reuses the canonical half-open coordinate owner. Source
  review is complete; essential regressions are prepared for GitHub, with no
  local Julia execution. Runtime, scientific and resource acceptance remain
  open. See the [model and limits](docs/SLAB_GYROKINETIC.md).

- `energy_budget.total` now composes original finite built-in reservoir
  expressions when normal-range/sign bounds cannot certify the associated
  native sum. Component values stay unchanged; the complete total can be finite
  when separately rounded components overflow. Existing computed J, exact
  particle/field/log/power/CGL owners and validated stored-pressure gamma are
  reused. Coarse polytropic bounds precede optional sharp refinement. This is a
  source-work reduction, without a measured speed claim; historical V6 timing
  failures and fixed gates remain intact. No tests added/run; numerical,
  compiler and cost acceptance remain open. See the
  [current source contract](docs/HYBRID_DIAGNOSTIC_RANGE_PLAN.md#current-complete-budget-source-2026-11-02).

- `MixedShockOperators` adds concrete D2/3 SBP-x/periodic-transverse gradient,
  curl, divergence and represented Euclidean transpose with cached Fourier
  work. `mixed_hall_rhs!` reuses canonical isothermal pressure/Ohm/Faraday owners
  and reports exact stored-source face work. Explicit constrained field/rate
  projection uses structural normal-face elimination, cold MINRES-QLP and fresh
  equation/constraint/work gates before publication. The reduced shock curl
  shares its geometry with a whole-file source inverse; its SAT divergence
  qualification is corrected. No advancing Hall shock or new checkpoint is
  supplied by this prerequisite. No tests added/run; native symmetry,
  convergence, compiler and cost acceptance remain open. See the
  [scope and projection contract](docs/MIXED_HALL_OPERATORS.md).

- `ImplicitHybridStepper` adds a separate synchronous CPU periodic isothermal
  collection family. Joint original-factor Boris-map/Faraday residuals feed
  bounded Newton–GMRES with fresh inner/outer gates and explicit typed failure;
  complete endpoints and represented pressure/continuity/work-defect reports
  prepare before publication. Distinct complete checkpoint/archive preserves
  full public fields, particles, fixed policy, clock and report. Canonical
  physical owners are reused; shared GMRES storage leaves the old whole ShockCN
  source recoverable by literal expansion. No tests added or run; convergence,
  accuracy, restart, compiler and resource acceptance remain open. See the
  [scope and equations](docs/IMPLICIT_HYBRID.md).

- `distributed_collisional_empic` composes the canonical actual-pencil EM
  candidate with local Coulomb scattering at new positions/carried-half
  velocities. All eleven fields and preceding drift current persist unchanged;
  one matching-layout checkpoint generation retains complete physical/RNG
  state. Cell, graph, RNG, busy and prepared-generation owners are shared with
  collisional ES; bare ES/EM/CES source paths retain whole-file inverses.
  Fixed sub-CFL duration, explicit normalization and cell ownership are required.
  No tests added or run; runtime/statistical/continuation/compiler/resource
  acceptance remains open. See the [guide](docs/DISTRIBUTED_COLLISIONAL_EMPIC.md).

- Actual-pencil isothermal hybrid gains complete matching-layout checkpoint
  source: eighteen carried fields, ordered particles and actual half-velocity/
  integer-field phase persist in distinct checksummed generations. Incoming
  state detaches and fully prepares on all ranks before publication; stale
  destination values/counts are replaceable. No deposit, Ohm, predictor or phase
  reconstruction runs during restore. Existing numerical and file mechanics
  remain unchanged. No tests added or run; runtime continuation, liveness,
  compiler and resource acceptance remain open. See the
  [checkpoint guide](docs/DISTRIBUTED_HYBRID_CHECKPOINT.md).

- `distributed_hybrid` adds actual-pencil positive-ion CIC collections with
  massless isothermal electrons, strict density, canonical prime/Boris/midpoint
  motion, frozen-moment RK4 subcycles and carried endpoint prediction. All
  eighteen fields, particles, counts and phase stage before publication.
  Diagnostics combine exact scalar mechanical sums and local log enclosures;
  no global field/log dictionary is gathered. Existing numerical owners retain
  literal whole-file inverses. Original-source differentiated Ohm range recovery
  is explicitly unsupported. Complete checkpoint is recorded separately.
  Source-reviewed only; no tests added or run. Actual-rank numerical, restart,
  compiler and resource acceptance remain open. See the
  [hybrid guide](docs/DISTRIBUTED_HYBRID_IMPLEMENTATION.md).

- Actual-pencil ES PIC gains fixed-step local Coulomb composition and complete
  matching-layout physical/full-RNG checkpoint. Owned logical cells use original
  global dx and canonical collision scheduling; an explicit cell-owner gate
  precedes publication. Per-rank stock RNGs stage with particles/fields and the
  sole inner clock/source ledger; metadata cannot alias writable recurrence
  storage. Zero logs skip both RNG copies. Canonical Coulomb/graph/ES seams are
  shared with source-preservation checks. No tests added or run; actual-rank,
  statistical/restart/compiler/resource acceptance remains deferred. See the
  [distributed collision guide](docs/DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md).

- Actual-pencil EM PIC gains distinct matching-layout checkpoint/restart with
  checksummed rank generations and complete eleven-field/current/endpoint phase
  history. All ranks validate detached source state and prepare same-type leaves
  plus derived counts before publication; stale destination values/counts can
  be replaced. Existing ES/Maxwell/shared file mechanics stay unchanged. No tests
  added or run; continuation/liveness/compiler/resource acceptance is deferred.
  See the [EM checkpoint guide](docs/DISTRIBUTED_EMPIC_CHECKPOINT.md).

- Actual-pencil electromagnetic PIC adds all-mobile nonrelativistic CIC
  collections with canonical original-factor Boris/prime/midpoint motion,
  complete remote transport and migration, represented-mode current correction,
  carried magnetic/velocity phase and staged eleven-field publication. Exact
  global diagnostics separate physical/grid charge and phase. Shared serial
  expressions and ES transport owners retain literal source inverses. No tests
  added or run; runtime/collective/compiler/cost acceptance is deferred and
  complete EM checkpoint is recorded separately below. See the
  [EM source guide](docs/DISTRIBUTED_EMPIC_IMPLEMENTATION.md).

- `YeePML` adds source-free causal directional ADE PML for serial D1–3 CPU
  Float32/64, separate target-lattice memories, original-duration midpoint KDK
  and signed exchange/field-defect ledgers. Distinct checkpoint/archive retains
  actual advanced memories; zero-prehistory init is explicit. Canonical Yee
  spatial/field owners are reused. Source-reviewed only; no tests added or run.
  Reflection/stability, runtime/restart/cost, general CPML and particle/source
  integration remain open. See the [PML guide](docs/YEE_PML.md).

- Actual-pencil electrostatic PIC gains initialized-boundary matching-layout
  checkpoint/restart. A distinct checksummed generation retains complete local
  particles, density/E and variable-step phase; all ranks validate detached
  sources and prepare same-type replacements before publication. Stale destination
  values and changed local counts are supported by the source contract. Seven
  file/checksum helpers are shared with Maxwell while its schema/orchestration
  stays unchanged. No tests were added or run; restart/liveness/compiler/resource
  acceptance remains deferred. See the [checkpoint guide](docs/DISTRIBUTED_ELECTROSTATIC_CHECKPOINT.md).

- `YeeConductiveLayer` adds a graded scalar electric/magnetic conductivity layer
  for D1–3 CPU Float32/64 Yee vacuum fields. Fixed-step midpoint damping surrounds
  the unchanged KDK candidate, with exact actual loss, source work and signed
  physical-field defect ledgers staged before publication. Distinct nested
  checkpoint/archive preserves fields and budget; derived caches rebuild.
  This is not PML. Source-reviewed only; no tests were added or run. Reflection,
  stability, restart/compiler/resource acceptance and particle coupling remain
  open. See the [layer guide](docs/YEE_CONDUCTIVE_LAYER.md).

- Generic D1–3 `EMPIC` collections gain explicit dividing species kick schedules
  and `set_empic_schedule!`. They share the specialized schedule/history owner,
  retain every-fast drift/current and full-vector physics, and stage histories
  with the complete operation. Default synchronous state/schema 1 is preserved;
  opted-in state has complete schema-2 restart and validated schema-1 migration.
  Collisional EMPIC explicitly retains its unscheduled contract. Source-reviewed
  only; no tests were added or run and runtime/resource acceptance is deferred.
  See the [schedule guide](docs/EMPIC_SPECIES_IMPLEMENTATION.md).

- Distributed ES diagnostics now include global three-component particle momentum.
  Original mass/weight/velocity products and rank sums remain exact until one final
  component rounding. Momentum retains the stored velocity half phase; no new
  conservation or measured cost claim is made. No tests or executions were added.

- `distributed_electrostatic_pic` adds a rank-owned periodic 3D CPU Float32/64
  CIC collection using actual physical pencil ownership, complete remote stencil
  transport, direct particle migration and the canonical leapfrog mover. Exact
  original-factor charge contributions route before density rounding; explicit
  `grid_dc=:project` and separate physical/grid charge ledgers define neutrality.
  Particles, density/E and clock publish after collective success. No tests were
  added or run; runtime/compiler/resource acceptance is deferred. Matching-layout
  checkpoint source is now supplied separately. See the [guide](docs/DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md).

- `CollisionalElectrostaticPIC` composes generic CPU Float32/64 CIC species
  collections with canonical local Coulomb scattering at new positions and
  carried half-time velocities. It owns the sole physical/source clock and
  stages fields, particles and RNG before publication. A distinct full-RNG
  checkpoint/archive retains the nested ES state. Fixed dt, zero background,
  required neutrality and equal positive permittivity/count conversion define
  the first scope. No tests were added or run; runtime/compiler/resource
  acceptance is deferred. See the [guide](docs/COLLISIONAL_ELECTROSTATIC.md).

- Serial `MaxwellVacuum` gains field-only checkpoint/archive overloads for D1–3
  CPU Float32/64 periodic/PEC Yee grids. A distinct schema records synchronized
  fields, exact configuration and clock; complete incoming validation/detachment
  precedes in-place copies, allowing stale destination values. Existing field
  arithmetic is retained. No tests were added or run; restart/compiler/resource
  acceptance is deferred. See the [guide](docs/MAXWELL_VACUUM_CHECKPOINT.md).

- `distributed_poisson_field!` adds a static rank-owned periodic 3D CPU Float32/64
  electrostatic field solve using the canonical serial Poisson mode and existing
  distributed transforms. Exact represented neutrality or an explicit compensator
  is required; three fields publish after collective success, with an exact charge
  ledger. Scalar charge reduction shares the existing Maxwell energy transport
  owner. No field is gathered. Source-reviewed only; runtime/compiler/resource
  acceptance remains open. The distributed ES model now composes this field owner;
  see the [field guide](docs/DISTRIBUTED_POISSON_FIELD.md).

- `CollisionalEMPICStepper` composes all-mobile CPU Float32/64 generic EMPIC
  collections with canonical local Coulomb scattering at the new positions and
  carried half-time velocities. Fields, particles, source clock and owned RNG
  stage before publication; a distinct schema preserves the complete RNG and
  nested field/particle state. Fixed dt, nonrelativistic motion, `n_sub=1`, zero
  background, required neutrality and equal positive collision permittivity/count
  conversion define this first scope. No tests were added or run; runtime and
  resource acceptance are deferred. See the [guide](docs/COLLISIONAL_EMPIC.md).

- Physical-pencil vacuum Maxwell state gains `save_mpi_checkpoint` and
  `load_mpi_checkpoint!` overloads with a distinct matching-layout schema, streamed
  rank-file SHA-256 digests, a checksummed manifest, guarded generation publication
  and detached all-rank restore preflight. Plans and communicators are omitted;
  changed-layout restart and archives remain separate work. Source-reviewed only;
  restart/liveness/compiler/resource acceptance is deferred. See the
  [checkpoint guide](docs/DISTRIBUTED_MAXWELL_CHECKPOINT.md).

- `distributed_maxwell_vacuum` adds rank-owned periodic 3D CPU Float32/64 field
  evolution over the existing physical PencilArrays and collective curl owners.
  Synchronized KDK uses complete scalar affine updates, a strict downward spectral
  CFL bound and staged six-field publication. Global energy exchanges exact scalar
  contributions; diagnostics retain divergence residuals. No full field is
  gathered. Runtime/compiler/resource acceptance is deferred during the current
  implementation pass; matching-layout restart is source-implemented, while
  coupled particle evolution remains next work. See the [field guide](docs/DISTRIBUTED_MAXWELL_IMPLEMENTATION.md).

- `YeeAbsorbingLongitudinalPIC` adds a nonrelativistic CPU Float32/64 finite-cloud
  source candidate on a closed longitudinal 1D Yee grid. It admits zero captures
  or one uniquely isolated certified-passive instantaneous capture per external
  step, using KDK with original forces and rational local durations. Full death
  current, dynamic exact free plates and separate capture/work/energy ledgers
  accompany stable Vector survivor compaction and distinct checkpoint/archive.
  Survivor carry handles the force jump separately from the reflection contract.
  GitHub numerical/convergence, rollback/restart, JET/compiler, formatting/CI
  and complete-step resource acceptance are deferred during the current
  implementation-first pass; existing authored checks remain unrun. General
  material/sheath laws, multiple captures, transverse walls and backends remain
  open; see the [model guide](docs/YEE_ABSORBING_LONGITUDINAL.md).

- `YeeImpactLongitudinalPIC` adds a nonrelativistic CPU Float32/64 finite-cloud
  reflection source candidate on a closed 1D Yee grid. Each external step admits
  zero contacts or one connected certified cluster, with original-factor forces,
  rational-duration KDK chunks and conservative carried physical remainder checks.
  Wall ledgers, diagnostics and distinct checkpoint/archive state stage before
  publication. Exact-revision GitHub numerical, convergence, rollback/restart,
  JET/compiler, formatting/CI and complete-step cost acceptance remain pending.
  Broader wall laws, multiple clusters, transverse walls and backend integration
  remain open; see the [model guide](docs/YEE_IMPACT_LONGITUDINAL.md).

- Private canonical Yee raw-sum/candidate-conversion, exact source-row and
  unrounded source-work owners support the reflection candidate without copying
  deposition or force kernels. Existing callers retain their expressions,
  signed-zero behavior and public routing; GitHub verification remains pending.

- Specialized `EMPIC1D` collections add dividing per-species kick strides,
  every-fast drift/current and independent staged velocity histories. The atomic
  `set_empic1d_schedule!` setter preserves old phase; diagnostics report nominal
  offsets and copied schedule/history. Schema 2/policy 1 explicitly migrates
  validated S1 state. Essential remote numerical/compiler/format and cost gates
  remain pending; see the [S1/S2 tracker](docs/EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md).

- `EMPIC1D(g, species)` adds synchronous CPU collections in the specialized
  Ex/Ey/Bz planar polarization, with optional relativistic motion, physical
  diagnostics, complete staged publication and distinct collection checkpoint/archive.
  Legacy scalar/two-species numerical owners are shared; exact-revision GitHub
  runtime/compiler/format and cost acceptance remain pending in the
  [S1/S2 tracker](docs/EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md).

- The existing manual Yee EM essentials suite records warmed complete-step
  timing, allocation traffic and GC metadata for small ordinary and original-force
  range profiles. Measurements run on GitHub; peak memory and performance
  acceptance remain open in the [bounded plan](docs/BOUNDED_PARTICLE_COUPLING_PLAN.md#complete-step-cost-observation).

- `YeeElectromagneticPIC` supplies a bounded nonrelativistic CPU 1D3V cloud
  source candidate with both transverse polarizations, compatible path currents,
  synchronized Boris/Yee stages, scoped rest/impact guards, diagnostics and a
  distinct checkpoint/archive schema. Exact-revision GitHub acceptance is pending;
  see the [contract](docs/BOUNDED_PARTICLE_COUPLING_PLAN.md#event-free-1d3v-electromagnetic-source-candidate).

- Fixed-step `CollisionalHybridStepper` composes `ConductingElectrons` with
  ion-ion scattering through its existing field/pressure and force/RNG owners.
  The inner conduction policy is required at restart; quotient scratch must
  be independent of the outer force snapshot. The [composition guide](docs/COLLISIONAL_HYBRID.md)
  records the first-order scope and pending GitHub acceptance.

- `ConductingElectrons(γ; κ)` adds constant isotropic heat flux from original
  stage `Q(pe)/Q(n)` through bounded two-term preparation, the full Nyquist
  Fourier operator and exact per-axis line-mean projection. The complete signed
  pressure rate retains its range and rollback guards. Required conduction
  policy 2 rejects policy-1 checkpoints; FFT accuracy and complete-step cost
  acceptance remain pending. The serial periodic CPU scope is in the
  [electron guide](docs/EVOLVING_ELECTRONS.md#constant-isotropic-conduction-source-candidate).

- `YeeLongitudinalPIC` supplies event-free CPU 1D longitudinal electrostatic
  coupling with an explicit physical-width cloud, synchronized Verlet stages,
  original-source background/free-plate initialization, energy/Gauss diagnostics
  and distinct checkpoint/archive state. Strict impact and rest-mode guards
  reject unsupported steps. Initialization builds Ex from exact cumulative
  nodal charge before final rounding. Broader bounded EM and wall-model
  acceptance remain open in the [contract](docs/BOUNDED_PARTICLE_COUPLING_PLAN.md).

- Optional `YeeTriangularCloud(a)` extends closed 1D Yee static charge, endpoint
  current and path voltage with a finite physical-width cloud and explicit
  particle-bound endpoint tails. `yee_longitudinal_field` requires that shape
  and returns its stationary work-adjoint field. Exact stored-input polynomials
  retain the original width precision and stored-dx work metric; point/CIC
  defaults remain. The cloud has a declared distinct wall response; mechanical
  wall dynamics belong to separate coupled source models. Full bounded EM
  remains open.

- A manual GitHub workflow runs existing implementation essentials on minimum
  and latest Julia, supplementing the unchanged full CI matrix. Runtime checks
  remain remote; no additional numerical tests are introduced by this workflow.

- Periodic CPU `EMPIC` collections support arbitrary finite charge/positive mass,
  mobile and fixed populations, one charge-conserving current/Maxwell update,
  synchronous fast velocity phases, explicit neutrality/background policy,
  optional relativistic motion, physical relativistic particle diagnostics,
  shared-field energy and distinct schema-1 checkpoint/archive. Saved mode and c
  restore together with their validated particles. Existing electron/ion and
  `EMPIC1D` APIs remain. Scoped verification, species multirate schedules and
  remaining intersections are recorded in the
  [implementation tracker](docs/EMPIC_SPECIES_IMPLEMENTATION.md).

- `total_momentum_relativistic` and collection `kinetic_energy_relativistic`
  complete physical moments across finite CPU Float32/Float64 populations with
  original mass, weight and velocity operands. Existing residual balls certify
  final rounding; unresolved cases use the existing exact radical owner.

- Fixed-step `CollisionalHybridStepper` accepts constant Ohmic/hyperresistive
  `TransportElectrons` with `η >= 0`, `ηH >= 0`, `de2 = 0`. An optional outer force snapshot
  restores the independently carried ion force alongside pressure and RNG;
  schema-1 restart retains the canonical required transport force and policy.
  Essential composition, malformed-state and continuation checks await GitHub
  CI. Field stages own the heating; ion-ion collisions leave pressure fixed.

- Fixed-step `CollisionalHybridStepper` accepts ideal `EvolvingElectrons`.
  The existing field stages evolve electron pressure, ion-ion collisions hold
  it fixed, and pressure and RNG share whole-step rollback and restart. The
  schema-1 reader checks agreement between its two stored pressure copies.
  Essential composition and continuation checks await GitHub CI; see the
  [scoped contract](docs/COLLISIONAL_HYBRID.md).

### Changed

- The full test suite is sharded into eight contiguous groups selected by
  `HPIC_TEST_GROUP`; CI runs a group x Julia matrix while an unset variable
  still runs everything in order.

- Fixed `EMPIC` and `EMPIC1D` collection populations share original stored-zero
  velocity admission, independent of arithmetic zero comparison under FTZ.
  Signed zeros remain accepted; essential lifecycle and restore checks await
  GitHub. This changes admission only, with numerical kernels and schemas retained.

- Closed Yee models retain concrete candidates with ordinary vector storage
  through the shared PIC staging owner. Public particle views and restart destinations
  retain their identities. Compiler and complete-step cost verification remains
  pending on GitHub.

- The event-free longitudinal and 1D3V EM source candidates retain complete
  original-geometry gathers, original q/m and the exact full-dt half through
  canonical Boris, rounding only each final half-kick velocity. Intermediate
  gathers can lie outside the grid type’s range when the completed velocity is
  finite. This changes ordinary velocity bits and requires model checkpoint
  schema 2 with an explicit force policy; schema 1 is rejected. Static IEEE-valued
  gathers and archive container schema 4 retain their contracts. GitHub
  runtime/compiler/format and complete-step resource acceptance remain pending.

### Fixed

- Boris native certification normalizes original quotient equality before its
  FMA residual test, preserving inexact ratios whose residual would underflow.
  The DD tier retains the original native range and arithmetic-mode requirements,
  checks its coefficient/norm domain, and uses the canonical full low-word
  outward distance. Exact-zero multiplication treats either operand symmetrically.
  Physical formulas, exact replay and four-ulp tolerances are unchanged; remote
  regression/compiler/resource acceptance remains required.

- `ImplicitEMPICStepper` now publishes divergence-consistent endpoint fields
  inside every residual: trial E is projected onto the endpoint charge and
  trial B is projected divergence-free before midpoint fields, gathers, curls
  and constraint gates are evaluated. Previously the endpoint Gauss defect and
  magnetic divergence equalled the divergence of the residual, bounded only by
  the nonlinear tolerance, so the 256eps gates rejected converged nonuniform
  Float64 steps in 1D–3D. Newton rows and the Jacobian are unchanged.

- Completed-clock admission checks the preceding duration at the existing live
  and checkpoint owners for closed Yee, pencil PIC, open Hall, legacy hybrid,
  explicit EM, serial collisional and conductive/PML families. Explicit EM uses
  its fast interval. Numerical clock additions, legacy initial origins and
  scalar zero-step behavior are retained; exact-revision CI remains pending.

- Absorbing Yee history validation selects the same rational/error fields
  explicitly, preserving all predicates while avoiding heterogeneous dynamic
  field arithmetic. Its restart/rollback snapshot helper compares every capture
  field by value before detachment. Numerical maps, checkpoint schema and all
  assertions/JET gates remain unchanged. The actual minimum-Julia parent run
  retains 22 failures; exact-revision remote acceptance of this repair is pending.

- Implicit hybrid ready/checkpoint validation rejects `last_dt > time` even
  when its solve report matches the malformed history. Canonical nonnegative
  clock addition requires this bound; varying durations remain supported.
  The physical solver and IO publication are unchanged. Essential remote phase
  regressions are prepared; exact-revision CI remains required.


- Impact longitudinal validation now checks certificate buffers directly and
  projects only the relevant history/grid values. Checkpoint population wrappers
  retain their concrete storage, and required-key validation avoids an
  intermediate boolean tuple. All validation conditions and JET gates remain.
  The overflow fixture doubles its finite tracer charge because the finite-cloud
  wall gather is half the stored field; no finite production result is rejected
  to satisfy that fixture. Source review is complete; remote acceptance remains
  pending in the [audit tracker](docs/DEEP_DEBUG_TRACKER.md#current-progress).


- Electrostatic and electromagnetic collection preflights now require injective
  particle leaves,
  rejecting repeated-cell strided views before staging or checkpoint restore.
  Existing injective storage retains its admission. Source review only; runtime
  validation is deferred during the current implementation-first pass.

- The absent conduction-quotient validator now specializes on `HybridStepper`,
  resolving a dispatch ambiguity that blocked nonconducting evolving-pressure
  initialization, stepping and restart preflight. Real quotient-grid validation
  is unchanged; existing regression suites await GitHub revalidation.

- Original IEEE factor, targeted collection sign and timestep-zero admission use
  stored bits after existing conversion and finite checks. Native zero charge
  returns its signed zero before division; charged native quotient limits remain.
  Scalar Hybrid/CAM-CL checkpoints now reject nonfinite charge and nonpositive or
  nonfinite mass through the stepping factor owner, preserving finite signed
  charge and extreme original ratios. Other legacy snapshot value policies and
  schemas are unchanged. Essential actual-mode and atomicity checks await GitHub.

- `apply_antenna!` stages complete finite-stored-curl updates before changing B,
  using shared exact conversion for affine overflow classification and finite
  rescue. Nonzero steps require nearest/gradual mode and finite B/materialized
  curl; curl unrepresentability rejects even if dt could later rescue it.
  Ordinary native bits, converted signed-zero no-op and input sharing remain.
  Workspace/field, shared-B and writable-metadata conflicts reject before curl.
  Essential GitHub runtime/compiler/format and cost checks remain pending.

- Hybrid collection checkpoint essentials now follow CPU original charge/mass
  admission: finite factors may have an overflowing or zero-rounded quotient.
  Nonfinite mass still rejects with complete persistent-state and storage
  checks; production validation and checkpoint schema are unchanged.

- Electrostatic collection restore detaches complete incoming state before
  publication, preserving inputs that refer to live destination buffers. Stored
  fields restore verbatim; particle and owner capacity is reserved beforehand.
  The focused input-cycle case and existing ES suite await GitHub verification.

- Periodic ES/EM collections and both closed Yee PIC models share rejection of
  read-only `Base.Experimental.Const` and atomic-only particle storage through
  standard wrappers before publication. All leaves, fixed populations and
  checkpoint destinations are covered. Writable host arrays, views and ordinary
  reinterpretations retain their contract; exact GitHub runtime and compiler
  verification is pending.

- The closed Yee arithmetic-context control uses a pre-stored minimal subnormal
  and verifies its bits under each mode before checking structural zero. The
  production helper is unchanged; exact-revision GitHub verification is pending.

- Collisional model binding compares original scalar coefficient values at the
  same concrete model/closure type, allowing matching BigFloat-backed state to
  resume after copying or serialization. The bound model is detached once;
  changed original values still reject before publication.

- HybridStepper B-only and B/pressure RK4 stages share original-duration affine
  range completion. Coefficient, product and weighted-sum range events retain
  complete original h/rate operands before final conversion; ordinary admitted
  expressions retain their operation order. Physical-rate storage and existing
  candidate/rollback policies remain; [scope](docs/TRANSPORT_SOURCE_RANGE_PLAN.md#original-duration-rk4-affine-combinations).

- Shared timestep preflight rejects a positive full step whose magnetic
  substep `dt/NB` rounds to zero before evolution. True zero full timesteps and
  resolved subnormal field steps retain their existing admission.

- Metal gathered-field publication has a finite-sample status scan and guarded
  stores, with explicit native field storage admission even for empty batches.
  Disjoint updates call the shared native formulas; particle aliases retain the
  ordered parent broadcast schedule. Source review is complete; actual remote
  device, compiler and unchanged cost gates remain open in the
  [preflight plan](docs/METAL_GATHERED_VALIDATION_PLAN.md).

- Route EMPIC collection checkpoint payloads through their schema owner at the
  shared validation entry, following the existing electrostatic layout policy.
  Legacy electron/ion payloads retain their validator and cross-layout loads
  reject before publication.

- Bound normal IEEE relativistic Boris kick error at the scalar and both EMPIC
  batch owners, reusing the canonical momentum expression and residual arithmetic.
  Supported native tuples are retained when certified; bounded unresolved rows
  attempt certified DD repair. The [work order](docs/RELATIVISTIC_NORMAL_ACCURACY_PLAN.md)
  records the corrected witness/compatibility results and remaining extreme
  replay, priming, compiler and resource gates.

- Scalar relativistic kinetic-energy contributions at high speeds require a
  physical error enclosure, including reference residuals, before retaining
  historical arithmetic. The point bound is four adjacent spacings; low-speed
  operation order, nonfinite-factor propagation and complete range replay remain.
  The existing scalar sum retains its summation order. Runtime acceptance is
  pending GitHub alongside the collection diagnostic tests.

- Original-factor relativistic priming defines its typed coefficient method
  after the carrier type, removing the observed include-time undefined name.
  The coefficient test and numerical expressions are unchanged; the failed
  initial remote observation is retained in the audit evidence.

- IEEE bit publication binds each integer word width to a concrete scalar
  format, preserving the existing packing and rounding decisions. Completed
  numerical CI passes exposed an abstract-target compiler boundary; the new
  publication methods await full GitHub compiler and regression checks.

- Coupled ideal and Ohmic scalar pressure stages request the existing complete
  signed source owner when derivative or product range loss is detected. Exact
  advection/compression/heating cancellation remains admissible; the complete
  required rate retains explicit range checks. Ordinary stage order and public
  rollback remain. Essential regressions await GitHub verification; see the
  [pressure range plan](docs/SCALAR_PRESSURE_RATE_RANGE_PLAN.md).

- The carried S1 Ohmic field completes `E_force + ηJ` before range loss can
  erase finite cancellation or a subnormal rounding tie. Ordinary order, the
  `η=0` copy and independent force storage remain; essential both-type,
  in-place and separate-output controls await GitHub CI.

- CGL Ohm evaluation shares the scalar complete signed point-expression owner,
  applying its represented tensor force on all three components. Ordinary
  arithmetic/write order and tensor generation remain unchanged; public
  both-type D1–3 cancellation controls await GitHub CI.

- Scalar CPU Hybrid Faraday stages with `de2=ηH=0` retain original represented
  constitutive operands on finite product-range events, including Ohmic
  sources. The existing dyadic FFT owner preserves sub-storage products until
  their derivative is composed; physical E publication and ordinary native
  paths are retained. Normal hidden-variation admission, general FFT
  conditioning and other closures/operators remain open in the
  [source work order](docs/IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md).

- Scalar Ohm evaluation combines complete signed finite represented factors on
  intermediate overflow or subnormal products, preserving ordinary operation
  order and raw nonfinite propagation. This repairs motional/Hall cancellation
  and inverse-density-rescued Hall fields. Broader original-product derivative
  admission remains open; essential checks run on GitHub.

- CPU hybrid/CAM-CL and electrostatic collection validation retains finite
  original charge and positive mass even when their derived quotient overflows
  or underflows. Existing sign, storage, neutrality and arithmetic output checks
  remain enforced. Essential full-model and restart cases await GitHub CI.

- Native Metal uniform and gathered pushes validate original charge/mass and
  reject nonfinite, zero or subnormal derived quotients for nonzero charge
  before launch/mutation. True zero charge and normal ratios remain supported.
  Device factor recovery and priming remain open; see the
  [factor range plan](docs/PARTICLE_CHARGE_MASS_RANGE_PLAN.md).

- Relativistic EMPIC priming and kicks retain original species charge and mass
  through range repair and finite-output publication, including overflowing or
  subnormal derived quotients. Both model owners keep C1 rejection and whole-step
  rollback. Essential phase/restart regressions await GitHub CI; full normal
  relativistic accuracy remains open. See the
  [factor range plan](docs/PARTICLE_CHARGE_MASS_RANGE_PLAN.md).

- Particle kick and coordinate batch bypasses require known injective writable
  storage. Repeated-index views retain sequential guarded evaluation without
  narrowing the accepted particle containers. Uniform and gathered ordering
  remain distinct; focused exact-dyadic public controls run in GitHub CI.

- Euler and electrostatic priming/prediction use joint range and normwise
  accuracy admission, including a finite exact-result margin. Exact cold replay
  resolves nearest-overflow ties and finite neighbors that widened arithmetic
  can merge, retaining valid native bits under the unchanged four-ulp term-scale
  tolerance. Priming overflow errors and predictor infinities keep their
  existing policy. Overlapping velocity/field storage and repeated-index views
  retain sequential scalar guards when writes could invalidate a batch prepass.
  See the [CPU particle tracker](docs/BORIS_ACCURACY_IMPLEMENTATION.md).

- Complete represented transport sources before final range classification.
  Active hyperresistivity uses shared normalized Fourier channels for the
  ideal/Ohmic/hyperresistive field, signed pressure equation and source-aware
  Faraday rate. Required physical E, current Laplacian and rates retain explicit
  representability checks. See the [S2 contract and acceptance status](docs/TRANSPORT_SOURCE_RANGE_PLAN.md).

- Nonrelativistic leapfrog and electrostatic priming, EMPIC ion priming and the
  hybrid Euler predictor retain stored species charge and mass through range
  guards and exact cold replay. A quotient rounded to zero no longer certifies
  a zero physical kick. Normal native arithmetic, priming overflow errors,
  predictor infinities and leapfrog phase/history are retained. Verification
  status is recorded in the
  [CPU particle tracker](docs/BORIS_ACCURACY_IMPLEMENTATION.md).

- Nonrelativistic CPU Boris kicks retain original species charge and mass in
  public pushes and model kick owners. Joint range and accuracy certificates
  guard native batch bypasses; unresolved finite kicks use exact replay with
  the existing four-ulp tolerance and nearest-even overflow rule. Device factors
  and full normal relativistic accuracy remain separate follow-ons.
  See the [CPU Boris contract and tracker](docs/BORIS_ACCURACY_IMPLEMENTATION.md).

- Distributed FFT CI prepares one environment and runs core and metadata
  regressions in separate MPI processes. Launcher errors propagate as failures;
  both collections retain their existing tests and uploaded logs.

- Exact cold-path Float16/32/64 conversions assemble the rounded IEEE bits
  directly and extract stored inputs with integer arithmetic. Subnormal results
  and inputs no longer depend on the caller's flush-to-zero setting; rounding
  decisions and BigFloat precision/rounding semantics are retained.

- Finite particle drifts with an exactly zero displacement factor retain native
  IEEE bits without exact range work. The measured near-boundary scalar and
  complete public profiles now allocate zero Julia bytes; signed zeros and
  nonfinite propagation are unchanged.

- CPU particle drifts and gathered midpoints classify nearest-even overflow from
  the complete stored expression `x + dt*v`. This repairs both finite native
  results whose exact result overflows and nonfinite native results whose exact
  result is finite. Overlapping particle storage uses per-particle guards when a
  batch bound can become stale. Ordinary finite bits and nonfinite-input
  propagation are retained; unresolved finite inputs use exact arithmetic with
  documented cold-path costs. See the [coordinate contract](docs/PARTICLE_KINEMATICS_RANGE_PLAN.md#coordinate-overflow-follow-on).

- Relativistic Boris and leapfrog priming now reject rounded velocities outside
  the existing supported subluminal domain, including native batch bypasses.
  Public relativistic EMPIC advances restore particles, fields, gathers, FFT
  scratch and clocks on error using a lazy reusable undo workspace; initialization
  rejects unsupported required-species velocities before mutation. Generic-precision radiation
  stages three velocity vectors before checked publication. Previously returned
  out-of-domain velocities now raise `ArgumentError`; supported native output
  bits and the explicit checkpoint/archive schemas are retained. See the
  [publication contract and verification status](docs/RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md).

### Added

- `yee_longitudinal_voltage` evaluates the closed Yee path-work contraction
  using the same CIC overlap stencil as endpoint impulse and the existing
  stored-dx Maxwell metric. It rounds the complete exact expression once,
  preserves read-only inputs and rejects nonfinite rounded voltage. Force
  selection, coupled stepping and wall ledgers belong to separate coupled
  model contracts.

- Constant hyperresistive `TransportElectrons` (CAP-011 S2) for serial periodic
  CPU Float32/64 grids in one to three dimensions. Each stage pairs
  `ηJ − ηH∇²J` in Faraday with positive Ohmic and Nyquist-completed electron
  heat; ions keep the carried ideal force. The zero-`ηH` path and policy-1
  restart meaning are retained. Requires SpectralOperators 0.1.7. Current-main
  numerical, resource and compiler acceptance remains tracked in the linked S2 plan.

- Rank-owned periodic vector operators and filters: `gradient!`, `divergence!`,
  `curl!`, `laplacian!`, `project_divfree!`, `exp_filter!` and
  `dealias_two_thirds!` methods on `DistributedFourierGrid` reuse the
  SpectralOperators 0.1.3 multiplier, Laplacian, projection and filter kernels
  with no global gather. Vector fields are three-tuples of rank-owned fields;
  the serial signatures, alias rules, Nyquist handling and range scaling are
  kept. Validation, filter parameters and every local step are agreed
  collectively, and all outputs are published together after the last
  agreement. The grid now allocates two more spectral and two more real
  publication buffers at construction. Operand validation no longer
  allocates, which also removes 5.3 KB from each warmed `deriv!` call.

- `TransportElectrons(γ = 5/3)` (CAP-011 S1) adds constant Ohmic transport to
  the stored electron pressure of serial periodic `HybridStepper` runs. Ohm's
  law is split into the ion force `E_force` and `E = E_force + ηJ`: Faraday
  uses `E`, the first-step priming, Boris kicks and endpoint predictor use the
  carried `E_force`, and each RK4 stage heats the pressure by `(γ−1)η|J|²` from
  its own current, so the periodic continuum energy total closes (fourth order
  with fixed ions, second order coupled). `η = 0` reproduces `EvolvingElectrons`
  bitwise. In this S1 slice, density below `nfloor`, unsupported `ηH`/`de2`, scalar `q ≠ 1`,
  non-finite heat and invalid force storage are rejected with full rollback.
  Checkpoints and archives store `Eforce` with `transport_policy_version = 1`
  (archive schema unchanged). MPI, CAM-CL and Hall-MHD reject the closure.

- Bounded Maxwell vacuum/cavity fields (CAP-002A): `YeeGrid`, `YeeFields`,
  `MaxwellVacuum` and `step_maxwell!` evolve staggered Yee E/B in 1D–3D with
  periodic or homogeneous perfectly conducting axes, synchronized
  B-half/E-full/B-half steps, a downward-rounded CFL bound and transactional
  rejection that leaves physical fields and the clock unchanged. Stages are
  native stencils with per-axis multipliers computed once per stage (about
  6 / 11 ns per cell and step for Float32 / Float64 on 48³, zero allocations);
  out-of-range multipliers use Float64 or exact power-of-two scaling and
  overflowing entries are completed exactly, under the normwise accuracy
  contract. `yee_curl!`, `yee_divergence!`, `apply_yee_pec!`,
  `em_field_energy`, `yee_poynting` and `yee_poynting_flux` complete the API;
  see `examples/maxwell_cavity.jl`.

- Prescribed-current 1D Yee fields: `step_maxwell!(model, dt, impulse)` advances a
  1D `MaxwellVacuum` with a caller-supplied integrated current on the electric
  lattices (validated before staging; a zero impulse reproduces the vacuum step
  bitwise). `yee_electric_charge!` reports all-node charge including PEC surface
  charge, and `yee_source_work`/`yee_modified_energy` give the fixed-step work
  balance `Hnew - Hold = -work`; see `examples/sourced_maxwell_cavity.jl` and
  `docs/SOURCED_YEE_FIELDS.md`.

- Closed PEC particle charge/current primitives: `YeeChargeWorkspace`,
  `yee_particle_charge!` and `yee_longitudinal_impulse!` provide extensive CIC
  nodal charge and compatible longitudinal path impulse in CPU 1D Float32/64.
  All populations are validated before scratch writes and output is published
  only after complete native/exact evaluation. Broad collections and strided
  storage remain supported. Coupled force, reflection, transverse current,
  wall ledgers and restart use separate source models; GPU/MPI integration
  remains open. See the [coupled contract](docs/BOUNDED_PARTICLE_COUPLING_PLAN.md).

- `CollisionalHybridStepper(inner, species; dt, rng, coulomb_log, epsilon0,
  weight_to_number, u_floor)` composes `HybridStepper` with cell-local
  `collide_coulomb!` at a fixed step: after the push and magnetic advance, the
  half-level velocities scatter for the full step at the new positions, and the
  carried E is then formed from the scattered velocities (first order; bitwise
  equal to step, collision and E recomputation). A failed step restores
  particles, all public fields, clock, source ledger and the owned RNG.
  `save_checkpoint`, `load_checkpoint!` and `archive_run` persist the advanced
  RNG, so restarted runs continue bitwise in the same environment. Zero `dt` is
  a no-op, other changed steps are rejected, and zero logarithms reproduce
  `HybridStepper` exactly. CPU, D=1–3, Float32/64, NGP/CIC/TSC, ideal massless
  isothermal electrons.

- Rank-owned periodic Fourier derivatives: `distributed_fourier_grid(n, L;
  comm)`, `distributed_field(grid)` and `deriv!`/`deriv` methods on the
  distributed grid reuse SpectralOperators 0.1.3 kernels on PencilFFTs pencils
  with no global gather. Configuration, axis and local validation are agreed
  collectively, and output is published only after every rank succeeds.

- Cell-local multispecies Coulomb collisions: `collide_coulomb!(species, grid,
  dt; coulomb_log, epsilon0, weight_to_number, u_floor, rng, workspace)` pairs
  self and cross-species particles within each periodic cell using physical
  charges, masses, macro-weights and cell volume, with the shared Takizuka–Abe
  angle law and Higginson–Holod–Link weight rejection. `CoulombCollisionWorkspace`
  makes repeated calls allocation-free; `coulomb_normalization(units)` returns
  the permittivity and weight-to-count factor for normalized hybrid units.
  Velocities publish only after the whole call succeeds. The homogeneous
  `collide_coulomb!(ps, gcoeff, dt)` operator is unchanged.

- Distributed FFT plans accept `real_type=Float32` or `Float64` (the default).
  Precision participates in collective configuration agreement before any plan
  is created; unsupported values and rank mismatches fail explicitly.

- Oblique Hall-MHD production validation (case 35) evolves oblique modes from
  the unchanged Hall-MHD solver across angles, temperatures, branch signs and
  grid/time refinements, fitting frequency and damping from production
  histories against the independent dispersion relation. Its history fitter
  has permanent tests; the 60-profile campaign is a manual validation run.

- Production cold parallel-wave validation evolves both hybrid and Hall-MHD
  branches across 40 profiles and six timestep-order histories. Signed frequency
  and growth/damping fits retain raw data, residuals, conditional regression
  errors, window sensitivity and hashed source/environment receipts. The
  validation fitter handles long phase histories with compensated accumulation.

- Pinned NHDS and Hybrid-VPIC reference builds share a provenance owner with
  tool/input/output checksums, verified cache reuse, retained failure logs and
  atomic receipt publication. External validation reports missing tools and
  failed builds explicitly; builders no longer install global tools or delete
  legacy caches.

- Periodic CPU `ElectrostaticPIC` accepts arbitrary kinetic species in 1D–3D
  with explicit neutrality/background policy, one shared Poisson field and
  synchronized leapfrog phase. Collection updates stage particles and fields
  before publication; diagnostics count the shared field once. Checkpoints and
  archives preserve species order, policy and phase, including writable views.
  Scalar electron constructors retain their existing API and checkpoint schema.

- A source-backed completeness plan and live implementation tracker separate
  existing solvers, partial building blocks and missing coupled capabilities.
  Each stage records owners, dependencies and correctness, restart, hardware,
  performance and CI gates; planned capabilities are not advertised as complete.

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

- Require registered SpectralOperators 0.1.8 for certified directional axis-null
  conditioning. The dependency passes all five CI jobs; signed mixed-axis
  pressure and complete-Faraday consumer controls await HPIC GitHub checks.
  General active-mode conditioning remains outside this certificate.

- Moment deposition computes each particle's stencil once for density and all
  momentum or current components, wraps periodic indices without integer `mod`,
  and replaces per-stencil-point range checks with a per-particle normwise gate
  plus a final finiteness check. Ordinary outputs are bitwise identical; a 2D
  hybrid step with 131k particles runs about 2.3× faster. Cases whose smallest
  stencil product would underflow, or whose accumulated cells overflow, still
  take the exact replay. The documented `work` scratch argument of `momentum!`
  and `current!` is no longer written.

- Linear whistler packets with Float16/32/64 timesteps and fixed-width step
  counts reuse one CN or Euler multiplier per mode, preserving each multiplication
  and avoiding a multiplier cache allocation. Positive BigInt counts that fit
  Int now use that same iterator, removing repeated arbitrary-precision counter
  allocation. Generic real timesteps and other counters retain the original
  step-major evaluation; zero steps perform no multiplier calls. Existing
  allocation gates remain unchanged, with a new bounded-BigInt cost regression.

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

- `push_gathered!` rejects gathered fields overlapping particle positions or
  velocities before publication on CPU and Metal. Read-only E/B and weight
  sharing remain supported; reversed field views can no longer observe
  velocities changed earlier in the same CPU push.

- BigFloat-backed `recommended_dt` calls enclose the complete logarithmic
  recommendation with directed bounds and round a positive lower result at the
  caller's precision. Low ambient precision no longer enlarges the bound;
  stored input precision and task settings are retained. BigInt subcycles and
  Rational{BigInt} parameters use the same core. Uncertifiable MPFR exponent
  contexts raise ArgumentError; native arithmetic and CN's unused NB remain
  unchanged. The separate concurrent MPFR shared-cache runtime gap remains open;
  directed enclosure assumes correctly rounded primitives.

- Relativistic Boris and leapfrog priming resolve the Lorentz gap from stored
  velocity components near the light cone. Priming uses residual intervals to
  certify a rounded result before high-precision replay. Well-conditioned
  certified inputs retain historical bits; other priming results can change.
  The guarantee is local normwise accuracy, with near-light trajectory and
  representable output-cone limits described in `docs/PARTICLE_PRIMING_INTERVALS.md`.

- `push_gathered!` rejects midpoint outputs that overlap particle arrays,
  gathered fields or other midpoint components before mutation on CPU and
  Metal. Disjoint buffers and shared read-only E/B inputs retain their behavior.

- Distributed Fourier operators reject operands overlapping read-only
  wavenumber, derivative, integer-mode or Nyquist metadata through the shared
  collective admission check, before changing public or private arrays.

- `recommended_dt` bounds combined Hall-wave and resistive/hyperresistive
  RK4 modes with a certified inner stability polygon. Separate axis limits
  previously admitted a production Fourier mode that grew by 46.7% per step.
  Native undamped arithmetic and native non-RK4 recommendations are unchanged.

- Require SpectralOperators 0.1.5 for complete inverse-FFT headroom. A finite
  cutoff-kernel convolution can grow inside the unnormalized inverse even when
  every final output is representable; the dependency bounds that complete
  operation before transforming. Independent shifted and anisotropic MPI
  regressions cover Float32 and Float64.

- Infinite-strength exponential filtering preserves the DC mean and removes
  every other mode even when a large accepted power underflows. The shared
  repair is required through registered SpectralOperators 0.1.4, with real MPI
  consumer regressions for Float32/64 and both infinite and overflowing strength.

- Radiation drag retains finite damping when a small force beside a cancelling
  large field is amplified by the coefficient/timestep product. Scaled
  evaluation includes that amplification in its error certificate; unresolved
  cases use exact original operands and task-local integer arithmetic. Ordinary
  calls retain their arithmetic and zero allocation; see
  `docs/RADIATION_FORCE_SPAN_PLAN.md` for the error derivation and cold costs.

- `resistive_dissipation` preserves finite Float32/Float64 Ohmic integrals when
  current squares, coefficient products or the cell volume overflow or underflow
  before the complete result is formed. Zero resistivity with large finite
  currents returns zero; ordinary inputs keep their summation order.

- Checkpoint configuration ownership no longer rejects a BitArray view disjoint
  from state-owned bit views on Julia 1.13.3, where `Base.dataids` reports every
  view as sharing the packed chunks; overlapping views, the whole array and raw
  chunk references remain rejected on every Julia version.

- Distributed FFT roundtrip diagnostics now support in-place complex, real and
  real-to-real transforms. They compare the container's input view and reuse the
  input as the default output, as required by PencilFFTs. An explicit separate
  output still rejects. Real MPI regressions verify the global maximum error
  and restored input for both transform storage modes.

- Distributed FFT construction validates sizes, transform configuration and
  periodicity collectively before creating topology. Mismatched ranks and empty
  pencil partitions fail explicitly; a recoverable local planning error reaches
  every rank. Mixed integer types in a valid size tuple no longer recurse into
  a stack overflow. Serial-reference and injected-failure tests run on real
  2-, 4-, and 8-rank MPI in the existing CI jobs.

- Finite Gaussian Coulomb draws no longer reject when rounded variance just
  above one overflows the tangent after the exact regime test selects Gaussian
  scattering. The shared sampler retains the tangent's binary scale and
  completes the physical velocity update. Nonfinite draws still reject before
  publication, with the same RNG order and rollback contract.

- Prescribed-rate Coulomb collisions accumulate normalized population weights
  with compensation in Float64 for hardware-float particles and retain that
  mean through variance construction. BigFloat accumulates with local
  count-dependent guard bits and rounds the mean at its active precision.
  Large or sparse Float16 and coarse-BigFloat populations no longer
  select the wrong scattering regime from a stalled sum or overflowing count.
  Zero weights remain included in the population count; validation, RNG ordering
  and transactional velocity publication are preserved.

- Mode-history least squares centers observations before reducing their mean
  and evaluates residuals in centered coordinates. Constant histories now
  have zero rate and residual; large amplitude or phase offsets no longer
  enter residual reconstruction. Existing validation tolerances are unchanged.

- Exactly sampled identity rotations in prescribed-rate Coulomb collisions
  preserve original velocities, including tiny components and signed zeros.
  Angle validation and unequal-weight acceptance draws retain their order.

- Prescribed-rate Coulomb scattering retains original rate factors through
  intermediate overflow or underflow and applies physical velocity increments
  directly. Small nonzero kicks, subnormal pair-weight ratios and cancellation
  near a velocity exchange are checked against independent complete-expression
  references. Historical floating-point trajectories change; same-version
  seeded determinism and RNG event order remain covered.

- Validation frequency extraction includes the first Fourier bin and both
  rotation directions, removes DC before peak interpolation, and rejects
  histories without usable oscillation evidence. Component centering and
  amplitude normalization preserve small oscillations beside large DC values.
  Invalid sampling intervals raise `ArgumentError`; case 20 retains its
  original physical acceptance thresholds.

- Linear Crank–Nicolson multipliers retain finite phase components across
  overflowing/underflowing products, integer inputs and near-quadrature
  cancellation. Ordinary BigFloat evaluation forms exact private intermediates;
  extreme exponents remain separate from significands. Both preserve caller
  rounding settings. Exact rational results
  widen to `Rational{BigInt}` when the former integer type is insufficient.

- Whistler reference frequencies combine the guide field, density and squared
  wavenumber before exceptional rounding. Complete frequencies outside finite
  Float64 range raise `ArgumentError` before packet construction. Narrow hardware
  coefficients use Float64 division, and unsigned grid counts retain negative
  band-limited modes.

- Electrostatic field energy retains finite complete-expression results when
  intermediate squares or cell-volume products overflow or underflow.

- Unit conversion evaluates complete physical expressions when intermediate
  scales or final operations lose floating-point range. The shared exact-root
  owner separates significands from exponents, preserving signed endpoints and
  working BigFloat precision/rounding without changing numerical settings.
  All six stored `PlasmaUnits` fields now require finite positive values,
  including positional construction. Ordinary arithmetic order is retained.
  CUDA/Metal adapters keep device layout and backend through exceptional host
  staging; ordinary device inspection synchronizes. Native Metal boundary,
  storage, complex-value and flag-lifetime regressions pass on Julia 1.12.7 and
  1.13.2; native CUDA execution remains unverified. See the
  [normalization plan](docs/NORMALIZATION_RANGE_PLAN.md) for costs and open gates.

- Shared exact BigFloat fallbacks round the complete finite rational even when
  its integer factors exceed the active MPFR exponent range. Guard/sticky bits
  preserve the caller's precision and rounding mode through compact overflow
  and underflow endpoints; production code does not change numerical settings.
  Public AMR cell-center and gyroaverage regressions cover retained-precision
  inputs. Independent rounding-cell, concurrency, resource and fixed-cost
  checks pass on Julia 1.12.7 and 1.13.2.

- Every `BackendMemoryStatus` construction form now rejects negative,
  noninteger and out-of-range byte counts through the same validation owner.
  Keyword and full-positional construction previously bypassed the checks.
  Valid integer/missing counters, metadata conversions and independent telemetry
  samples retain their behavior; pressure clamping is unchanged.

- `PerpShock3D` Crank–Nicolson field steps now solve the frozen-moment induction
  equation with restarted GMRES and a finite true-residual check. This replaces
  four unchecked Picard iterations that could amplify resistive decay. The
  complete step commits only after all subcycles and final fields succeed;
  failure preserves caller state. `ShockCNWorkspace` reuses trial/solver
  storage, and explicit tolerances/iteration controls expose solve acceptance.
  CN requires CPU Float32/Float64 and ordinary particle `Vector`s; RK4 remains
  the default. The reduced induction model and checkpoint state are unchanged.

- Coulomb velocity updates follow the normwise accuracy standard: widened
  arithmetic replays only range loss (nonfinite intermediates, subnormal
  tangent/variance, mass fractions or pair scale) instead of every
  near-cancelling component, which stays within a few ulps of the pair scale.
  The cell-local operator hoists its rate coefficient once per cell (self and
  triplet) and once per distinct reuse value (cross). Per-pair throughput
  roughly doubles in both precisions; seeded runs differ from before only in
  cancellation cases, by a few ulps of the pair scale.

- External NHDS and Hybrid-VPIC validation READMEs describe the actual failure
  behavior, comparison scope and acceptance limits. Moving reference revisions
  and unchecked cached output are documented; historical numbers are not
  presented as current verified results.

- NGP nearest-node selection preserves values immediately below a half-cell
  boundary and uses upward ties consistently on CPU and Metal. CIC and TSC
  retain their existing rounding policies.

- Metal uniform-field particle pushing uses a device kernel calling the shared
  Boris update, removing recursive dispatch. Invalid charge, mass, timestep
  and uniform field values reject before particle mutation, including empty sets.

- Legacy electromagnetic PIC retains complete longitudinal current coefficients
  when intermediate modal division loses floating-point range.

- Checksum subprocess tests keep the resolved active package environment and
  verify source identity when run from an external test environment.

- Float32 field-energy reduction checks accumulated rounding error and retains
  complete finite integrals when the ordinary sum cannot certify accuracy.

- Ordinary interval validation creates indexed error labels only on rejection,
  preserving outputs and random-number state while removing measured allocations.

- Particle rank lookup and migration use the same half-open periodic seam as
  position wrapping. Migration prepares replacement storage before publishing
  it and coordinates storage/count preparation failures across MPI ranks.

- Ionization retains threshold, cooling and optical-depth decisions when
  intermediate floating-point arithmetic loses range or endpoint accuracy.

- The electron-inertia energy reservoir retains complete finite current squares,
  its coefficient and cell spacings through range-sensitive integration.
  Ordinary results and allocations are preserved. Zero magnetic/current states
  avoid exact-arithmetic traversal when cell-volume or coefficient products
  overflow, while nonzero underflow and nonfinite inputs retain their handling.

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

- Radiation reaction no longer returns an unchanged or zero velocity where the
  frozen-rate substep is finite but an intermediate overflows or underflows at
  extreme light speeds, fields or coefficient–timestep products (RR-002), and
  weak positive drag no longer increases a stored speed through the momentum
  round trip (RR-003). Float16/32/64 particles are scaled by
  `1/√(1 + (1 − |v|²/c²)(e^{2α} − 1)) ∈ [0, 1]`, evaluated natively with a
  compensated `1 − |v|²/c²` and an a-posteriori rounding bound; flagged
  particles are re-evaluated in Float64 or exponent-scaled double-double
  arithmetic, without BigFloat and without allocation. Results are normwise
  accurate to a few ulps of `|v|`, so ordinary outputs differ from the previous
  arithmetic at the rounding level. BigFloat particles keep the ambient-precision
  momentum update.

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
