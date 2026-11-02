# Plasma simulator completeness and implementation plan

Initial source assessment: **2026-10-23**, revision
`bf2a8bc615f514f8d2868934849bf44814140bbb`. This document plans work; it does not
certify the proposed implementations. Source status was reconciled on
**2026-11-02** against the IEEE, CPU particle-factor and transport
implementations, accepted C1 `1bc618b` and the bounded reflection/capture source candidates;
earlier line references remain pinned to the initial assessment. The
[audit tracker](DEEP_DEBUG_TRACKER.md) owns accepted fixes, measured results and exact-revision CI status.

**The package is not yet a complete implementation of the expanded simulator
goal.** It contains implemented solver families and an extensive test suite,
but several advertised domains still lack an end-to-end evolution path.
Source presence, a passing primitive test and a complete simulation are separate
milestones. No current evidence establishes that the package is globally bug-free,
fully optimized, or stronger than every external plasma code.

**Recorded-milestone snapshot (2026-11-02):** 14 checked and 65 open
checkbox items out of 79 (17.7% checked). This is an unweighted count of
this plan's recorded milestones. Checked source milestones can still await
GitHub acceptance; the count is neither weighted physics/runtime completion
nor a measured percentage of the full goal. All eleven CAP areas retain open
work.

**Current validation policy (2026-11-02):** combine source/caller review with
exact-revision GitHub CI for every push. Package, JET/compiler, formatting and
MPI checks run remotely; local Julia/device workloads remain excluded. Earlier
source-only publications retain pending acceptance until their corresponding
checks pass. Failed checks must be diagnosed and repaired without weakening
assertions or acceptance limits. The [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress)
records the active runs and failures.

The canonical [complete energy-budget source](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#current-complete-budget-source-2026-11-02)
now composes finite built-in original reservoirs when normal-range/sign bounds
do not certify the native sum. Scalar component contracts and true total Inf,
custom/nonfinite and CGL singular behavior remain. No acceptance milestone is
advanced: numerical/compiler/concurrency/cost gates and historical V6 rejection
remain open pending exact-revision remote and resource evidence.

The original scope recoverable from the repository is the periodic hybrid model,
its companion solvers and the open completeness roadmap: [README](../README.md),
[roadmap](DEEP_DEBUG_TRACKER.md#completeness-roadmap),
[distributed FFT decision](distributed_fft_evaluation.md),
[multispecies design](MULTISPECIES_HYBRID_PLAN.md),
[electron energy contract](HYBRID_ENERGY_MODEL_PLAN.md) and the two
[ray-tracing](superpowers/specs/2026-07-05-ray-tracing-design.md) /
[Raycon](superpowers/specs/2026-07-05-raycon-port-notes.md) specifications.
The current module docstring links directly to this plan as the implementation
checklist. The repository contracts above provide the basis for this expanded
work; this plan does not reconstruct an unavailable original specification.

**Scope distinction:** CAP-001–007 expand gaps already present in the repository's
roadmap into implementable work. CAP-008–011 are proposed model-scope expansions
identified by this assessment. Their absence is not a defect in the explicitly
limited existing models, and this document does not silently promise every
possible plasma regime. These expansions are included in the active goal; each requires a concrete
public model contract and the staged evidence below before implementation acceptance. Collision/field
composition under CAP-010 is also a dependency if the new local collision operator
is to be offered as a restartable collisional simulation. No additional user
approval is presumed necessary for work already authorized by the active goal.

## Assessment boundaries and evidence

This pass reads source, call sites, public contracts, tests and primary method
references. **The initial assessment ran no Julia, GPU, MPI or timing workload**,
to preserve another track's measurement window. “Verified” in its source matrix means
the stated source/contract fact was inspected, not that all numerical behavior
was rerun. Subsequent runtime findings are recorded in the live audit tracker
and stage status. Unmeasured risks remain hypotheses with explicit settling tests.

Line references below are anchored to the revision above unless a later accepted
change is explicitly identified. Recheck them when a dependent slice changes the owner. Search methods, exports and callers as well
as filenames before changing a capability's status.

Source names below resolve to these canonical files:
[hybrid stepper](../src/integrators/extrapolated_leapfrog.jl),
[hybrid fields](../src/fields/field_state.jl),
[ES](../src/models/electrostatic.jl), [EM](../src/models/full_pic.jl),
[Hall-MHD](../src/models/hall_mhd.jl),
[MPI](../src/parallel/mpi.jl),
[Pencil extension](../ext/HybridPlasmaPICPencilFFTSExt.jl),
[CUDA extension](../ext/HybridPlasmaPICCUDAExt.jl),
[Metal extension](../ext/HybridPlasmaPICMetalExt.jl),
[linear CN](../src/integrators/semi_implicit.jl),
[shock3d](../src/initial_conditions/reflecting_wall_shock3d.jl),
[shock2d](../src/initial_conditions/reflecting_wall_shock2d.jl),
[curvilinear mesh](../src/meshes/curvilinear_mesh.jl),
[AMR](../src/meshes/adaptive_mesh_refinement.jl),
[gyrokinetics](../src/models/gyrokinetics.jl),
[instability drivers](../src/verification/instability_sweep.jl),
[checkpoint](../src/io/checkpoint.jl),
[case 27](../validation/27_nhds_dispersion_comparison/run.jl),
[case 28](../validation/28_hybridvpic_perp_shock/run.jl) and
[case 32](../validation/32_magnetic_reconnection/run.jl).

| Capability | Existing implementation and checked evidence | Completion state |
| --- | --- | --- |
| Periodic hybrid, positive ion collections | `HybridStepper` / `CAMCLStepper`; README:78–93,156–169; `extrapolated_leapfrog.jl:67–120,359,766`; collection and pressure tests | Implemented within documented closure, geometry and time-centering restrictions; inherited audit findings remain open |
| ES/EM PIC | Periodic `ElectrostaticPIC`, `EMPIC` and `EMPIC1D`; ES collection implementation `d46029e` | Periodic CPU ES collections are landed; CPU EMPIC collections with opt-in D1–3 dividing schedules and specialized EMPIC1D S1/S2, both with optional relativistic motion are source-implemented, exact-revision GitHub acceptance pending. Broader schedules, bounded coupling and resident GPU/MPI integration remain open |
| Hall-MHD | `HallMHDModel`, `HallMHDState`, `step_hall_mhd!`; `test_hall_mhd.jl`, `test_io_hall_mhd.jl` | A periodic solver exists; it is not a general geometry/closure solver |
| Coulomb collisions | Canonical homogeneous and local collection `collide_coulomb!`; `CollisionalHybridStepper` | Local self/cross CPU operator and fixed-step isothermal hybrid/RNG restart exist; ideal evolving-scalar, S1/S2 transport and S3b conduction extensions await CI, bounded generic and actual-pencil collisional ES/EMPIC sources are implemented with deferred acceptance; broader compositions and parallel/device support remain |
| Nonperiodic EM | Periodic `EMPIC`; separate Yee vacuum/prescribed-current fields; closed 1D PEC charge/impulse/path voltage plus optional finite physical-width cloud static coupling | B2a, optional cloud, event-free longitudinal/1D3V EM, one-cluster reflection and one isolated passive longitudinal capture/restart source are implemented. Exact-head GitHub numerical/JET/CI and resource acceptance, general capture laws, multiple clusters/captures, transverse walls, measured absorber behavior/stability, general CPML and backends remain CAP-002 work; separate scalar-conductivity and causal directional ADE PML field layers with full budget/restart are source-implemented |
| Distributed evolution | Replicated-field MPI reference, rank-owned Fourier operators, static Poisson/vacuum Maxwell, coupled ES/EM and isothermal hybrid physical-pencil source | Ownership, complete CIC transport, direct migration, canonical movers/field expressions and scalar diagnostics are source-implemented; Maxwell, ES, EM and complete eighteen-field hybrid matching-layout restart are source-implemented; exact-head runtime/restart/resource acceptance, broader closures/EM families and cluster scaling remain CAP-003 work |
| GPU | Metal extension has deposit/gather/push methods; CUDA extension:1–68 has storage/telemetry hooks; `HybridFields` uses `Array` at field_state:12–29 | Partial kernels, not complete resident solvers, CAP-004; device availability is not numerical acceptance |
| Implicit evolution | Linear Fourier reference and residual-checked reduced shock CN; separate [synchronous coupled hybrid source](IMPLICIT_HYBRID.md) | CPU periodic isothermal collections now have joint velocity/B residual, bounded Newton–GMRES, work-defect report and complete persistence source. Runtime/method acceptance, evolving pressure/orbit averages, EM and backends remain CAP-005 work |
| Geometry | `ToroidalGrid`, metric gradient/divergence in curvilinear_mesh:20–49,90–235 | Geometry operators exist; no coupled toroidal plasma stepper, CAP-006A |
| Gyrokinetics | Local `GuidingCentre` utilities plus [self-consistent slab source](SLAB_GYROKINETIC.md) | Full-f gyrocentre collections, reciprocal ring-TSC, full-k polarization/modified adiabatic electrons, collective RK4, physical diagnostics and complete synchronous checkpoint/archive are source-implemented and reviewed. Exact-revision runtime/scientific/resource acceptance and broader geometry/physics remain CAP-006B work |
| AMR | `AMRGrid` is a single 1D scalar level, adaptive_mesh_refinement:8–23; prolong/restrict/refine:322–366 | Transfer primitives; no adaptive plasma hierarchy, particle crossing or coupled timestep, CAP-006C |
| External validation | Pinned external references plus case 34 parallel and case 35 oblique Hall production histories | Bounded production validation exists; complete manual campaigns, broader regimes, full profile replay and matched cost/accuracy comparisons remain CAP-007 work |
| Kinetic species diversity | Arbitrary ES collections landed; CPU EMPIC collections with opt-in D1–3 dividing schedules and specialized EMPIC1D S1/S2, both with optional relativistic motion are source-implemented and reviewed | Exact-revision GitHub and resource acceptance, broader schedules and collision/boundary/partition/device composition remain CAP-008 work |
| Multidimensional shocks | Existing reduced shocks plus separate [OpenHallHybrid source](OPEN_HALL_HYBRID.md), built on real mixed operators | Bounded isothermal numerical reservoirs, full-Hall projected field trials, SBP CIC, synchronous particles/sources, energy and full RNG restart are source-implemented/reviewed. Runtime/physical-boundary/resource acceptance, evolving pressure and broader boundaries/backends remain CAP-009 work |
| Coupled sources and collisions | `CollisionalElectrostaticPIC` / `CollisionalEMPICStepper` source with deferred acceptance; `CollisionalHybridStepper` and full RNG checkpoint/archive state (`2b2420b`); ideal evolving, S1/S2 transport and S3b conduction pressure/force extensions are CI-pending | Fixed-step first-order CPU composition plus bounded actual-pencil ES/EM physical/RNG source families; runtime acceptance, anisotropic/inertial models, electron-ion heat exchange, other sources, schedules and backends remain CAP-010 work |
| Electron transport | `EvolvingElectrons` covers ideal massless pressure; `TransportElectrons` provides constant Ohmic transport (S1), with complete-source hyperresistive S2 and shared scalar pressure range routing. `ConductingElectrons` adds an original-quotient S3b source candidate with bounded preparation and projected full-Nyquist conduction | [S2 gates](TRANSPORT_SOURCE_RANGE_PLAN.md), [scalar pressure verification](SCALAR_PRESSURE_RATE_RANGE_PLAN.md) and [S3b GitHub verification](EVOLVING_ELECTRONS.md#constant-isotropic-conduction-source-candidate) remain open; broader conditioning, complete-step cost, conduction models, equilibration, anisotropy, inertia and backends remain CAP-011 work |

The geometry integration finding also follows from a repository-wide symbol
search: outside their defining files, `ToroidalGrid`, `AMRGrid`, `GuidingCentre`,
`push_guiding_centre!` and `gyroaverage` occur in `src/` only as exports.
`test_fusion.jl` checks the primitives; those tests do not construct a coupled
fusion simulation.

## Completion contract shared by every slice

A slice becomes complete only after all applicable boxes below have evidence.
Unsupported combinations must fail at construction/preflight, before mutation;
they are recorded in the compatibility matrix rather than silently approximated.

- [ ] State equations, units, normalization, spatial/temporal centering, admissible
  inputs and model limits. Specify which conservation identities are discrete,
  convergent, statistical or inapplicable because a boundary/source exchanges them.
- [ ] Name canonical owners and every affected caller. Preserve documented
  constructors and checkpoint readers or provide an explicit schema/API migration.
  Share physics kernels where the equations agree; distinct discretizations need
  distinct state, not duplicated physical policy or a universal speculative framework.
- [ ] Freeze baseline source, resolved dependencies, machine configuration and
  inputs before editing. Retain failures and rejected candidates. Complete the
  pre-edit brief and at least three wrong-patch/mutant challenges.
- [ ] Test exact/manufactured solutions, independent discrete oracles and held-out
  invariants; vary dimensions, precision, signs, empty populations, shape order,
  species, grid parity, crossings and timestep histories. Use timestep/grid/particle
  refinement separately so noise does not masquerade as truncation error.
- [ ] Exercise nonfinite/out-of-range inputs, aliases, invalid configuration,
  resource limits, incomplete solves and I/O failure. Define the commit point;
  before it, rejection restores public physical state, clocks and persistent phase.
  State precisely whether RNG/scratch may advance. Distributed transport failure
  is distinct from recoverable collective validation failure.
- [ ] Verify uninterrupted versus restarted evolution, corrupted/versioned input,
  ownership, object cleanup and same-environment reproducibility. Cross-backend or
  changed-rank continuation needs numerical/statistical equivalence unless a
  deterministic reduction/order contract has actually been implemented.
- [ ] Run minimum-supported and latest-stable Julia in exact-source environments;
  record actual versions. Run focused tests, neighbors, Aqua, concrete JET and full
  package JET trees. Classify every changed report; a lower count is not a diagnosis.
- [ ] Run real devices/ranks when claimed. A skip, emulation, `COMM_SELF` run or
  CPU conversion is not a hardware/cluster pass. Required unavailable checks stay open.
- [ ] Complete fixed performance gates below, then independently audit the whole
  dataflow from a different entry point, including restart/error paths and tests'
  ability to reject wrong implementations.
- [ ] Update examples, model support table, diagnostics and restart schema. Commit
  and push each accepted implementation/fix separately to `main`; require all six
  exact-SHA GitHub jobs green, or retain explicit evidence of a billing block.
  Preserve hardware/large validation results separately from ordinary GitHub CI.

### Performance acceptance, fixed before candidate measurement

“Fully optimized” is not a measurable stopping rule. Each implementation must
instead identify workloads, dominant costs and a concrete target, and retain an
open profile ledger for remaining bottlenecks.

1. Freeze a workload matrix covering empty/small/production-sized populations,
   all supported dimensions/precisions, shapes, closures, species ratios and
   ordinary/exceptional numerical paths. Include complete steps and diagnostics,
   not only isolated kernels. Record setup/compilation separately from steady work.
2. For behavior-preserving refactors, require unchanged complete outputs and
   unchanged warmed allocation where already allocation-free. A physics-method
   change instead uses identical physical error, conservation and uncertainty
   targets; different arithmetic order cannot be called bitwise equivalent.
3. Pre-register the baseline, repetitions and rejection rule. For unchanged paths,
   the existing audit's screening ceiling is **1.5** for both ratio of medians and
   median of paired ratios, not a claim that a 50% slowdown is an optimization.
   Record every profile; do not average away failures. An optimization needs a
   lower predeclared target at its measured hotspot and unchanged neighboring gates.
4. Warm, specialize and calibrate identical-work and deliberately doubled-work
   controls; use alternating paired samples in an exclusive window. For GPU work,
   synchronize the measured stream and report kernel plus full-step wall time,
   transfers, device/host allocations and peak live storage. For MPI, time the
   slowest rank and report communication volume and imbalance.
5. If no working parent exists, freeze a correct serial implementation or trusted
   external algorithm as the numerical baseline. Establish throughput and memory
   measurements; do not invent a speedup against a failing or nonterminating path.
6. Before a production/scaling claim, record time-to-solution at matched error,
   strong/weak scaling, peak memory and failure recovery on named hardware.
   Select quantitative budgets after the baseline campaign and before changing
   the algorithm. An unset budget leaves that acceptance gate open.

## Sequencing and integration ownership

The [parallel implementation workflow](IMPLEMENTATION_WORKFLOW.md) owns current
worker assignments, exclusive file paths, handoffs and the next ready work.
Independent source work continues while GitHub checks run; integration
acceptance still requires the relevant exact-revision results.

Finish inherited correctness/resource gates before using the repaired owner as a
new baseline. CAP-001 can proceed independently of field redesign. CAP-002,
CAP-003 and CAP-004 share storage/operator requirements and need a short joint
interface review, but they must not wait for a large abstract-framework rewrite.
Implement the smallest storage/operator change needed by the first real solver,
then extend it with the next backend while retaining the serial oracle.

| Order | Coherent work | Dependency / merge boundary |
| --- | --- | --- |
| 0 | Existing arithmetic, NGP, GPU push, budget and compiler findings | Exact-source local gates and CI; held radiation work remains separately owned |
| 1 | CAP-001 local collision kernel; CAP-007 reproducible references; CAP-005 existing shock solve audit | Independent tracks; no simultaneous timing acceptance |
| 2 | CAP-002A/B closed-boundary EM; CAP-003A/B local ownership/operators; CAP-004A/B one resident solver | Shared array/operator contracts reviewed once; separate physics changes and pushes |
| 3 | CAP-008 species collections; CAP-010 composed evolution and restart; CAP-011 dissipative energy | Depends on phase/particle growth/diagnostic contracts, not on unrelated hardware |
| 4 | CAP-002 open boundaries, CAP-003 coupled distributed step, CAP-004 additional solvers | Compose already accepted kernels; rerun intersection cases |
| 5 | CAP-005 coupled implicit hybrid/EM; CAP-009 full-Hall shocks | Reuse accepted operators, boundary ledger, species and nonlinear solve contracts |
| 6 | CAP-006A geometry, CAP-006B gyrokinetics, CAP-006C adaptive plasma hierarchy | Distinct physics models; staged analytic and independent-code validation |
| Every stage | CAP-007 validation and independent final audit | A capability closes only with integration, failure, restart and cost evidence |

### Implementation tracker

“Planned” means source assessment and an acceptance outline exist; it does not
mean code or runtime verification is complete. “Scope proposal” marks the
additional CAP-008–011 model expansion. Root owns integration, tracker consistency,
commit/push and CI. A named subsystem below is a code owner; an unassigned
implementation agent is explicitly unassigned, not implicitly running.
Current assignments (2026-11-02): static composite ES implements in its own
worktree; root integrates frozen implicit EM after periodic curvilinear source, while
a freed worker maps resident CUDA. Hall and slab candidates are source-complete
and await remote acceptance. Root owns shared registrations, current docs, dependency
publication, source integration and remote CI. Earlier audit assignments and
quota failures describe historical work; the linked workflow is authoritative
for active ownership.

| Slice | Implementation owner / canonical code | Dependencies | State and checked evidence | Next concrete action |
| --- | --- | --- | --- | --- |
| 001 local self/cross collisions | Canonical collision/sorting owners | Accepted sampler and particle arithmetic | CPU self/cross operator, normwise scheduling and fixed-step isothermal composition landed | MPI/device pairing and broader cost campaign; [live evidence](DEEP_DEBUG_TRACKER.md#current-progress) |
| 002A/B bounded Maxwell + particles | Canonical mesh/operator, field, boundary, coupling and diagnostic owners | Accepted vacuum/prescribed-current fields; full particle/wall/DC contract required | Closed 1D PEC nodal CIC charge, longitudinal endpoint impulse and its exact scalar path-work adjoint are implemented with explicit IEEE scalar boundaries and exact mode decisions. An optional triangular-cloud family adds static stationary work-adjoint fields with an explicit distinct wall shape. Event-free longitudinal/1D3V EM, one-cluster reflection and one isolated passive capture source candidates include diagnostics, staged publication and distinct restart | [Coupled contract](BOUNDED_PARTICLE_COUPLING_PLAN.md), [reflection](YEE_IMPACT_LONGITUDINAL.md) and [capture](YEE_ABSORBING_LONGITUDINAL.md): exact-head numerical/JET/CI and resource acceptance remain pending; broader wall laws, multiple clusters/captures, transverse walls and backends remain open |
| 002C/D open boundaries + production integration | Canonical Yee field, conductivity, boundary and I/O owners | 002A/B; 003/004 for parallel/device claims | Scalar-conductivity and causal directional ADE PML layers, exact signed budgets and distinct complete restart source implemented; no coupled open PIC path | Accept reflection/convergence/stability/restart/cost; derive particle/source charge-flux integration, general CPML and backend extensions |
| 003A/B local storage + distributed operators | Distributed audit track; MPI/Pencil/SpectralOperators owners | Collective precision, shared kernels and owned storage | Rank-owned storage and derivative/vector/filter operators landed; `4724947` exact-SHA CI/compiler gates pass | Shared filter/inverse repairs are registered as 0.1.4/0.1.5; consumer and metadata repairs are pushed with local integration gates passed. Their exact-SHA CI is tracked separately; next complete [resources and integration](DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md) |
| 003C/D distributed coupled solvers | MPI orchestration + canonical solver owners | 003A/B and accepted serial solver | Rank-owned periodic Maxwell, coupled CIC ES/EM and isothermal hybrid source with global scalar diagnostics implemented | Maxwell, ES, EM and hybrid matching-layout checkpoint source implemented; broader closures/EM families and runtime/rank/restart validation remain open |
| 003E cluster scale/restart | Agent unassigned; benchmark/MPI/I/O | 003C/D; real multi-node allocation | Hardware campaign pending | Freeze strong/weak scaling and changed-layout restart workloads before requesting/running resources |
| 004A GPU primitive closure | Root; extensions and coupling | Current kernel gates and exact-SHA CI | Metal uniform push and NGP have green CI. Current gathered preflight source adds guarded canonical native publication and explicit field-storage admission; historical V1 allocation, V2–V5 latency and V6 calibration failures remain retained | Accept current source on GitHub and actual Metal, including native bits, compiler and both fixed cost matrices. The 1.5 latency and 8 KiB additional-allocation limits remain unchanged; [historical results do not accept this source](METAL_GATHERED_VALIDATION_PLAN.md) |
| 004B/C resident solver families | GPU architecture design agent; fields, integrators/models, extensions | 004A; actual device storage/operator/API contract | Concrete first-family contract under source review; no resident solver implemented yet | Freeze actual CUFFT/kernel/range/phase/private-public/transfer/full-persistence contract before source implementation; hardware/runtime acceptance remains open |
| 004D multi-GPU | Agent unassigned; MPI + device extensions | 003C/D, 004B/C; real devices/network | Hardware/integration pending | Establish buffer ownership/staging contract and matched device/rank oracle |
| 005A existing CN-labelled paths | Shock3d and integrator stability owners | Independent linear-system oracle | Shock CN `ce52836` and [linear CN](LINEAR_CN_RANGE_PLAN.md) `03e4eb3` have green exact-SHA CI | Continue with 005B/C; preserve separate resource regressions and corrections |
| Magnetic timestep correctness | Root; `recommended_dt` and both RK4 stage owners | Independent frozen Hall/damping mode and exact stability polynomial | [Certified polygon repair](RK4_MIXED_STABILITY_PLAN.md) passes both runtime campaigns and full compiler/resource gates | `a7fb3c8` retains its failed packet-allocation CI; the later containing `fd0a168` is green. The separate BigFloat repair has the revision-specific outcome below |
| Relativistic stored-output closure | Particle, priming, EMPIC, radiation and persistence owners | Unchanged supported cone; independently rounded references and whole public-step undo | [C1 local gates](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md) accepted, with original compiler failures retained under complete matched disposition | `1bc618b` passes all six exact-SHA CI jobs and both complete compiler reviews; [C2/C3](RELATIVISTIC_VALIDATION_PLAN.md#planned-domain-and-state-extensions) remain unimplemented |
| Original species charge/mass factors | Canonical kick, validation, prime/model/backend owners | Original q,m,dt and complete stored-input references | CPU public and nonrelativistic model kick threading, scalar/batch accuracy certificates and exact replay are implemented with the shared IEEE converter. Nonrelativistic priming/prediction include range/accuracy admission and exact endpoints. Both relativistic EMPIC owners retain factors through range replay and C1 publication; essential boundary, phase and checkpoint regressions are included | GitHub CI and complete resource/compiler acceptance remain open. [CPU implementation](BORIS_ACCURACY_IMPLEMENTATION.md); [CPU verification, normal relativistic accuracy and device follow-ons](PARTICLE_CHARGE_MASS_RANGE_PLAN.md) |
| Coordinate exact-range closure | Canonical coordinate owner and full/midpoint/native-batch callers | Exact stored-input overflow classification plus unchanged four-ulp/finite-bit contract | Canonical CPU correction passes 4,786 new and 46,725 legacy assertions per runtime; 108 paired resource profiles retain all 36 ordinary 0 B profiles; full compiler/Aqua, formatting and 3,543 actual-Metal assertions per runtime pass | [Contract, proof and costs](PARTICLE_KINEMATICS_RANGE_PLAN.md#coordinate-overflow-follow-on): both broader model lanes pass 117,819 / 117,820 assertions; correctness commit `a29da03` passes all six exact-SHA jobs and complete compiler review (1,246,963 / 1,246,964 assertions on actual Julia 1.12.7 / 1.13.3; whole forests byte-equal accepted C1). The separate zero-product optimization passes 7,366 focused assertions and all 108 paired profiles per runtime, including all 36 affected profiles at 0 B; actual Metal passes 3,543 and legacy 46,725 per runtime; complete compiler/Aqua captures match both parents and composition preserves every tested source byte. Optimization `8d7473f` passes all six exact-SHA jobs and complete compiler review: 1,249,543 / 1,249,544 assertions, with whole forests byte-equal accepted C1 |
| BigFloat timestep contexts | Root; canonical recommendation core | Independent positive-rate geometry, exact endpoint controls and production-mode witness | [Directed context repair](BIGFLOAT_TIMESTEP_CONTEXT_PLAN.md) passes isolated minimum/latest package, native bit/allocation, complete compiler/Aqua and formatting gates | `36d1794` CI completed: minimum Julia/formatting/MPI pass, latest Julia fails the separate packet allocation gate. Containing `fd0a168` is green with complete compiler review. NUM-001 runtime integration remains open |
| 005B/C implicit hybrid | Root integration; canonical integrator + moment/Ohm owners | Explicit represented Boris-map/Faraday and work-defect contract | Full bounded synchronous CPU isothermal collection core/report/checkpoint/archive source implemented; [guide](IMPLICIT_HYBRID.md) | Exact-revision numerical, residual/convergence, failure/restart, compiler and cost acceptance; broader pressure/orbit-averaged models remain open |
| 005D implicit EM / parallel geometry | Root integration of frozen EM/field/nonlinear/IO source | Actual Boris/current-continuity/Maxwell and nonlinear owners | Complete bounded coupled source reviewed/integrated, three legacy inverses and all eight frozen file hashes checked; essential remote regressions prepared | Exact-revision numerical/constraint/work, failure/restart, compiler/resource gates; broader conservative/relativistic/backend methods remain open |
| 006A curvilinear solver | Root integration of frozen mesh/coupling/model/IO source | Compatible geometry/component/volume conventions | Complete bounded periodic-shear core, diagnostics, loader and physical IO source reviewed/integrated | Exact-revision runtime, field/work/refinement, restore/compiler/resource gates; toroidal axis/walls remain separate |
| 006B self-consistent gyrokinetics | Root integration of frozen population/core/coupling/IO source | Reviewed full-f slab model and registered SO 0.1.8; 006A remains required for toroidal claims | Complete bounded source and canonical callers reviewed; essential remote regressions prepared | Exact-revision runtime, field/reciprocity/convergence, failure/restart, compiler and resource acceptance; broader geometry/physics remain open |
| 006C adaptive plasma hierarchy | Static composite-ES implementation worker; AMR/field/coupling/model/IO | Real patch/covered/leaf ownership and composite field method | Complete 373-line contract approved, all 33 source pins checked; implementation active | Complete source/IO review and remote field/crossing/self-force gates; dynamic refinement remains later |
| 007A/B reference provenance + production dispersion | Canonical validation/reference owners; next campaign agent unassigned | Source/network/toolchain; accepted production solver | Pinned provenance, parallel production histories and oblique Hall case 35 landed; CI state is in the audit tracker | Complete manual campaigns and broaden thermal/kinetic/EM regimes, uncertainty and matched-cost comparisons; [production plan](PRODUCTION_DISPERSION_PLAN.md) |
| 007C–E profiles, new physics, comparative costs | Agent unassigned; validation/benchmark owners | Accepted corresponding capability; datasets/hardware | Existing scalar comparisons inspected; campaigns unrun here | Freeze matched inputs/diagnostics/uncertainty and preserve raw reference outputs |
| 008 ES/EM arbitrary species | Canonical electrostatic/EM collection owners | Existing deposition, spectral Maxwell and phase owners | ES containing `4724947` passes all six jobs/compiler review. [EMPIC collections with generic/specialized dividing schedules](EMPIC_SPECIES_IMPLEMENTATION.md) implemented with shared ordinary-store admission, staged publication and persistence; new GitHub acceptance pending | Accept current sum/continuity, per-species phase, physical diagnostics, restart/migration and matched cost; nondividing/adaptive schedules and backends remain open |
| 009 full-Hall multidimensional shocks | Root integration of frozen Hall A and B source | Mixed operators and reviewed bounded field/reservoir/time/phase contract | Complete [numerical-reservoir core](OPEN_HALL_HYBRID.md) and [full physical/RNG IO](OPEN_HALL_HYBRID_CHECKPOINT.md) source reviewed; essential remote regressions prepared | Exact-revision runtime, boundary/convergence, failure/restart, compiler and resource acceptance; broader physics remains open |
| 010 composed sources and RNG restart | Integrator/checkpoint owners; independent collision audit track | 001 and phase/growth/checkpoint contracts | Fixed-step first-order isothermal CPU `CollisionalHybridStepper` landed; `2b2420b` exact-SHA CI/compiler gates pass. Ideal evolving, S1/S2 transport and S3b conduction extensions implemented, CI pending | Bounded generic collisional ES/EMPIC source is implemented; accept current extensions, then extend schedules, other closures/models/sources and backends; preserve full pressure/force/RNG/phase/rollback contract |
| 011 electron transport/energy | Electron/integrator owners and transport track | Term-specific balance and carried force contract | S1 is landed. The registered SO 0.1.8 consumer floor and bounded scalar original-product prerequisites are integrated; the dependency's five passing checks do not accept HPIC. Earlier 0.1.7 and rejected-candidate evidence remains historical | Exact-revision HPIC CI, transport resource/peak-memory and restart checks, general original-product admission and S3–S8 remain open |

For every row, its detailed section below supplies the required conservation,
accuracy/convergence, error/rollback, restart, parallel, hardware, performance and
independent-audit gates. Update this table and the top-level tracker together when
an owner starts, a candidate fails, a dependency changes or acceptance completes.

## CAP-001 — spatially local interspecies collisions

Implementation owner: the collision track. Canonical code remains
`particles/collisions.jl`, using periodic cell classification from
`particles/sorting.jl`; do not introduce a second scattering-angle algorithm.

The proposed first public slice is a collection/grid overload of
`collide_coulomb!` with explicit Coulomb logarithm, normalization, RNG and reusable
workspace. Retain the existing homogeneous single-set API. First scope is CPU,
nonrelativistic, periodic Cartesian 1D–3D, Float32/64; rejection of unsupported
backends is part of the contract. Equal macro weights conserve pair momentum and
energy using physical species masses; unequal weights require the selected
method's statistical conservation, rate and reuse-count proof.

The standalone CPU slice landed as `1b7b8b7` and passed all six exact-SHA CI
jobs and full compiler-tree review. The canonical rate/rotation owners include
the earlier null-angle, small-angle, population-mean and finite-draw repairs;
`cda63f1` subsequently adds normwise updates and hoists per-cell/reuse rate
coefficients. Its historical CI failure is retained; containing revision
`4724947` passes all six jobs and complete compiler review. The
[detailed local tracker](LOCAL_MULTISPECIES_COULOMB_PLAN.md) retains rejected
candidates and numerical/resource evidence. No global cost optimum follows.

- [x] Cell-local self/cross pairing with empty/singleton/odd populations,
  arbitrary species counts/weights, shared physical-rate arithmetic, alias
  rejection and staging before velocity publication.
- [x] Independent drift/diffusion, equilibration-rate, conservation, scaling,
  cross-cell isolation, fixed-RNG replay and legacy regression suites are part
  of the accepted CPU tests. Ordinary reusable-workspace calls are allocation-free.
- [ ] Broader complete-call cost campaign, extreme occupancy/species imbalance,
  production workload and independently controlled timing acceptance.
- [x] CAP-010 supplies fixed-step first-order isothermal CPU composition and full
  advanced-RNG/phase restart.
- [ ] Accept the implemented ideal evolving-scalar, S1/S2 transport and S3b conduction extensions
  through exact-revision GitHub CI; ion-ion collisions leave pressure fixed,
  while existing field stages provide the applicable Ohmic/hyperresistive heat and conduction.
- [ ] Accept actual-pencil [collisional ES](DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md)
  and [EM sources](DISTRIBUTED_COLLISIONAL_EMPIC.md): migrated populations pass an explicit canonical collision-cell ownership gate,
  use owned logical bins/global dx and persist full per-rank RNG/physical state.
  Actual-rank/statistical/restart acceptance and broader MPI families remain open;
  CAP-004 must still supply real device pairing/streams.

## CAP-002 — nonperiodic electromagnetic PIC

**Contract and owner.** Add a genuine bounded Maxwell discretization under the
existing EM model orchestration, with mesh operators in the mesh/operator owner,
trajectory current deposition in `coupling/`, physical boundaries in
`boundaries/`, and persistent boundary state in `io/`. The existing periodic
spectral constructor/algorithm stays available. A boundary flag on `FourierGrid`
that merely clips particles would not implement nonperiodic Maxwell physics.

**Proposed first discretization:** a Cartesian staggered Yee field layout and
charge-conserving trajectory current. Derive the compatible discrete divergence,
curl and particle interpolation before coding. PML requires actual auxiliary
field evolution and has finite discretization reflection; this follows from the
[WarpX PML equations and discretization](https://warpx.readthedocs.io/en/latest/theory/boundary_conditions.html).
Use existing boundary particle utilities only when their ownership and charge
ledger satisfy the coupled method.

- [x] **002A: Maxwell vacuum/cavity.** Dimension-parametric staggered state,
  periodic cross-check and conducting walls; explicit centering/CFL contract;
  compatible curl/divergence and field-energy/Poynting diagnostics. Gate cavity
  frequencies, analytic pulse transport, reflection phase and spatial/time order
  at multiple modes, incidence directions and corner configurations.
- [x] **002B prerequisite: closed longitudinal charge/current.** Extensive nodal
  CIC charge and endpoint path impulse are implemented for CPU Float32/64 with
  exact range replay and atomic output publication. Their scalar path voltage
  shares the same overlap geometry and stored-dx work metric. Optional finite
  physical-width triangular-cloud static coupling adds a continuous stationary
  work-adjoint field, with explicit particle-bound endpoint tails and distinct
  wall response. Native voltage/cloud arithmetic and resource acceptance remain
  open. GitHub verification is pending; see
  [the scoped contract](BOUNDED_PARTICLE_COUPLING_PLAN.md).
- [ ] **002B: bounded particle coupling.** Deposit the complete trajectory before
  removal/injection; split multi-cell and multi-face crossings, preserve IDs and
  species charge, account for wall surface charge/current. Gate discrete local
  continuity including boundary flux, Gauss propagation, particle-field work,
  reflecting-wall energy, exact free-flight crossing maps and empty populations.
  An event-free longitudinal electrostatic `YeeLongitudinalPIC` source candidate
  now supplies synchronized cloud Verlet coupling, exact background/free-plate
  initialization, diagnostics and distinct checkpoint/archive state. Its
  [scoped contract and GitHub gates](BOUNDED_PARTICLE_COUPLING_PLAN.md#event-free-longitudinal-coupled-source-candidate)
  do not accept transverse EM coupling or wall events.
  `YeeElectromagneticPIC` now adds a distinct event-free 1D3V source candidate:
  both polarizations, exact compatible transverse currents/PEC mask, synchronized
  Boris/Yee stages, scoped combined-rest/physical-H exit guards, diagnostics and
  complete schema-2 restart. [Contract and pending GitHub gates](BOUNDED_PARTICLE_COUPLING_PLAN.md#event-free-1d3v-electromagnetic-source-candidate)
  do not accept general nonlinear stability, wall events or broader bounded EM.
  Both coupled source candidates now retain complete rational original-geometry
  gathers, q/m and full-dt half through canonical Boris, rounding only final
  half-kick velocities. Their schema-2 force policy rejects the former schema-1
  rounded-gather map. Exact-revision runtime/compiler/format and complete-step
  resource acceptance remain pending; ordinary native bits are not preserved.
  The distinct [`YeeImpactLongitudinalPIC` source candidate](YEE_IMPACT_LONGITUDINAL.md)
  now implements nonrelativistic CPU Float32/64 finite-cloud reflection on a
  closed longitudinal 1D grid: zero contacts or one connected certified cluster per external step,
  original-factor forces and at most two rational-duration synchronized KDK
  chunks. Conservative carried state/source/time bounds check the full remainder
  from the original external-step state; ledgers, diagnostics and distinct restart
  stage before publication. Exact-revision GitHub numerical, independent
  force-driven convergence, rollback/restart, JET/compiler, formatting/CI and
  complete-step resource acceptance remain pending. Broader wall laws, multiple
  clusters, transverse walls, relativity and backend integration remain open.
  The distinct [`YeeAbsorbingLongitudinalPIC` source candidate](YEE_ABSORBING_LONGITUDINAL.md)
  adds zero or one uniquely isolated certified-passive instantaneous full-cloud
  capture per external step in the same closed nonrelativistic CPU Float32/64
  longitudinal 1D domain. It retains original-force rational-duration KDK,
  full death current, dynamic exact plate and capture ledgers, stable Vector
  compaction and distinct restart. Its survivor carry handles the finite force
  jump separately from reflection's fixed-source contract. Exact-head GitHub
  numerical/convergence, rollback/restart, JET/compiler, formatting/CI and resource
  acceptance remain pending. General material/sheath laws, multiple captures,
  transverse walls and backends remain open; **002B is not complete**.

- [ ] **002C: absorbing/radiating boundaries.** Add a measured absorbing layer
  with configurable thickness/profile and accounted energy removal; inject
  external fields/particles with explicit energy/charge budgets. Gate reflection
  spectra versus angle/frequency/resolution and long-time stability, including
  face/edge/corner layers and particle exits through them.
The separate [graded scalar-conductivity source](YEE_CONDUCTIVE_LAYER.md) now
provides D1–3 CPU fixed-step midpoint E/B damping around canonical Yee KDK,
exact actual damping loss, source work and signed field-map defect ledgers,
staged publication and a distinct complete nested restart.
The scalar layer is distinct from the new [causal directional ADE PML](YEE_PML.md),
which supplies separate target-lattice memories, original-duration midpoint KDK,
signed exchange/field-defect ledgers and complete actual-history checkpoint/archive.
Zero-prehistory init is explicit; general CPML tuning is not supplied.
Reflection/convergence/long-time stability, runtime/restart/compiler/resource
acceptance and coupled particle/source charge integration remain open. No tests
were added or run; 002C remains unchecked.

Serial [Yee vacuum checkpoint/archive](MAXWELL_VACUUM_CHECKPOINT.md) is now
source-implemented for synchronized field/clock state on matching grids. Its
restart/compiler/resource acceptance is deferred; this is a field-only part of
002D, with no absorbing-layer or particle-state substitution.

- [ ] **002D: production integration.** Relativistic and subcycled species,
  restart including layer/surface state and half-time fields; CAP-003 halos,
  CAP-004 resident kernels and external bounded-plasma comparisons under CAP-007.

Invalid boundary combinations, incomplete injection data and oversized CFL steps
must reject before mutation. Midstep recoverable failure needs full particle,
field, charge-ledger and phase rollback. Performance gates measure full steps,
trajectory splitting, absorber storage and boundary load imbalance; ordinary
periodic paths retain their baseline. Embedded curved conductors, dielectric
material models and sheath emission are subsequent explicitly specified slices,
not claims obtained from the first Cartesian wall implementation.

## CAP-003 — domain-owned distributed stepping

**Contract and owner.** `parallel/domain_decomposition.jl`,
`parallel/pencil_decomposition.jl`, `parallel/mpi.jl` and the Pencil extension own
partition/transport policy. Physics remains in existing moments, Ohm and
integrator owners. Each rank stores owned cells plus required halos, not the full
global field. The replicated reference remains a small-case oracle.

The [operator implementation plan](DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md)
specifies precision, shared spectral kernels, physical storage and staged operator
acceptance. It keeps distributed transport in HPIC and numerical Fourier policy
in SpectralOperators. Kernel availability in a locally developed dependency does
not close reproducible package resolution.

The [periodic vacuum Maxwell source](DISTRIBUTED_MAXWELL_IMPLEMENTATION.md)
composes the existing operators into actual rank-owned field evolution. It uses
synchronized spectral KDK, exact scalar affine rows over finite represented
curls, a strict spectral CFL bound, staged six-field publication and collective
energy/divergence diagnostics. Source review adds no runtime or scaling evidence.
[Matching-layout checkpoint](DISTRIBUTED_MAXWELL_CHECKPOINT.md) is source-implemented
with a distinct checksummed generation schema and detached all-rank preflight;
restart acceptance is deferred. The [static periodic Poisson field](DISTRIBUTED_POISSON_FIELD.md)
now also supplies a source implementation using the canonical mode owner, exact
represented charge admission and staged fields; runtime acceptance is deferred.
The [distributed ES source](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md) now
composes actual-pencil CIC anchor ownership, complete remote gather/deposition,
direct migration and the canonical leapfrog mover with that field solve. Exact
charge contributions route before density rounding; required explicit
`grid_dc=:project` and separate physical/grid ledgers avoid a hidden neutrality
assumption. Particles, density/E and clock stage before collective publication.
[Matching-layout ES checkpoint](DISTRIBUTED_ELECTROSTATIC_CHECKPOINT.md) now
retains complete particles, density/E and variable-step phase. Checksummed
generations and detached all-rank preparation precede publication; old destination
values may be stale and local counts may change. No field reconstruction occurs.
Its restart/liveness/compiler/resource acceptance is deferred. A separate
[distributed EM source](DISTRIBUTED_EMPIC_IMPLEMENTATION.md) now adds all-mobile
nonrelativistic CIC/Boris, corrected midpoint current, carried B/velocity phase,
exact scalar diagnostics and complete staged eleven-field publication.
[Complete matching-layout EM checkpoint source](DISTRIBUTED_EMPIC_CHECKPOINT.md)
retains those actual fields/current/endpoints and phase with detached all-rank
restore, deriving count cache from incoming leaves. The
[isothermal hybrid source](DISTRIBUTED_HYBRID_IMPLEMENTATION.md) now composes
positive-ion collections, strict density, canonical midpoint motion, frozen-moment
RK4 magnetic subcycles and endpoint prediction, with all eighteen physical fields
and particles staged before publication. Its scalar free-energy diagnostic keeps
canonical local log maps on each rank and reduces exact enclosure endpoints.
[Complete hybrid checkpoint source](DISTRIBUTED_HYBRID_CHECKPOINT.md) retains
all eighteen carried fields, full particles and actual clock/half-phase without
reconstructing predictor moments or E. Detached incoming state and complete
replacement/count preparation precede the final all-rank gate; stale destination
values and changed local counts are supported. Broader closures/EM families,
collision compositions, backends, runtime/restart acceptance and cluster scaling remain open.
No global field/population gather occurs in these source paths; CAP-003 remains
incomplete and its acceptance checkboxes remain open.

The accepted dependency decision is PencilFFTs for periodic 3D transforms. Its
[documented local/global index and permuted-output layouts](https://jipolanco.github.io/PencilFFTs.jl/dev/tutorial/)
must be respected when applying wavenumber multipliers; a serial wavenumber index
cannot be reused as a local pencil index without a mapping.

- [ ] **003A: layout/storage.** Explicit owned/ghost regions, cell volumes,
  staggered locations where relevant, communicator lifecycle and physical-boundary
  ownership. Gate unequal partitions, thin/empty allowed partitions, all directions,
  periodic seams and shape-radius halos; reject impossible layouts collectively.
- [ ] **003B: distributed operators.** Reuse FFT plans and implement derivative,
  curl, Laplacian, projection/filter and normalization in actual pencil layouts.
  Gate analytic modes, DC/Nyquist, divergence/curl identities, transpose/inverse
  parity and preallocated plan reuse on 2/4/8 real ranks. Mixed nonperiodic
  SBP/Fourier operators are a separate accepted backend, not a periodic shortcut.
- [ ] **003C: moments and full hybrid step.** Deposit local/ghost contributions,
  sum ghost moments once, apply species/rank summation before floors, then solve
  local fields and gather with complete halos. Migrate across any crossed domains
  before the next owner-dependent operation. Gate serial-reference trajectories,
  budgets, pressure, variable-step phases and collections across rank boundaries.
- [ ] **003D: production families.** Extend accepted composition to ES/EM,
  full-Hall shocks, collisions and open boundaries; attach rank-local diagnostics,
  partition-aware checkpoint and restart to each actual solver family.
- [ ] **003E: scaling.** Real multi-node strong/weak scaling and load imbalance
  campaigns, plus restart on supported changed layouts. Use no full-field gather
  in the measured hot path; account separately for explicitly requested output.

Collective validation must make all ranks take compatible communication paths;
test one-rank invalid inputs and asymmetric I/O failures without hanging. A lost
MPI process is not recoverable unless a separately implemented fault-tolerant
transport contract exists. Gate particle/field ownership, no duplication/loss,
global conservation and checkpoint publication. Measure maximum per-rank memory,
bytes communicated, reductions/transposes, migration volume and time-to-solution;
verify local storage scales with owned cells plus halos. Single-machine ranks
cannot close the cluster scaling gate.

## CAP-004 — complete GPU-resident solvers

The shared gathered-push [midpoint storage repair](GATHERED_MIDPOINT_STORAGE_PLAN.md)
passes all six exact-SHA CI jobs and complete compiler review at `a79a07b`. The
separate [field-input storage repair](GATHERED_FIELD_STORAGE_PLAN.md) is pushed
as `0f78536`; [all six exact-SHA jobs and complete compiler review pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40411920570). These correctness
gates do not establish GPU latency or complete device residency.

**Contract and owner.** Extend concrete storage parameters in `HybridFields`,
stepper/model workspaces and operator plans only as required by a real backend.
Reuse canonical physics and input validation; device methods belong in CUDA/Metal
extensions. `copy_particles_to_backend` alone does not move solver fields or
make host scalar loops valid device kernels.

- [ ] **004A: kernel contract closure.** Audit deposit/gather, push, wrapping,
  moments, pressure, sorting/growth and reductions on real hardware; include shape
  boundaries, empty/block-edge counts, invalid input/rollback and large index
  ranges. Existing Metal tests cover primitives; CUDA currently lacks analogous
  package particle/coupling kernels. Preserve honest precision support.
- [ ] **004B: one complete periodic solver.** Build storage/operator support and
  one complete hybrid stepping path with scalar algebraic closure, including
  initialization, diagnostics and restart. No per-step field/particle host copy,
  host scalar indexing or hidden CPU solve. If a backend cannot supply the chosen
  device FFT/operator, leave that solver/backend unsupported until implemented.
- [ ] **004C: full family coverage.** Add evolving/CGL closures, collections,
  CAM-CL, ES/EM including relativistic/current-correction paths, collisions,
  boundaries and all persistent phase state as separate composed slices.
- [ ] **004D: multi-GPU.** Combine CAP-003 ownership with measured direct-device
  or explicitly reported host-staged communication, device/RNG ownership, rank
  migration and restart. Hardware support must be detected and tested per path.

Gate CPU/device numerical parity at stated tolerance, analytic trajectories,
charge/divergence identities, time/grid convergence and ensemble collision rates.
Atomic deposition order needs a documented tolerance or a separately measured
deterministic mode, not an unsupported bitwise claim. Error/rollback tests must
synchronize at the relevant transaction boundary. Long repeated-step tests must
measure live device memory, pool growth, kernel/transfer time and complete-step
throughput with scalar indexing disabled. Real Apple GPU and NVIDIA GPU jobs are
distinct required gates for their respective claims; neither substitutes for the
other. Optimization requires matched physics accuracy and includes diagnostic/I/O
overhead. Multi-node GPU efficiency remains open until actually measured.

## CAP-005 — nonlinear coupled implicit integration

**Existing scope matters.** The analytic CN multiplier is a linear reference.
The [shock CN solve](SHOCK_CN_SOLVE_PLAN.md) now uses restarted GMRES for the
frozen-moment reduced induction equation, recomputes its true linear residual
and rejects failed convergence before public state publication. This replaced
the former four-iterate approximation. Those two field-only paths retain their
limited scope. A separate [coupled implicit hybrid source](IMPLICIT_HYBRID.md) now
implements synchronous CPU isothermal positive-ion collections in D1–3, using
canonical moments/Ohm/Boris/Faraday, bounded joint velocity/B Newton–GMRES, fresh
residual gates, staged endpoint/report publication and complete persistence.
Its represented work-defect report makes no exact energy theorem; numerical,
convergence, continuation, compiler and cost acceptance remain unexecuted.
The [Stanier–Chacón–Chen method](https://arxiv.org/abs/1803.07158) provides a
primary-source design reference for compatible cell-centered spatial operators,
implicit midpoint, ion subcycling/orbit averages and conservation. Its properties
must be re-derived for the selected implementation rather than inherited by name.

- [ ] **005A: audit existing implicit-labelled paths.** The shock path now
  supplies a converged, true-residual-checked frozen-moment CN solve with explicit
  tolerances/iteration limits and staged whole-step publication. Its
  [scoped source and verification record](SHOCK_CN_SOLVE_PLAN.md) retains the
  former approximation as history. Complete current-revision residual,
  amplification/timestep-range and assembled-system acceptance; align
  `recommended_dt` with the actual coupled resolution limits. This broader
  milestone remains open.
- [ ] **005B: residual and transaction.** The new synchronous family implements
  joint velocity/B represented residuals, frozen dimensionless block scales,
  absolute/relative tolerance, bounded Newton/GMRES/probe/backtracking budgets
  and private candidate ownership. Fresh final gates precede full publication;
  typed failures preserve public state and prior report by source dataflow. No
  timestep retry occurs. Independent numerical and transaction acceptance remain
  open; see the [complete contract](IMPLICIT_HYBRID.md).
- [ ] **005C: coupled hybrid.** Periodic isothermal CPU collections, full-time
  fields/particles and complete fixed-policy/report persistence are source-
  implemented. Strict density and canonical original-factor Boris/moment/field
  owners are reused. Evolving pressure, orbit averaging, selected dissipative
  terms and preconditioning remain open. Accept the actual spatial/particle
  residual, phase/order, work and momentum errors and restart before advancing
  this milestone; the represented defect report is not exact conservation.

Current bounded EM source: [ImplicitEMPICStepper](IMPLICIT_EM.md) supplies a
real joint velocity/E/B map with fresh endpoint charge/current in every residual,
separate physical constraints, four-piece work ledger and full synchronous IO.
It shares the bounded nonlinear owner with hybrid and preserves old numerical
owners by literal extraction. Exact-neutral all-mobile CPU CIC/odd periodic
scope is explicit; remote runtime/scientific/resource acceptance remains open.
This continuity-corrected gather/current pair does not inherit an exact energy
theorem. Conservative orbit averages, broader geometry and backends remain work.

- [ ] **005D: implicit EM and broader geometry.** EM requires its own coupled
  Maxwell/charge-current equations; it is not an automatic hybrid mode. Integrate
  compatible boundaries, domain ownership and device execution only after the
  corresponding explicit operator contracts pass.

Acceptance: independent Jacobian/residual checks; manufactured equilibria;
second-order time convergence; conservation error versus solver tolerance;
stiff whistler/ion-acoustic/EM dispersion including phase accuracy; moving cold
beam and multi-cell orbit tests; positivity/admissibility; deliberately exhausted
linear/nonlinear iterations; bad preconditioner and NaN residual rejection;
restart and variable-step equivalence. Include near-vacuum/floor activation and
identify where the physical model's conservation proof no longer applies.

Measure residual evaluations, particle pushes, linear iterations, preconditioner
setup/reuse, allocations and time-to-solution at the same physical error as the
explicit method. A larger stable timestep alone does not establish efficiency
or resolved kinetic dynamics. Freeze both stiff and nonstiff workloads; preserve
the explicit path's cost and report where the implicit method loses.

## CAP-006 — geometry, gyrokinetics and adaptive plasma evolution

These are three related implementation tracks, not one interchangeable mesh flag.
They share normalized units, diagnostics and persistence, but retain distinct
equations and invariants.

### CAP-006A: coupled toroidal/curvilinear plasma solver

Current bounded source: [ShearedPeriodicHybridPIC](SHEARED_PERIODIC_HYBRID.md)
evolves real full-f ions and common fields on one periodic nonorthogonal J=1
shear, with explicit physical/logical components, canonical complete TSC,
true continuous Lorentz rates, collective RK4/endpoint refresh, physical work/
energy/flux diagnostics and full synchronous physical IO. Source review included
the source-risk/Faraday, alias, phase and restart paths and repaired a receipt
validation gap. Essential GitHub checks are prepared. This does not establish
runtime/scientific/resource acceptance or the toroidal/axis/wall/equilibrium,
variable-J/imported geometry and backend requirements below.

- [ ] Define physical versus logical particle coordinates, vector/covector
  components, Jacobian-weighted densities, axis/outer-boundary regularity and
  equilibrium representation. Extend `curvilinear_mesh.jl` with the compatible
  curl/field solve and volume/face measures required by the selected hybrid model.
- [ ] Implement geometry-aware deposit/gather with a proved discrete work
  identity; evolve actual particle populations and fields on `ToroidalGrid`.
  Retain the existing metric utilities as tested components, not the stepper.
- [ ] Add equilibrium loading, boundaries, full restart and volume-weighted
  energy/flux/particle diagnostics. Then add CAP-003/004 ownership and devices.

Gate coordinate/inverse-map consistency, metric identities, uniform physical
fields, axis convergence, global particle number, weighted integration by parts,
orbit invariants, a genuine toroidal wave/equilibrium evolution and independent
equilibrium/orbit comparisons. Momentum conservation must follow actual geometric
symmetries: the [curvilinear implicit hybrid paper](https://arxiv.org/abs/2110.03886)
distinguishes energy conservation from geometry-dependent momentum conservation.
Measure metric caches, mapping/gather costs, memory and complete timesteps at
matched physical resolution. General stellarator equilibria and shaped separatrix
meshes require their own mapping/import/boundary slices before being advertised.

### CAP-006B: self-consistent gyrokinetic solver

Current bounded source: [SlabGyrokineticPIC](SLAB_GYROKINETIC.md) evolves actual
full-f gyrocentre collections and a common potential in a constant-B periodic
3D slab. It includes ring-TSC deposition/potential derivatives, full-k
polarization with modified adiabatic electrons, all collective RK4 stages and
endpoint field refresh, energy diagnostics, and full physical checkpoint/archive
continuation. Canonical source/callers were reviewed; the rounded-upper periodic
coordinate case and public archive restoration were completed during review.
Essential remote regressions are prepared, with no local Julia execution.
These source changes do not close the runtime/scientific/resource gates below
or supply toroidal geometry, evolving weights, collisions or electromagnetic
gyrokinetics.

- [ ] Choose and document a complete model before implementation: proposed first
  slice is electrostatic gyrokinetic ions with adiabatic electrons, explicit
  equilibrium, ordering and long-wavelength/polarization approximation. Introduce
  a distribution/marker population with owned weights, parallel velocity,
  magnetic moment and equilibrium parameters; the existing single `GuidingCentre`
  remains a local-field utility.
- [ ] Implement marker evolution, gyroaveraged deposition/interpolation and the
  matching quasineutrality/polarization field equation, with gauge/nullspace
  handling and a free-energy diagnostic derived from the same approximation.
  Full-f versus delta-f is an explicit model decision, not a runtime fallback.
- [ ] Add geometry and collision integration, restart of weights/equilibrium/RNG,
  then kinetic electrons and electromagnetic potentials with their own energy
  balance. Couple Raycon heating only through CAP-010's accounted source contract.

Gate uniform-field drifts, mirror/bounce orbits, finite-Larmor-radius ring response,
analytic linear modes, particle/free-energy balances, field-equation residual,
weight evolution, and refinement in time, space, velocity sampling and ring
quadrature. Run a matched Cyclone linear/nonlinear case against an independent
solver; [GX's published benchmark suite](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/gx-a-gpunative-gyrokinetic-turbulence-code-for-tokamak-and-stellarator-design/2C4BB81955E7E749B95B8B8141E997FA)
is a checked reference for such comparisons, not evidence this package passes.
Measure marker/ring/field costs and time to a fixed uncertainty in transport flux.
Reject invalid equilibrium or failed field solve before committing the population.

### CAP-006C: adaptive plasma hierarchy

The complete first static 1D composite-ES contract is approved after actual
source review and 33 source-pin checks. Implementation is active: true nested 2:1
patch/leaf ownership, one global mixed Neumann solve, fixed physical triangular
clouds, real coarse/fine crossings, signed work/rounding defects and complete
matching-hierarchy IO. The scalar-source helper moves byte-identically to its
shared validation owner; cloud/force/plate numerical owners stay in place.
This is an implementation stage, without runtime or adaptive-hierarchy acceptance.

- [ ] Add an owned hierarchy of patches, valid/ghost regions, refinement ratios,
  tagging/regridding and covered-cell accounting. Existing `AMRGrid` only stores
  one scalar 1D level; its transfer arithmetic remains the local oracle.
- [ ] First couple a static hierarchy to electrostatic PIC with a composite
  field solve, charge normalization, particle ownership and compatible gather.
  Then implement moving refinement, conservative remap and synchronized time
  evolution; do not claim adaptive evolution from a `refine` return value.
- [ ] Extend to hybrid/EM with flux/current synchronization, divergence and
  charge constraints, time interpolation/subcycling and coarse/fine correction
  derived for the actual solver. Add distributed partition/rebalance and device
  storage after serial hierarchy correctness.

Gate refinement/restriction conservation, manufactured composite solutions,
particle crossing and self-force at every patch face/edge/corner, no double
counting in covered cells, Gauss/divergence residuals, repeated refine/coarsen,
wave transmission and restart during each supported synchronization phase.
[WarpX's AMR analysis](https://warpx.readthedocs.io/en/latest/theory/amr.html)
documents particle self-forces and EM interface reflection as separate issues;
scalar interpolation tests do not settle either. Inject regrid/solve failures
and verify the old hierarchy remains usable. Compare cost and peak storage to a
uniform fine grid at the same physical error, including regrid and communication
cost; retain difficult interface trajectories in the fixed workload matrix.

## CAP-007 — independent validation and comparative evidence

Canonical owners are `validation/`, `verification/`, the independent test oracles
and reference metadata. Keep solver verification, model validity and external
performance comparisons separate in reported results.

- [ ] **007A: reproducible references.** Pin NHDS/Hybrid-VPIC source revisions,
  patches, compiler/build options, input decks, checksums, units, license and raw
  outputs. The [shared provenance implementation](EXTERNAL_REFERENCE_PROVENANCE_PLAN.md)
  replaces moving clones and unchecked cache reuse. Real pinned builds and all
  ten unchanged physics metrics pass on Julia 1.12.7 and 1.13.2; pushed as
  `95c5678`, with all six exact-SHA CI jobs and full compiler review complete.
  Missing/build-failed references stay
  failed/unavailable. This metadata does not establish a hermetic toolchain.
- [ ] **007B: production dispersion.** Keep the independent eigenvalue comparison
  in case 27, but also measure frequencies/damping/growth from actual initialized
  hybrid/Hall/EM time series and compare with the applicable independent theory
  or solver. Sweep modes, angles, temperature/mass ratios and refinement; freeze
  fitting windows, uncertainty estimation and branch selection before results.
  The [existing frequency helper repair](VALIDATION_FREQUENCY_PLAN.md) is a
  prerequisite: 3,153 registered checks, compiler/resource review and actual
  case 20 pass per runtime; correction `27372b2` has six green CI jobs and a
  complete compiler-tree review. This does not close the
  broader campaign. The [first production matrix](PRODUCTION_DISPERSION_PLAN.md)
  now passes all 324 metrics on both runtimes across 40 cold parallel profiles
  and six timestep-order histories. Corrected long-history fitting, actual test
  sandboxes, independent raw-data replay and CLI/plot integration pass locally;
  the slice is pushed as `839c109`. Minimum-Julia CI exposed twelve fitter
  residual failures; its other five jobs passed. The
  [observation-centering repair](MODE_HISTORY_CENTERING_PLAN.md) retains that
  evidence and all original thresholds. Oblique, thermal, kinetic and EM
  regimes remain open.
- [ ] **007C: nonlinear profiles.** Reproduce matched shock setups, including
  initialization, boundary conditions, closures and diagnostics. Compare full
  density/B/velocity/pressure profiles and time statistics, not just peak scalar
  agreement. Add full Preisser dataset replay with checksummed downloads and
  parameter mapping; do not infer it from the bundled summary roundtrip.
- [ ] **007D: expanded physics.** Benchmark bounded EM against an independent PIC
  solver, collisional equilibration against an independent kinetic reference,
  toroidal orbits/gyrokinetics against matched references, and AMR against uniform
  fine-grid solutions plus independent code. For reconnection, distinguish
  resistive versus collisionless mechanisms and measure reconnection rate,
  topology and energy partition beyond mode amplification.
- [ ] **007E: fair cost/accuracy campaigns.** Pin hardware/software and problems;
  compare accuracy/uncertainty, runtime, peak memory, scaling, setup and restart
  costs. Publish raw measurements and unsuccessful cases. Scope every comparative
  claim to the actual problems and machines; no single benchmark establishes
  worldwide superiority.

The current reconnection helper defaults to `η=0.005` and constructs an
isothermal model with that resistivity (`instability_sweep.jl:269,305–306`);
case 32 checks coherent-mode growth versus a uniform control (`run.jl:30–57`).
**Verified scope:** those checks alone do not establish collisionless reconnection
rates. **Hypothesis (untested):** the claimed mechanism is contaminated or driven
by explicit/numerical resistivity. Settle it with an independently diagnosed
generalized-Ohm balance, resistivity/electron-inertia scans, grid/time/particle
convergence and matched external cases. Do not state a causal result beforehand.

Case 27's [README](../validation/27_nhds_dispersion_comparison/README.md) and
`case_27_nhds_dispersion_comparison` now agree: unavailable references return a
failed availability metric, and malformed numeric data can raise a parsing
error. These failures are not skips. Source and prose were checked again
2026-10-29; this check does not rerun the external physics campaign. Keep
historical measured claims only with recoverable commands, source hashes and
outputs.

## CAP-008 — arbitrary kinetic species in ES/EM

This extends the documented current species restriction; it is not a defect in
that restriction. Generalize the canonical ES/EM charge/current and workspace
owners to ordered populations with arbitrary finite charge, positive mass and
independent storage. Preserve scalar/two-population constructors through the
same implementation, with explicit background/neutrality policy and per-species
time levels. Never normalize away nonzero physical net charge silently.

- [x] ES collections, shared Poisson source and charge/momentum/energy diagnostics.
- [ ] Accept [synchronous CPU EM collections](EMPIC_SPECIES_IMPLEMENTATION.md):
  summed charge-conserving currents, mobile/fixed populations, optional relativistic
  motion, physical diagnostics and restart are implemented; exact-revision GitHub
  checks and full compiler review remain.
- [ ] Accept [specialized EMPIC1D synchronous collections](EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md):
  S1 source is implemented and reviewed, with the actual Ex/Ey/Bz solver,
  invariant planar inputs, staged publication and distinct persistence.
  Exact-revision GitHub runtime/compiler/format acceptance remains pending.
- [ ] Accept per-species dividing integer schedules: specialized S2 and the
  explicit generic D1–3 opt-in share velocity histories, future setters and
  every-fast drift/current with independent magnetic history. Both have complete
  schema-2/policy-1 restart and explicit tagged schema-1 migration; default generic
  state/schema 1 stays. Specialized 59 and generic 33 text/source checks and
  root source/caller review are complete; runtime, physical refinement, compiler,
  formatter and resource gates remain pending. No new tests were added or run
  for the generic extension.
  Broader nondividing/adaptive/asynchronous schedules remain open.
- [ ] Accept new EM collection checkpoint/archive, including generic/specialized
  schedule phase/migration; extend CAP-001 interactions, CAP-002 boundaries, CAP-003 partitions and
  CAP-004 residency. Existing ES collection persistence remains landed.

Gate scalar/two-species compatibility, species permutation of physical totals,
pair-plasma and unequal-charge/mass dispersion, neutral/empty species, charge
continuity at every substep, total field energy counted once, initialization and
variable-step phases, alias rejection and complete restore. Compare storage and
complete-step cost versus species count/imbalance and existing two-species path.

## CAP-009 — full-Hall multidimensional open-boundary shocks

Current 2D/3D shock induction deliberately uses motional/resistive fields. Keep
that reduced-model behavior explicitly identified while adding a separately
specified full hybrid formulation in the canonical shock/Ohm/field owners.
For general transverse variation, derive the induction electric field and
compatible SAT/open-boundary energy flux; do not simply copy the particle field
into Faraday without checking the discrete constraints and boundary terms.

The [mixed Hall prerequisite](MIXED_HALL_OPERATORS.md) is now source-implemented:
D2/3 SBP-x/periodic-transverse cached operators, represented Euclidean transpose,
raw isothermal full-Hall induction and exact stored-source face-work ledger.
Explicit field/rate projection structurally removes normal face DOFs and uses
a consistent singular Q-metric multiplier construction with fresh equation,
direct constraint and signed work/energy gates. All represented nulls remain.
Native operator symmetry, convergence and cost remain unaccepted. That operator
slice owns no advancing state. The separate [OpenHallHybrid](OPEN_HALL_HYBRID.md)
source now supplies actual SBP-volume CIC coupling, both face reservoirs with
full-time newborn/first-exit transport, explicit SAT/projection work, bounded
paired RK4 trials, synchronous original-factor Boris velocities, signed energy
intervals and complete atomic physical/source/RNG persistence. Core and IO were
reviewed together, including last/cumulative per-face exchange records; essential
remote regressions are prepared with no local Julia execution. Runtime,
scientific-boundary, compiler and resource acceptance remain open. This is a
qualified perpendicular numerical-reservoir model, with no transparent Hall
outflow or general boundary-energy theorem. Evolving pressure, transport,
inertia and broader boundary/backend compositions remain separate work.
The existing reduced shock SAT remains outside curl; tensor derivative
commutation alone does not establish divergence preservation for that update.

- [ ] Derive and implement Hall/pressure/inertia terms supported by the declared
  model, compatible mixed operators and accurate timestep/implicit solve control.
  The massless isothermal raw RHS/operator/projector source is present; its
  runtime acceptance remains open; the separate bounded controller source is
  implemented. Pressure evolution,
  transport and inertia require separately declared boundary/model extensions.
- [ ] Couple upstream/downstream particle reservoirs and fields, pressure
  evolution, actual source/flux diagnostics and restart phase.
- [ ] Extend to distributed/device execution through accepted CAP-002/003/004
  operator and ownership contracts; validate multidimensional shock structure.

Gate 1D reductions, transverse whistler propagation, manufactured induction,
divergence behavior including SAT faces, upstream equilibrium, energy/charge/mass
flux balances, boundary-reflection sensitivity, rippling and refined shock
profiles. Include arbitrary transverse perturbations, not only uniform reductions.
Measure increased stiffness and error-matched throughput; the reduced formulation
is not a valid speed baseline without reporting its different physics.

## CAP-010 — composed sources, collisions and reproducible restart

The [bounded generic EMPIC composition](COLLISIONAL_EMPIC.md) now has a source
implementation: fixed-step field-first local Coulomb scattering at new positions
and carried half-time velocities, consistent unit field/collision normalization,
private field/particle/RNG staging and distinct complete restart. It retains the
canonical field and collision owners. No tests were added or run; runtime,
stochastic/physical, rollback/restart, compiler/CI and resource acceptance remain
deferred. [Collisional electrostatic PIC](COLLISIONAL_ELECTROSTATIC.md) now also
source-implements generic CIC ES composition with a sole wrapper clock, private
state/RNG staging and complete restart; its acceptance is deferred. Specialized
schedules, relativistic collisions, additional sources and broader backends remain
separate work. [Actual-pencil collisional ES](DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md)
now supplies fixed-step canonical owned-cell scattering with per-rank RNG staging,
sole inner clock/source ledger and complete one-generation physical/RNG restart.
Explicit cell ownership admission and equal field/collision normalization apply.
[Actual-pencil collisional EM](DISTRIBUTED_COLLISIONAL_EMPIC.md) now adds the same
local scheduler after the complete private EM candidate, preserving all eleven
fields and actual drift-current/charge history. A distinct complete physical/RNG
generation retains its sole clock and carried half phase. Cell/RNG/graph/busy and
prepared-generation owners are shared with ES; existing source paths reconstruct
by literal expansion. Actual-rank/statistical/rollback/restart/compiler/resource
acceptance remains open.
This does not complete CAP-010.

Add an explicitly owned evolution composition in the integrator layer; reuse
collision, ionization, radiation and antenna owners rather than embedding copies
in each model. Define the order and physical time of each source operation,
charge/energy/momentum exchange, particle growth, ID allocation and recomputation
of fields. First implement the accepted CPU collision/field split, then add other
sources one by one with their own conservation/accuracy contracts.

- [x] First slice (landed as `2b2420b`, exact-SHA CI/compiler review passed): `CollisionalHybridStepper`
  composes fixed-step CPU `HybridStepper` runs (ideal isothermal electrons) with
  cell-local collisions between the field advance and the carried-E formation,
  persists the advanced RNG and source ledger in checkpoints/archives and rolls
  back the whole step including the RNG. [Guide and verified state](COLLISIONAL_HYBRID.md).
- [ ] Accept the ideal `EvolvingElectrons` extension: implemented on the same
  source schedule, with pressure/RNG rollback and schema-1 consistency checks;
  focused component, zero-log, restart and failure tests pass in the
  [scoped remote run](COLLISIONAL_HYBRID.md#verification-status); full CI remains pending.
- [ ] Accept S1 `TransportElectrons` composition (`eta>=0`, `etaH=de2=0`):
  optional outer force snapshot, canonical inner policy-1 persistence and
  pressure/force/RNG regression tests are implemented; exact-revision CI pending.
- [ ] Accept the S2 admission extension (`eta>=0`, `etaH>=0`, `de2=0`): it reuses
  the same full-field/pressure/force/RNG owners and policy-1 persistence, with
  focused positive-hyperresistivity regressions passing in that scoped run;
  full exact-revision CI remains pending.
  Anisotropic/inertial composition and electron-ion heat exchange remain open.
- [ ] Persist complete RNG state/stream assignment, split phase, source clocks,
  ionization ID counters and external reservoir state. A saved initial seed alone
  does not specify continuation after random draws.
- [ ] Make rollback semantics explicit for state, counters, RNG and externally
  effectful callbacks. A source failure cannot commit only half a species pair.
- [ ] Add accounted Raycon power deposition only after specifying coupling time,
  receiving reservoir and conservation; standalone ray tracing remains usable.
- [ ] Integrate MPI cell ownership and real device RNG/particle growth when those
  backends are accepted. Specify deterministic versus statistical cross-layout
  reproducibility and test that distinction.

Gate splitting-order convergence, variable-step and zero-step behavior,
uninterrupted/restart equality in the same environment, source energy/charge
accounting, injected failures at each phase, unique IDs after migration/growth,
and bounded long-run allocation. Profile full composed steps, population growth
and output/checkpoint overhead; unused optional sources must add no per-particle
work. Broader reaction networks, differential cross-section secondary spectra and
full radiation models need separately derived physics contracts; the existing
reduced radiation operator is not relabeled as those models.

## CAP-011 — electron transport and complete energy accounting

Extend the current ideal evolving-pressure contract through explicit models,
starting with resistive/hyperresistive heating and their matching field losses,
then a specified heat-flux closure, then anisotropic pressure evolution and
selected electron-inertia terms. Canonical owners remain `electrons/`, hybrid
model/integrators and conservation diagnostics. Algebraic closures must retain
their present meanings; new dynamics require owned state and restart.

- [ ] Derive continuous and discrete reservoir/work/flux balances for each added
  term, including density floors, boundary work and source coupling.
- [ ] Implement coupled stage updates and admissibility checks with transaction
  rollback; add transport timestep/implicit residual limits and boundary data.
- [ ] Integrate collections, MPI/device operators, checkpoint/archive and
  time-centered diagnostics only for supported combinations.

Gate homogeneous relaxation/heating, manufactured transport, exact zero-term
reduction, thermal equilibration, positivity/admissibility, total-energy closure
versus refinement/solve tolerance and restart. Do not repair negative pressure by
undocumented clipping. Measure operator/stage/solver costs and transport scaling
against the frozen ideal baseline. A fully kinetic electron distribution remains
the ES/EM/GK model's responsibility, not an implicit property of a fluid closure.

**S1 status (constant Ohmic, serial CPU).** `TransportElectrons` splits Ohm's
law into the carried ion force `E_force` and `E = E_force + ηJ`, heats the
stored pressure by `(γ−1)η|J|²` from each RK4 stage current, rejects density
below `nfloor` and, in S1, unsupported `ηH`/`de2`/scalar charge, rolls back the carried
force with the pressure transaction and persists it as
`transport_policy_version = 1` (archive schema 4 unchanged). Tests cover
direct-DFT stage rates, the semidiscrete Ohmic exchange, analytic Beltrami
heating (fourth-order energy closure with fixed ions), second-order coupled
common-time energy, bitwise `η = 0` reduction, admissibility/rollback,
range-safe heating and restart ([guide](EVOLVING_ELECTRONS.md#constant-ohmic-transport)).
MPI, CAM-CL and Hall-MHD reject the closure. Thermal equilibration,
S2 acceptance, conduction (S3–S5), backends (S6), anisotropy (S7)
and inertia (S8) remain open. S1 landed as `a6d20d3` with all six exact-SHA
CI jobs and complete compiler review passed. `CollisionalHybridStepper` admits
ideal isothermal/evolving scalar electrons and S1/S2 transport or S3b conduction
with `de2=0`; the stored-pressure extensions await exact-revision CI. The transport wrapper preserves
the independent force in outer rollback and transport policy 1 at restart. Hyperresistive
and original-quotient conduction compositions reuse those same state owners;
their exact-revision CI remains pending.
The subsequent [Ohmic diagnostic range repair](RESISTIVE_DIAGNOSTIC_RANGE_AUDIT.md)
is pushed as `a356c6b` and passes all six exact-SHA CI jobs and complete compiler review.

**S2 current-main candidate.** The hyperresistive implementation and complete represented-source repair
are composed with the shared IEEE and CPU particle-factor changes without
replacing newer unrelated source or documentation. The shared normalized
Fourier API is supplied by SpectralOperators; its registered 0.1.8 consumer floor
adds certified directional axis-null conditioning. Exact-revision HPIC acceptance
is pending in the
[consumer tracker](SCALAR_PRESSURE_RATE_RANGE_PLAN.md#spectraloperators-018-consumer-adoption).
The prepared earlier HPIC
snapshot passed bounded 504-check and 48-profile ordinary-bit probes per runtime;
these historical results do not accept this composition. Minimum/latest Julia package/compiler, formatter and MPI validation runs in
GitHub, as recorded by the [S2 implementation tracker](TRANSPORT_SOURCE_RANGE_PLAN.md).
Public physical lapJ/E and complete required rates retain representability
requirements; normalized-dt RK storage remains unsupported. Raw Ohm partial-work
and coupled init/step rollback remain separate contracts. The
[original-factor ideal-force work order](IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md)
remains incomplete. Its scalar/CGL pointwise range prerequisite now combines complete
represented factors on intermediate overflow/subnormal products, with essential
GitHub checks pending. Bounded scalar Hybrid Faraday with `de2=ηH=0` also retains
original represented constitutive products on finite product-range events,
through the existing dyadic FFT owner. Normal swallowed variation, severe FFT
conditioning and broader source-aware operator/backend admission remain open.
Linear completion of represented ideal E alone cannot recover an original
product lost before differentiation. CAP-011 remains incomplete.

**S3a predecessor (represented temperature).** The initial source candidate
used `L(round_T(Q(pe)/Q(n)))` and conduction policy 1. Its dated source/CI
evidence remains history. S3b replaces that map and explicitly rejects its
checkpoint policy rather than reinterpreting it.

**S3b source candidate (original pressure/density quotient).** Active constant
isotropic conduction now retains exact original stage `Q(pe)/Q(n)`. Each axis
removes its own exact null background, prepares two normalized IEEE terms with
residual at most `u²*abs(input)`, applies the full stored-k/Nyquist operator through
the shared dyadic FFT owner, and subtracts the computed rational line mean.
Only the complete signed pressure source is rounded and range-classified.

The preparation and projection have explicit mathematical bounds; FFTW accuracy
still requires independent remote qualification. Severe cancellation and final
range thresholds retain that limitation. The represented gamma-minus-one owner,
finite stored rates/stages, pressure/density/rollback and scoped pure-conduction
guards remain. No positivity, nonlinear stability, exact energy or speed theorem
is inferred. See the [canonical model and checks](EVOLVING_ELECTRONS.md#constant-isotropic-conduction-source-candidate).

- [x] Implement exact quotient storage, fresh per-axis null removal, bounded
  dyadic preparation and post-transform rational projection in the canonical
  pressure path. General quotients are prepared before entering the unchanged
  dyadic-only transform; no arbitrary rational is silently accepted by it.
- [x] Compose both hyperresistive branches, scalar/collection and fixed
  collisional consumers. Preserve zero-κ reduction and derived-scratch retry;
  required conduction policy 2 rejects policy 1. Transport policy 1, collision
  outer schema 1 and archive container 4 retain their meanings.
- [x] Add compact independent original-variation/null/Nyquist, quotient-range
  rescue, complete signed-source cancellation, rejection/retry and policy/resume
  controls while retaining existing tolerances and unsuppressed compiler gates.
- [ ] Accept numerical/compiler/format and integration checks on the exact
  GitHub revision; measure complete-step BigInt/FFT work, allocation and memory.
- [ ] Extend conditioning qualification and broader transport models/backends.
  Normalized-duration rate storage remains a separate implementation gate.

## Additional risks requiring tests, not assertions

| Candidate | Verified source observation | Exact settling work |
| --- | --- | --- |
| Converged shock field solve may still leave coupled energy/accuracy limits | The current shock CN owner independently recomputes its true linear residual; the public CN path stages the full trial and checks finite final fields before publication. The reduced SBP/SAT and particle splitting remain separate contracts | Recheck the residual/rollback gates on the exact integrated revision, then assess whole-step temporal convergence, boundary work, diffusion energy and particle/field exchange; a converged linear solve alone is not a global energy-stability proof |
| Backend integration may introduce hidden host copies or unsupported arithmetic | fields/workspaces currently `Array`; GPU helpers cover only selected primitive calls | Real-device complete-step profile with scalar indexing disabled, transfer traces, allocation/memory ledger and all exceptional validation paths |
| Local weighted cancellation may remain despite correct global charge | Existing roadmap:690 records an unresolved local moment accuracy assessment | Independent exact per-cell deposition oracle across shapes/weights; serial/thread/device comparisons and coupled non-DC field response |
| Existing validation numbers may not replay on current dependencies | Source has historical numerical prose; this pass did not rerun campaigns | Immutable source/dependency replay with raw outputs, checksums and fitting scripts; classify unavailable artifacts rather than copying claims |
| Arbitrary exponent-span fallback resource cost | Existing exact-arithmetic plans retain unresolved large-range resource scope | Bound or measure complete public owners over declared supported exponent/precision ranges, including cancellation; no guessed resource ceiling |

## Progress and release evidence

- [x] Read immutable source, README, current roadmap and relevant design/test owners.
- [x] Separate implemented solver contracts, partial building blocks, missing
  integration and untested risks; identify CAP-008–011 as additional scope
  proposals beyond the coarse roadmap.
- [x] Expand CAP-001–011 into contracts, dependencies, implementation slices,
  correctness/failure/restart/parallel gates and fixed measurement strategy.
- [x] Root reviewed the complete plan and integrated its owners/stages with the
  live audit tracker and README; implementation acceptance remains separate.
- [ ] Execute each implementation slice and replace its unchecked gates with
  exact source, commands, outputs, hardware, cost matrix and CI evidence.
- [ ] Run the final compatibility matrix: model × geometry × closure × species ×
  source × precision × integrator × backend × restart. Cover each supported
  intersection and reject unsupported ones before state mutation.
- [ ] Independent full audit of all public paths and complete final validation;
  revisit the matrix after the last implementation, not merely the last failing test.

For each completed slice, record: parent/candidate/accepted SHA; affected public
contract and owners; baseline witness; independent oracle and mutants; focused,
neighbor, full-suite and compiler results; failure/restart/parallel evidence;
fixed performance results including rejected runs; real hardware configuration;
exact-SHA CI URL; remaining limitations. A plan checkbox, test count or green
CPU-only CI run alone cannot close a missing physics or hardware capability.
