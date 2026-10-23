# Plasma simulator completeness and implementation plan

Source assessment: **2026-10-23**, revision
`bf2a8bc615f514f8d2868934849bf44814140bbb`. This document plans work; it does not
certify the proposed implementations. The [audit tracker](DEEP_DEBUG_TRACKER.md)
owns accepted fixes, measured results and exact-revision CI status.

**The package is not yet a complete implementation of the expanded simulator
goal.** It contains implemented solver families and an extensive test suite,
but several advertised domains still lack an end-to-end evolution path.
Source presence, a passing primitive test and a complete simulation are separate
milestones. No current evidence establishes that the package is globally bug-free,
fully optimized, or stronger than every external plasma code.

The original scope recoverable from the repository is the periodic hybrid model,
its companion solvers and the open completeness roadmap: [README](../README.md),
[roadmap](DEEP_DEBUG_TRACKER.md#completeness-roadmap),
[distributed FFT decision](distributed_fft_evaluation.md),
[multispecies design](MULTISPECIES_HYBRID_PLAN.md),
[electron energy contract](HYBRID_ENERGY_MODEL_PLAN.md) and the two
[ray-tracing](superpowers/specs/2026-07-05-ray-tracing-design.md) /
[Raycon](superpowers/specs/2026-07-05-raycon-port-notes.md) specifications.
`src/HybridPlasmaPIC.jl:7` mentions an implementation checklist in the repository
root; no such file is present at this revision. This plan does not reconstruct an
unavailable specification from memory.

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

Line references below are anchored to the revision above. Recheck them when a
dependent slice changes the owner. Search methods, exports and callers as well
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
| ES/EM PIC | `electrostatic.jl:14–64`; `full_pic.jl:1747–1813`; `test_espic.jl`, `test_empic.jl`, `test_empic1d.jl` | Periodic CPU evolution exists; arbitrary kinetic species, nonperiodic Maxwell and resident GPU/MPI solvers are missing |
| Hall-MHD | `HallMHDModel`, `HallMHDState`, `step_hall_mhd!`; `test_hall_mhd.jl`, `test_io_hall_mhd.jl` | A periodic solver exists; it is not a general geometry/closure solver |
| Coulomb collisions | `collide_coulomb!` in `particles/collisions.jl`; single-population public contract | Local interspecies operator under CAP-001 implementation; composition/restart and parallel support are later gates |
| Nonperiodic EM | `EMPIC` constructor requires `FourierGrid`; `full_pic.jl:2199` wraps ions, electron path likewise periodic | Missing coupled field/particle boundary implementation, CAP-002 |
| Distributed evolution | `mpi.jl:934–960,1008–1027` explicitly deposits/reduces full grids and steps replicated fields; Pencil extension:31–60 wraps FFTs | Infrastructure and reference evolution exist; domain-owned solver missing, CAP-003 |
| GPU | Metal extension has deposit/gather/push methods; CUDA extension:1–68 has storage/telemetry hooks; `HybridFields` uses `Array` at field_state:12–29 | Partial kernels, not complete resident solvers, CAP-004; device availability is not numerical acceptance |
| Implicit evolution | `semi_implicit.jl:57–96` evolves linear Fourier amplitudes; shock3d:389–404 performs four Picard field iterations | Two limited paths exist; nonlinear coupled particle/field solve and solve-failure contract missing, CAP-005 |
| Geometry | `ToroidalGrid`, metric gradient/divergence in curvilinear_mesh:20–49,90–235 | Geometry operators exist; no coupled toroidal plasma stepper, CAP-006A |
| Gyrokinetics | `GuidingCentre`, local-field push and `gyroaverage` in gyrokinetics:559–645,823–876 | Single-particle/local sampling utilities; no self-consistent gyrokinetic distribution/field solve, CAP-006B |
| AMR | `AMRGrid` is a single 1D scalar level, adaptive_mesh_refinement:8–23; prolong/restrict/refine:322–366 | Transfer primitives; no adaptive plasma hierarchy, particle crossing or coupled timestep, CAP-006C |
| External validation | Cases 27/28 contain external-code adapters; case 27:85–88 evaluates a test eigenvalue oracle; case 28:94–155 runs the shock solver and compares scalars | Useful partial comparisons; production-solver dispersion, full profile replay and fair cost/accuracy comparisons required, CAP-007 |
| Kinetic species diversity | EMPIC docs fix electrons to `q=-1,m=1` and at most one `q=+1` mobile ion, full_pic:1751–1759 | Broader ES/EM collection evolution missing, CAP-008 |
| Multidimensional shocks | shock3d:9–19,337–359 and shock2d:230–236 intentionally exclude Hall induction while retaining Hall particle forces | Documented reduced shock model; full-Hall multidimensional branch required for broader claims, CAP-009 |
| Coupled sources and collisions | Core hybrid/CAM-CL/EMPIC steppers have no calls to collision, ionization or radiation operators; checkpoint:1 names deterministic state | Standalone operators exist; reproducible composed evolution missing, CAP-010 |
| Electron transport | `EVOLVING_ELECTRONS.md:3–10` restricts evolving scalar pressure to ideal massless parameters | Ideal scalar pressure exists; dissipative/transport-consistent energy extension missing, CAP-011 |

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

| Slice | Implementation owner / canonical code | Dependencies | State and checked evidence | Next concrete action |
| --- | --- | --- | --- | --- |
| 001 local self/cross collisions | `/root/local_multispecies_collisions`; `particles/collisions.jl`, `particles/sorting.jl` | Existing particle/moment arithmetic gates | Isolated implementation active; no acceptance claimed here | Complete normalization, independent rate, rollback, legacy compatibility and resource gates |
| 002A/B bounded Maxwell + particles | Agent unassigned; `models/full_pic.jl`, mesh/operator and coupling owners | Shared storage/layout decision; current periodic EM remains oracle | Planned; constructor and periodic calls checked | Derive staggered discrete operators and wall charge ledger; build vacuum/trajectory baseline before editing |
| 002C/D open boundaries + production integration | Agent unassigned; boundaries, EM, I/O | 002A/B; 003/004 for parallel/device claims | Planned; no coupled open Maxwell path | Pre-register pulse reflection/charge-flux suite and absorber state schema |
| 003A/B local storage + distributed operators | Agent unassigned; MPI/domain/pencil owners, Pencil extension | Storage/layout contract; existing real-MPI infrastructure | Planned; replicated moments and FFT wrapper inspected | Specify index/halo ownership and implement analytic distributed-operator tests |
| 003C/D distributed coupled solvers | Agent unassigned; MPI orchestration + canonical solver owners | 003A/B and accepted serial solver | Planned; full-grid allreduce reference inspected | Compose a periodic hybrid step without replicated grids and test 2/4/8 ranks |
| 003E cluster scale/restart | Agent unassigned; benchmark/MPI/I/O | 003C/D; real multi-node allocation | Hardware campaign pending | Freeze strong/weak scaling and changed-layout restart workloads before requesting/running resources |
| 004A GPU primitive closure | Root; extensions and coupling | Current kernel gates and exact-SHA CI | Metal uniform push and NGP pushed; gathered V1 passes behavior but fails allocation gate; V2 active | Complete public V2 hardware/compiler/cost gates, then enumerate remaining CUDA/Metal primitive paths |
| 004B/C resident solver families | Agent unassigned; fields, integrators/models, extensions | 004A; storage/operator contract; actual device operators | Planned; concrete host arrays inspected | Implement one complete resident hybrid path and trace all transfers/scalar access |
| 004D multi-GPU | Agent unassigned; MPI + device extensions | 003C/D, 004B/C; real devices/network | Hardware/integration pending | Establish buffer ownership/staging contract and matched device/rank oracle |
| 005A existing CN-labelled paths | `/root/completeness_architecture`; shock3d and integrator stability owners | Independent linear-system oracle | Confirmed residual/solution failures on both runtimes; repair active | Implement residual-checked solve and full-step rollback; retain direct-solve oracle |
| 005B/C implicit hybrid | Agent unassigned; integrator + moment/Ohm owners | Derived discrete conservation and solve contract | Planned | Freeze unknowns/residual, independent small-system oracle and failed-solve transaction |
| 005D implicit EM / parallel geometry | Agent unassigned; EM/parallel/geometry owners | 005B/C method evidence; accepted relevant operators | Planned | Derive EM-specific continuity/work residual before adding a public method |
| 006A curvilinear solver | Agent unassigned; mesh/coupling/model/I/O | Compatible geometry operators, chosen physical model | Partial primitives; no production callers found | Specify coordinate/component/volume conventions and coupled manufactured problem |
| 006B self-consistent gyrokinetics | Agent unassigned; gyrokinetics/coupling/field/I/O | Explicit model choice; 006A for toroidal claim | Partial primitives; no distribution/field integration found | Fix electrostatic model, marker/weight state and polarization equation |
| 006C adaptive plasma hierarchy | Agent unassigned; AMR/coupling/model/I/O | Composite field method; 002/003/004 only where claimed | Partial scalar-level primitives | Design static ES hierarchy and crossing/self-force oracle before dynamic refinement |
| 007A/B reference provenance + production dispersion | Agent unassigned; validation/reference owners | Source/network/toolchain; accepted production solver | Moving clones, unchecked cache reuse and eigenvalue-only comparison verified; READMEs reconciled | Pin references and define time-series fitting/refinement contract |
| 007C–E profiles, new physics, comparative costs | Agent unassigned; validation/benchmark owners | Accepted corresponding capability; datasets/hardware | Existing scalar comparisons inspected; campaigns unrun here | Freeze matched inputs/diagnostics/uncertainty and preserve raw reference outputs |
| 008 ES/EM arbitrary species | Scope proposal; root selects implementation agent | Species charge/time-level design; existing EM current owner | Current two-species restriction verified | Fix compatibility/background/multirate contract, then implement ES collection slice |
| 009 full-Hall multidimensional shocks | Scope proposal; root selects implementation agent | Full mixed-operator/field-boundary balance; 005 optional | Deliberate reduced induction verified | Derive full induction and SAT budget; build transverse-mode witness |
| 010 composed sources and RNG restart | Scope proposal; dependency of restartable collision composition | 001 and phase/growth/checkpoint contracts | Standalone operators/deterministic checkpoint inspected | Specify split order, saved RNG/counters, commit and rollback semantics |
| 011 electron transport/energy | Scope proposal; root selects implementation agent | Term-specific energy derivation and integrator contract | Ideal-only evolving-pressure restriction verified | Derive matched resistive heating/loss discrete balance |

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

- [ ] Cell-local self/cross pairing, empty/singleton/odd populations, arbitrary
  species counts/weights and validated staging before committing velocities.
- [ ] Independent Landau-generator drift/diffusion and equilibration-rate tests,
  density/charge/mass/timestep scaling, no cross-cell interaction, alias/error
  rollback, fixed-RNG replay and legacy differential checks.
- [ ] Reused-workspace allocation and complete-call cost matrix; cells with extreme
  occupancy and species imbalance must be measured, not excluded from pairing.
- [ ] CAP-010 owns automatic splitting and saved RNG/phase. CAP-003 must supply
  cell-complete local populations after migration; CAP-004 supplies real device
  pairing/streams. The local CPU slice alone closes none of these integration gates.

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

- [ ] **002A: Maxwell vacuum/cavity.** Dimension-parametric staggered state,
  periodic cross-check and conducting walls; explicit centering/CFL contract;
  compatible curl/divergence and field-energy/Poynting diagnostics. Gate cavity
  frequencies, analytic pulse transport, reflection phase and spatial/time order
  at multiple modes, incidence directions and corner configurations.
- [ ] **002B: bounded particle coupling.** Deposit the complete trajectory before
  removal/injection; split multi-cell and multi-face crossings, preserve IDs and
  species charge, account for wall surface charge/current. Gate discrete local
  continuity including boundary flux, Gauss propagation, particle-field work,
  reflecting-wall energy, exact free-flight crossing maps and empty populations.
- [ ] **002C: absorbing/radiating boundaries.** Add a measured absorbing layer
  with configurable thickness/profile and accounted energy removal; inject
  external fields/particles with explicit energy/charge budgets. Gate reflection
  spectra versus angle/frequency/resolution and long-time stability, including
  face/edge/corner layers and particle exits through them.
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

**Existing scope matters.** The analytic CN multiplier is a linear reference;
the shock option is a finite fixed-point approximation with no solve residual
check. Neither implements an implicit coupled kinetic-ion/fluid-electron step.
The [Stanier–Chacón–Chen method](https://arxiv.org/abs/1803.07158) provides a
primary-source design reference for compatible cell-centered spatial operators,
implicit midpoint, ion subcycling/orbit averages and conservation. Its properties
must be re-derived for the selected implementation rather than inherited by name.

- [ ] **005A: audit existing implicit-labelled paths.** Measure the four-iterate
  shock solve's residual, amplification and timestep range; compare a small
  assembled linear system solve. Decide whether to implement a converged solve
  or explicitly name/bound the approximate scheme. Align `recommended_dt` and
  documentation with the actual method. Preserve working resolved-step behavior.
- [ ] **005B: residual and transaction.** Define the nonlinear unknowns, discrete
  equations, normalization, residual norms, absolute/relative tolerances,
  iteration limits and candidate-state ownership. Failed convergence returns a
  diagnostic error with no physical state/clock commit; any timestep retry is
  explicit and bounded, not silent acceptance of the last iterate.
- [ ] **005C: coupled hybrid.** Start periodic scalar-electron hybrid with one
  population, then collections, evolving pressure, orbit averaging and selected
  dissipative terms. Put residual/linearization/preconditioner in integrators and
  reuse force/moment owners. Derive adjoint work and momentum cancellations for
  the actual spatial/particle discretization before claiming exact conservation.
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
  outputs. Current builders clone moving HEADs (`build_nhds.sh:18`,
  `build_vpic.sh:27`). Missing/build-failed references stay failed/unavailable.
- [ ] **007B: production dispersion.** Keep the independent eigenvalue comparison
  in case 27, but also measure frequencies/damping/growth from actual initialized
  hybrid/Hall/EM time series and compare with the applicable independent theory
  or solver. Sweep modes, angles, temperature/mass ratios and refinement; freeze
  fitting windows, uncertainty estimation and branch selection before results.
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

Review validation prose too: case 27's README says missing builds skip, whereas
its current runner:16–18,43–55 returns a failure metric. Keep historical measured
claims only with recoverable commands, source hashes and outputs; current-run
source inspection does not revalidate old numbers.

## CAP-008 — arbitrary kinetic species in ES/EM

This extends the documented current species restriction; it is not a defect in
that restriction. Generalize the canonical ES/EM charge/current and workspace
owners to ordered populations with arbitrary finite charge, positive mass and
independent storage. Preserve scalar/two-population constructors through the
same implementation, with explicit background/neutrality policy and per-species
time levels. Never normalize away nonzero physical net charge silently.

- [ ] ES collections, shared Poisson source and charge/momentum/energy diagnostics.
- [ ] EM collections, charge-conserving summed currents, mobile/immobile species,
  relativistic restrictions and explicit multirate schedule.
- [ ] Collection checkpoint/archive, CAP-001 interactions, CAP-002 boundaries,
  CAP-003 partitions and CAP-004 residency.

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

- [ ] Derive and implement Hall/pressure/inertia terms supported by the declared
  model, compatible mixed operators and accurate timestep/implicit solve control.
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

Add an explicitly owned evolution composition in the integrator layer; reuse
collision, ionization, radiation and antenna owners rather than embedding copies
in each model. Define the order and physical time of each source operation,
charge/energy/momentum exchange, particle growth, ID allocation and recomputation
of fields. First implement the accepted CPU collision/field split, then add other
sources one by one with their own conservation/accuracy contracts.

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

## Additional risks requiring tests, not assertions

| Candidate | Verified source observation | Exact settling work |
| --- | --- | --- |
| Converged shock field solve may still leave coupled energy/accuracy limits | The fixed-iteration residual defect is now reproduced on both runtimes; CAP-005A owns its repair. The reduced SBP/SAT and particle splitting remain separate contracts | Verify the accepted residual solve, then assess whole-step temporal convergence, boundary work, diffusion energy and particle/field exchange; a converged linear solve alone is not a global energy-stability proof |
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
