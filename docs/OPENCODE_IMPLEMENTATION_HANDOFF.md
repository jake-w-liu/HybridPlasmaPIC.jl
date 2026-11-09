# OpenCode implementation and acceptance handoff

Prepared 2026-11-09. **Codex stops after this handoff; OpenCode has not started.** The full simulator goal remains incomplete. This plan closes the versioned CAP-001–011 scope, all 79 existing milestones, and the compatibility, audit and resource requirements below. It cannot certify every conceivable plasma model or worldwide superiority.

Use the [79-row acceptance tracker](OPENCODE_ACCEPTANCE_TRACKER.md), [44 work-package records](../handoff/2026-11-09/WORK_PACKAGES.json), [priority issue register](../handoff/2026-11-09/ISSUE_REGISTER.json) and [portable packet index](../handoff/2026-11-09/README.md). The detailed source reviews supplement each package: [CAP-002–006](../handoff/2026-11-09/reviews/CAP002_006_REVIEW.md), [CAP-001 and CAP-007–011](../handoff/2026-11-09/reviews/CAP001_007_011_REVIEW.md), and [remaining source/transport slices](../handoff/2026-11-09/flux-proposal/REMAINING_CAP010_CAP011.md). These are frozen evidence, not proof that the proposed methods passed execution.

## Exact starting state

| Item | Checked state and consequence |
|---|---|
| Accepted-main snapshot | `9b4d363806722319c2241e4af04ed1ed1da00119`; [run 41813949753](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41813949753) passed all 26 jobs. This is historical CI acceptance, with science/resource work still open. Re-fetch before promotion. |
| New source baseline | `5aa628513fb854c96b2b387780edef1ea1faa66b`, branch `implement/cgl-integration-20261109`. The handoff commit adds documentation/packet files; its `src`, `ext`, `test`, `Project.toml` and workflow tree must match [CODE_TREE_BASELINE.txt](../handoff/2026-11-09/CODE_TREE_BASELINE.txt). |
| Current loader failure | **VERIFIED:** [job 129005385324](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42570140304/job/129005385324), raw log lines 765–771, reports a Julia `ParseError` at `src/models/empic_collection.jl:1430`. The same invalid form remains at line 1854. The [retained MPI artifact](../handoff/2026-11-09/ci/MPI_FAILURE_ARTIFACT.json) independently reports it. Fix P01 before accepting any new suite. |
| Required parser repair | Review both `es isa EMPIC ? :vector :planar` expressions; the explicit ternary separator is `es isa EMPIC ? :vector : :planar`. Search analogous sites and run minimum/latest remote package-load, formatter, host and MPI jobs. This handoff does not apply the repair. |
| Additional compiler failures | Whole-tree JET cap stays 86. Earlier `5fec151` produced 94/87, `7ea6ae8` 96/88; full old chains and additions were source-traced. `566ddc1` repairs two owners but its execution is pending. New `5703c95` latest reports 91 and `bc86550` minimum 99; their full added-tree disposition remains open. Package-located zero does not pass the whole-tree cap. |
| Other CI work | Formatting failed on several source candidates. `75f74bd` adopts a reviewed runner patch; later source needs fresh formatting after parser repair. `123c4a2` restores the missing pressure-essential IO prerequisite; rerun it. [CI snapshot](../handoff/2026-11-09/ci/STATUS_BEFORE_HANDOFF_PUSH.json) and [failure summaries](../handoff/2026-11-09/ci/LATEST_FAILURE_LOG_SUMMARIES.json) are dated observations, not current green status. No inspected evidence establishes a billing exception. |
| Implemented, unaccepted additions | Physical-pencil controls, CUDA first solver, toroidal annulus, implicit evolving pressure, grounded/periodic multidimensional AMR, persistent schedules, field-aligned conduction, NGP/CIC/TSC thermal exchange and ideal CGL composition are integrated as source candidates. Older focused successes do not accept the latest combined tree. |
| CGL at `5aa6285` | Ideal algebraic global-invariant CGL with ion-ion collisions; original Float16/Float32/Float64/BigFloat invariants, full anisotropic rollback and schema-4 persistence. Root reviewed the complete changed path/callers and 238-line essential file. [26 source/text checks](../handoff/2026-11-09/cgl-reviewed/ROOT_SOURCE_CHECKS.json) pass; all new executable/science/resource assertions remain unrun. Dissipative/evolving tensor CGL and scalar electron-ion thermal exchange are outside this slice. |
| Coupled controller | Complete nine-path [frozen patch](../handoff/2026-11-09/controller/coupled-empic-final-v1.patch), base `566ddc1`; **not adopted, not runtime accepted**. It applies mechanically to the source baseline, but retains the parser defect at frozen lines 1480/1904. Reconcile the parser fix during P30 review. Final 1,099-line model and other final files need complete integrator review; an earlier partial read does not approve them. |
| Flux-limited conduction | [Frozen read-only proposal](../handoff/2026-11-09/flux-proposal/PRE_EDIT_CONTRACT.md), **design choices not approved and no implementation**. P40 must settle the cold-front/mixed-zero domain, geometry and entropy/rounding decisions before code. Do not interpret this handoff as approval of that proposal. |
| Dependency release | SpectralOperators 0.1.9 registration, annotated tag and [release](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.9) are complete. HPIC requires it. Re-resolve and verify compatibility on actual minimum/latest Julia; no registration work remains from the earlier 0.1.4 request. |
| Progress count | 23 legacy checked / 56 open / 79 = **29.1% of unweighted recorded milestones**. All 11 CAP areas retain work. No reliable overall simulator-completeness percentage is established. |

## Takeover constraints

Work in an isolated checkout of the published handoff branch. Preserve every existing dirty checkout and inherited process. The original primary is `../HybridPlasmaPIC.jl` at `563888f`; its seven protected file hashes and the compiler ledger hash are in [MANIFEST.json](../handoff/2026-11-09/MANIFEST.json). Review that dirty radiation/forcing work independently before adoption; do not overwrite, reset, stash, delete or silently promote it. Preserve the 34-entry compiler ledger byte for byte. Do not terminate an inherited process or CI run; stop only a process proved launched for the active task in its own captured session.

Keep the shared PC free of Julia workloads. No local Julia, parser, compiler, formatter, test, numerical proxy, profiler, MPI or device computation. Read source locally; run required verification on GitHub or approved remote hardware. Lightweight Git, text, dependency-graph and checksum checks are allowed. Source-only review cannot establish runtime or scientific correctness.

Push each coherent fix/implementation to a validation branch, inspect exact-head failures, repair confirmed causes and re-run affected checks. Main receives the reviewed exact revision only after its full required CI passes. This validation-before-main order protects accepted main while satisfying the user's main-publication goal. A later merge/rebase changing source needs fresh exact-head acceptance. Never raise the JET cap, weaken filters, edit the historical ledger, loosen numerical tolerances, remove cases, skip failed jobs or convert errors to success to obtain green CI. A billing exception requires actual recorded GitHub billing evidence; it does not close unrun science/hardware gates.

## First actions and parallel execution

1. **P00:** read this plan and packet index, verify hashes, fetch refs, record dirty worktree state without mutation, and establish the support/evidence matrix. Assign one integration lead and exclusive exact-path owners. Preserve the baseline and every rejected candidate.
2. **P01 first:** repair both loader expressions; inspect the raw JET trees and all actual failed jobs. Run remote load/format preparation first, review the generated formatter patch, then run focused and full CI. Source work on independent paths may proceed while CI runs, but no package is accepted from the blocked baseline.
3. Freeze one small pre-edit contract per coherent slice: equations, admitted inputs, canonical owners/callers, phase/commit boundary, IO schema, independent oracle, resource targets and relevant matrix rows. Set error/order/noise thresholds from derivation or a pinned reference before seeing results.
4. Implement the complete chosen slice (evolution, diagnostics, rejection, persistence and docs), perform full caller/source review, then submit compact essential controls remotely. Avoid writing a new testing framework or broad abstractions. Never postpone a known failure to an unrelated final testing phase.
5. Reconcile each frozen worker patch against current owners, read the whole final result, inspect exact-head CI and add raw evidence. Update the original roadmap, crosswalk, package statuses and user-facing docs together. Continue another independent ready package while remote jobs execute; report status at least every minute during active Codex work.
6. Complete the applicable dependency joins, scientific campaigns and measured optimization. Run a fresh independent whole-code audit after the final implementation, then qualify the final integrated revision and publish it to main under P43.

Acceptance dependencies below govern when a package can close, rather than forcing unrelated design/source work into a serial queue. Use a ready queue; do not block all implementation on one hardware campaign. Split a blocked composite package into owned sub-slices with their own evidence, leaving the parent open. Every composed row additionally depends on its actual serial physics, geometry, schedule, IO and backend prerequisites even when the high-level graph omits that particular row.

| Parallel lane | Packages and initial boundary | Shared paths requiring lead coordination |
|---|---|---|
| CI/integration | P00–P01; central registrations and final promotion | `src/HybridPlasmaPIC.jl`, `Project.toml`, `.github/workflows/`, shared IO unions/codecs, global docs, `test/runtests.jl` |
| Collisions and sources | P02–P03, P35–P39 | `collisional_hybrid.jl`, collision/checkpoint/RNG owners and particle growth; serialize overlap with transport and schedules |
| Bounded EM/open shocks | P04–P07, P32–P34 | Yee field/boundary/coupling owners; coordinate implicit EM and source/current ownership |
| Distributed and device | P08–P15 | Keep pencil/CUDA/Metal paths separate where possible; sole writer for extension registration and shared array/operator interfaces |
| Implicit, geometry and AMR | P16–P24 | Separate family files; shared nonlinear/pressure/current/mesh owners require explicit transfer |
| Species and control | P29–P31 | Collection/full-PIC/kinematics/checkpoint paths; no concurrent controller and parser or source-calendar edits |
| Electron transport | P40–P41 | Energy-equation, closures, collision/integrator/checkpoint owners; coordinate source lane |
| Physical validation | P25–P28 alongside all lanes | Independent references/data may proceed without production edits; coordinate shared examples and case labels |
| Final audit/release | P42–P43 after dependency joins | No concurrent source changes during final acceptance snapshot |

An agent receives a sizeable independent slice, exact allowed paths and a handoff boundary. Never assign two writers to the same actual path across worktrees. Integration lead owns imports/exports/registrations/docs and merges; agents propose hunks for shared paths. Stop or finish overlapping writers before applying formatting or patches. A frozen handoff includes complete files/patch, base/head/hash, all callers, formulas, schema, source review and pending gates. All Codex workers are now finished; no ownership is implicitly active.

## Acceptance contract for every package

Record G1–G10 separately as pass, fail, pending or physically justified N/A. A source-only change may be `source reviewed`; it is `complete` only after every applicable gate passes. An unavailable required device, dataset or independent benchmark remains pending, rather than N/A. Do not demand nonsensical Cartesian products (for example adiabatic-electron GK with a fluid pressure module); derive applicability from the model and document each exclusion. Missing promised physics is never an applicability exclusion.

| Gate | Evidence required to pass |
|---|---|
| G1 — model | Actual equations, normalization/units, centering, geometry, boundary/source exchanges, admitted dimensions/shapes/precision and limitations; derive conservation as exact-discrete, convergent or statistical as appropriate. A defect ledger cannot replace physical accuracy. |
| G2 — architecture | Read canonical owners and every affected caller, constructor, diagnostic, optional extension and checkpoint reader. One physical policy owner where equations agree; preserve public API/schema or document a real migration. No fake fallback, silent host solve or unimplemented mode advertised as supported. |
| G3 — baseline and discriminators | Freeze parent source, resolved dependencies, input/RNG, environment and pre-edit contract. Retain rejected candidates. Exercise at least three relevant wrong-method challenges (for example sign error, omitted phase/source term, wrong metric/weight) using compact controls that would fail those methods. Use existing essential harnesses; do not build a mutation framework. |
| G4 — independent numerical/scientific result | Use derived/manufactured/assembled independent oracles and pinned independent data/code for production claims. Separate timestep, mesh and particle-number/noise refinement. Declare expected order, fit range, noise/uncertainty, tolerance and held-out observables before results; examine where order is valid and breaks down. Actual fields and particles must evolve self-consistently for solver claims. |
| G5 — failure/publication | Enumerate aliases, invalid/nonfinite/range inputs, empty/boundary cases, budget exhaustion, incomplete solves and IO failure. Define the commit point. Before it, failed operations preserve physical state, time, phase, persistent reports/history, IDs and required RNG. Distinguish recoverable collective validation from fatal transport failure; never promise MPI recovery without implementing it. |
| G6 — persistence | Cover raw/archive forms, versions/migrations, corruption/stale compatibility, independent destination objects and fresh-process continuation. Retain every physical half-step/controller/source/layer/surface/RNG/mesh state. Require same-environment deterministic behavior only where implemented; changed rank/device environments need declared numerical/statistical equivalence. |
| G7 — language/compiler/CI | Actual Project minimum Julia and floating latest stable, all focused/neighbors, full tests, Aqua, concrete inference and complete whole-package JET leaf+frame trees, formatting and all required exact-head workflows. Record actual Julia and dependency versions from logs. Existing cap/filter/ledger and numerical/resource gates stay intact. No green primitive run substitutes for composed-family CI. |
| G8 — actual backends | Execute the declared path on real 2/4/8 ranks, the relevant CUDA and Metal hardware, and multiple nodes/devices for those claims. Cover ownership, migration, halos, reduction/FFT normalization, device failures and restart. Emulation, hosted extension admission and one-rank calls do not establish device/distributed acceptance. |
| G9 — resources and independent audit | Predeclare paired baseline/candidate workloads and accuracy targets; measure whole initialization, evolution, diagnostics, save/restore, allocation, peak memory and transport at required problem sizes/edge regimes. Meet existing absolute and relative gates. Derive new budgets from a recorded remote baseline before editing; do not invent ceilings afterward. Re-profile hotspots and have a fresh reviewer trace final formulas/callers and rejected inputs. |
| G10 — user-facing truth/release | Update README, guides, examples, API/schema, support matrix, roadmap and acceptance records at the same source revision. Mark limitations/pending backends explicitly. Publish only the exact accepted source and retain raw logs/data plus reviewable links. |

For unchanged ordinary paths preserve documented bits, zero-allocation promises and fixed absolute budgets. A paired slowdown screen of 1.5× can flag regressions where already applicable; it is not an optimization target or permission to relax a stricter existing gate. New hotspot targets must be predeclared from measured baseline evidence and accuracy requirements. Optimizing only a primitive while full-step time/memory regresses does not pass G9. Explain every retained necessary cost. No finite audit proves global bug absence or a global speed maximum.

## Support matrix and evidence records

Create versioned rows for actual meaningful combinations: model/family; geometry/dimension; particle shape and kinematics; scalar/evolving/CGL/electron-inertia closure; collision/source law; species/calendar/controller; field/material/open boundary; serial/pencil/CUDA/Metal/multi-device backend; precision; raw/archive schema. Include new/empty/unequal species, odd grids, seams, event/regrid/migration phases and variable steps. Every promised applicable combination must map to a package and case; separate unsupported/unimplemented from validated. The broad Cartesian product is not automatically promised, but the versioned required scope cannot shrink to make the tracker green.

Each evidence row stores:

- Stable milestone and package IDs, owner, source/API/physics domain and applicable gates.
- Parent, candidate and accepted commit SHA; full command/workflow input; resolved dependency manifest, actual Julia/runtime, threads/ranks/device/node/driver and input/seed hashes.
- Predeclared independent oracle, expected order/error/noise and resource budgets; baseline/candidate raw outputs and rejected candidates.
- Pass/fail/pending per gate, raw artifact URL plus checksum, reviewer disposition and unresolved limitations. Record code/source review separately from executable/science/resource acceptance.
- For each CI failure: exact job, first actual failing operation and raw error; checked source owner, fix, affected callers and rerun. Causal statements require this evidence. Unmeasured explanations are labeled hypothesis with the decisive test.

Re-fetch actual latest stable Julia when work starts and at final release; the matrix's `1` selector must resolve to it and logs must show the actual version. The previously observed 1.12.7/1.13.3 versions are evidence snapshots, not a permanent definition of latest. Verify dependency compat/registration on the resolved environments, not just Project text.

<a id="packages"></a>
## Detailed work packages

Every package below inherits G1–G10 and the precise milestone text in the acceptance tracker. The canonical paths identify current owners to read; new files are allowed only after the architecture review. Scope includes the detailed derivation/acceptance tables in the frozen appendices. Source present does not imply fully implemented or accepted.

<a id="p00"></a>
### P00 — Freeze and reconcile takeover

**Scope:** CAP-GLOBAL · **Acceptance dependencies:** none (startup).

**Current state:** Handoff prepared; operator startup and dirty-primary review remain. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `Project.toml`, `docs/DEEP_DEBUG_TRACKER.md`.

**Deliver:** Inventory refs/worktrees/source packets, preserve the seven primary hashes and compiler ledger; establish explicit model/support/evidence matrices.

**Close only when:** A clean isolated operator checkout, verified artifact hashes, no lost dirty work, exact source/dependency/hardware baseline and named owners; every required combination has a row.

<a id="p01"></a>
### P01 — Close current CI and compiler regressions

**Scope:** CAP-GLOBAL · **Acceptance dependencies:** [P00](#p00).

**Current state:** VERIFIED package-load ParseError remains at two EMPIC ternaries; new 91/99-report JET failures need complete disposition; prior repairs are not current-head acceptance. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `.github/workflows/HybridPlasmaPIC-CI.yml`, `.github/workflows/implementation-essentials.yml`, `test/test_jet.jl`, `src/models/empic_collection.jl`.

**Deliver:** First repair both invalid symbol-valued ternaries at src/models/empic_collection.jl:1430/1854 and reconcile the frozen controller. Then triage every actual failed job, including new 91/99-report JET trees, formatter artifacts, focused prerequisites and real MPI failures. Preserve all numerical/compiler/resource gates.

**Close only when:** Every failure has raw-log/source disposition; actual minimum/latest Julia, focused/neighbors/Aqua/concrete and whole-tree JET plus all required exact-head jobs pass. Actual billing evidence is the only CI exception; unexecuted science/hardware stays open.

<a id="p02"></a>
### P02 — Complete serial collision resources and regimes

**Scope:** CAP-001 · **Acceptance dependencies:** [P01](#p01).

**Current state:** CPU pairing/scattering exists; broader cost/range campaign open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/particles/coulomb_local.jl`, `src/particles/collisions.jl`, `src/particles/coulomb_scattering.jl`.

**Deliver:** Retain canonical self/cross/reuse physics; complete occupancy/species/weight/exponent-span and full-operation optimization campaigns.

**Close only when:** Independent drift/diffusion/rate/isolation and seeded unequal-weight statistics; empty/singleton/odd cells; paired init/bin/shuffle/scatter/diagnostic/restore time, allocation and peak-memory targets frozen before optimization.

<a id="p03"></a>
### P03 — Complete rank/device collision ownership

**Scope:** CAP-001 · **Acceptance dependencies:** [P02](#p02), [P08](#p08), [P13](#p13).

**Current state:** Pencil ES/EM source exists; broader/device composition missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/pencil_collisional_pic.jl`, `ext/pencil_collisional_empic.jl`, `ext/HybridPlasmaPICCUDAExt.jl`.

**Deliver:** Keep complete physical collision cells together or derive cross-rank pairing; implement device pairing/stream/reuse ownership and necessary solver compositions.

**Close only when:** Actual 2/4/8 ranks and real devices: unique pair coverage, empty peers, imbalance, rates/conservation-in-expectation, rank-local/device rejection and full advanced-stream continuation; no host/emulation substitution.

<a id="p04"></a>
### P04 — Accept current bounded Yee families

**Scope:** CAP-002 · **Acceptance dependencies:** [P01](#p01).

**Current state:** Bounded vacuum and restricted particle maps exist. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/yee_electromagnetic_pic.jl`, `src/models/yee_impact_longitudinal_pic.jl`, `src/models/yee_absorbing_longitudinal_pic.jl`.

**Deliver:** Validate existing finite-cloud PEC, event-free EM and restricted longitudinal impact/capture maps before extending them.

**Close only when:** Independent Yee/cavity and trajectory/source oracles; continuity/Gauss, actual particle/grid/plate work, both polarizations, force-driven refinement, late undo and full phase raw/archive continuation. Ledger closure alone does not pass physical accuracy.

<a id="p05"></a>
### P05 — Implement general material particle boundaries

**Scope:** CAP-002 · **Acceptance dependencies:** [P04](#p04).

**Current state:** D2/3/general repeated events remain missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/yee_electromagnetic_pic.jl`, `src/models/yee_impact_longitudinal_pic.jl`, `src/models/yee_absorbing_longitudinal_pic.jl`.

**Deliver:** Add complete D2/3 crossed-cell/face/edge/corner paths, repeated impact/capture/injection, transverse/material/sheath/emission laws and surface reservoirs/IDs.

**Close only when:** No lost path remainder or duplicated charge/IDs; independent local continuity including wall flux, impulse/energy receiving reservoirs, simultaneous-event termination/refinement, late failure/retry and complete pending-event restart.

<a id="p06"></a>
### P06 — Accept absorbers and implement open plasma coupling

**Scope:** CAP-002 · **Acceptance dependencies:** [P04](#p04), [P05](#p05).

**Current state:** Vacuum conductivity/ADE PML source exists; coupled/open scope incomplete. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/yee_pml.jl`, `src/models/yee_conductive_layer.jl`, `src/io/yee_pml_checkpoint.jl`.

**Deliver:** Measure existing layers; derive charge-compatible PIC/source exits and injection; add general CPML parameters/conditions only with a concrete operator contract.

**Close only when:** Reflection versus angle/frequency/polarization/thickness/resolution, corner overlap and late-time behavior; actual layer memories/work and restart; plasma charge/Gauss/boundary flux and long-run stability tested together.

<a id="p07"></a>
### P07 — Compose nonperiodic production physics

**Scope:** CAP-002 · **Acceptance dependencies:** [P05](#p05), [P06](#p06), [P10](#p10), [P14](#p14), [P31](#p31).

**Current state:** Broader relativistic/scheduled/backend boundaries missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/yee_electromagnetic_pic.jl`, `src/models/full_pic.jl`, `ext/pencil_empic.jl`.

**Deliver:** Compose relativistic and scheduled species, curved material fields/thermal boundaries, domain ownership and residency with the accepted boundary/event model.

**Close only when:** Each supported combination has phase/current/impulse/work and full persistent state; independent bounded-plasma comparison, whole-operation rejection and actual-rank/device acceptance at matched error.

<a id="p08"></a>
### P08 — Close distributed layout and operator contract

**Scope:** CAP-003 · **Acceptance dependencies:** [P01](#p01).

**Current state:** Periodic operators exist; broader layout/conditioning scope open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/parallel/pencil_decomposition.jl`, `ext/pencil_fourier.jl`, `ext/pencil_particles.jl`.

**Deliver:** Qualify owned/permuted/thin/empty layouts and all directions; implement required nonperiodic/wider-shape/staggered support without a global hot-path gather.

**Close only when:** Independent local/global operator/moment oracles at actual ranks, original-factor/Fourier null/DC/Nyquist checks, collective invalid-input rejection, bounded owned-plus-support memory and complete communicator/plan cleanup.

<a id="p09"></a>
### P09 — Accept actual-pencil solver lifecycles

**Scope:** CAP-003 · **Acceptance dependencies:** [P08](#p08).

**Current state:** Charged ES/EM/hybrid/Maxwell and restart source exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/pencil_electrostatic.jl`, `ext/pencil_empic.jl`, `ext/pencil_hybrid.jl`, `tools/ci_distributed_fft.jl`.

**Deliver:** Complete numerical, collective, physical-phase and matching-layout persistence acceptance for every existing actual-pencil family.

**Close only when:** Nonuniform charged multiple species, real domain/seam migration and empty peers; serial observable comparisons, full saved arrays and variable-step phase; asymmetric rejection/I/O finishes collectively without partial publication or hang.

<a id="p10"></a>
### P10 — Implement broader distributed families

**Scope:** CAP-003 · **Acceptance dependencies:** [P09](#p09), [P35](#p35).

**Current state:** Broader closures/shapes/relativity/open/source families incomplete. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/pencil_hybrid.jl`, `ext/pencil_empic.jl`, `ext/pencil_collisional_pic.jl`.

**Deliver:** Extend admitted pressure/transport/CGL, shapes/schedules/relativity, source and mixed/open families using accepted serial physics and actual owned operators.

**Close only when:** Summation before floors, source-range recovery, charge/current/energy and real migration/rollback/full restart for each matrix row; reject unimplemented ownership cases before collective work.

<a id="p11"></a>
### P11 — Implement repartition and cluster scaling

**Scope:** CAP-003 · **Acceptance dependencies:** [P09](#p09), [P10](#p10).

**Current state:** Changed layout/rank and real multi-node closure open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/pencil_checkpoint.jl`, `ext/pencil_particles.jl`, `src/parallel/pencil_decomposition.jl`.

**Deliver:** Add rebalance, changed-rank/layout physical and RNG redistribution and generation/recovery policy; run real multi-node scaling.

**Close only when:** IDs/weights/charge/streams neither lost nor duplicated; same-layout deterministic and changed-layout statistical contracts separated; strong/weak scaling, slowest rank, communication, imbalance and peak storage on named nodes.

<a id="p12"></a>
### P12 — Close actual GPU primitive contracts

**Scope:** CAP-004 · **Acceptance dependencies:** [P01](#p01).

**Current state:** CUDA/Metal primitives exist; actual device/range gates open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/HybridPlasmaPICCUDAExt.jl`, `ext/HybridPlasmaPICMetalExt.jl`.

**Deliver:** Verify native lowering, arithmetic/precision/mode/range, gather/deposit/push/wrap/reduction and failure/stream-context contracts on each real device.

**Close only when:** Scalar indexing disabled; independent primitive/trajectory/source oracles including empty/block edges, alias/range failures and subsequent reuse; actual synchronized kernels, transfers, allocations and context restoration.

<a id="p13"></a>
### P13 — Accept and extend complete resident first solvers

**Scope:** CAP-004 · **Acceptance dependencies:** [P12](#p12).

**Current state:** CUDA D1/F64/CIC/isothermal lifecycle exists; real acceptance and other variants open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/cuda_hybrid.jl`, `ext/cuda_hybrid/types.jl`, `ext/cuda_hybrid/evolution.jl`, `ext/cuda_hybrid/checkpoint.jl`.

**Deliver:** Accept the current CUDA family, then implement D2/3/F32/other shapes/variable intervals/range scope and a genuine Metal resident field-plus-particle lifecycle.

**Close only when:** Actual-device init/prime/repeated coupled steps, independent dispersion/constraints/work, late async rejection and changed-count raw/archive continuation; measured absence of per-step full field/particle host solves/copies.

<a id="p14"></a>
### P14 — Implement resident physics families

**Scope:** CAP-004 · **Acceptance dependencies:** [P13](#p13), [P35](#p35).

**Current state:** Broader closures, solvers, sources and geometry/device combinations missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/cuda_hybrid.jl`, `ext/cuda_hybrid/`, `ext/HybridPlasmaPICMetalExt.jl`.

**Deliver:** Port each required accepted serial family and closure/source/phase to resident storage/operators, including ES/EM, transport/tensors, collisions, implicit/geometry/GK/AMR where valid.

**Close only when:** Each advertised device row runs complete numerical/scientific/failure/restart and resource gates; atomic-order tolerance is declared, host replay or import-only success cannot close it.

<a id="p15"></a>
### P15 — Implement and accept multi-GPU

**Scope:** CAP-004 · **Acceptance dependencies:** [P11](#p11), [P14](#p14).

**Current state:** Missing complete multi-device composition. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/parallel/pencil_decomposition.jl`, `ext/HybridPlasmaPICPencilFFTSExt.jl`, `ext/HybridPlasmaPICCUDAExt.jl`.

**Deliver:** Compose domain-owned arrays/migration with direct-device or explicitly accounted host staging, ordering/streams, local RNG and repartition/restart.

**Close only when:** Actual multiple devices/nodes, verified transport path, unequal occupancy, arbitrary migration, late failure/stream separation, supported restart and matched-error scaling/memory/transfer evidence.

<a id="p16"></a>
### P16 — Close implicit solver audit and common numerical owner

**Scope:** CAP-005 · **Acceptance dependencies:** [P01](#p01).

**Current state:** Residual-controlled sources exist; advice/preconditioning/cost gates open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/integrators/implicit_pic_nonlinear.jl`, `src/integrators/gmres_storage.jl`, `src/initial_conditions/shock_cn.jl`, `src/integrators/timestep_bounds.jl`.

**Deliver:** Qualify linear/reduced CN versus coupled models, true residual/scaling/Jv/GMRES/budgets; resolve timestep advisory contract and optimize measured repeated FFT/Krylov work.

**Close only when:** Independent assembled equations/actions/amplification, fresh final residuals, budget/domain/backtracking failures with reuse; no hidden explicit/smaller-step/direct fallback; matched-accuracy full-call residual/iteration/allocation cost.

<a id="p17"></a>
### P17 — Accept and extend coupled implicit hybrid

**Scope:** CAP-005 · **Acceptance dependencies:** [P16](#p16).

**Current state:** Isothermal and evolving-pressure joint source exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/integrators/implicit_hybrid.jl`, `src/electrons/energy_equation.jl`, `src/io/implicit_hybrid_checkpoint.jl`.

**Deliver:** Accept joint v/B/pe and old schema1/new schema2; implement derived orbit averaging/subcycled paths and method-specific preconditioning; compose later physics only after its contracts.

**Close only when:** Independent nonlinear pressure/velocity/B equations, true residual/work/positivity and zero-pressure/Jv limits; nontrivial order/dispersion/tolerance sweeps, late rollback, fresh/stale/legacy raw/archive continuation and time-to-error cost.

<a id="p18"></a>
### P18 — Accept and extend coupled implicit EM

**Scope:** CAP-005 · **Acceptance dependencies:** [P16](#p16), [P04](#p04).

**Current state:** Spectral and compatible periodic Yee maps exist; broader forms missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/integrators/implicit_em.jl`, `src/integrators/implicit_yee_em.jl`, `src/io/implicit_em_checkpoint.jl`.

**Deliver:** Accept current spectral/Yee source-current-field/work equations, then implement D2/3 compatible paths, wider admitted orbits/spacing and boundary/relativity/schedule/backend compositions.

**Close only when:** Independent coupled solves and all published equation/continuity/Gauss/work defects, transverse and multi-cell/seam controls, order separated from mesh/particle error, exact phase restart and late failure. Tighter nonlinear tolerance cannot erase spatial correction work.

<a id="p19"></a>
### P19 — Accept coupled shear and toroidal annulus

**Scope:** CAP-006A · **Acceptance dependencies:** [P01](#p01).

**Current state:** Genuine shear/annulus coupled source exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/sheared_periodic_hybrid_pic.jl`, `src/models/toroidal_annulus_hybrid_pic.jl`, `src/fields/toroidal_annulus_hybrid_operators.jl`.

**Deliver:** Qualify coordinate/inverse/J/metric/source-work and all-stage ions/beta/pe/contact/reservoir evolution before broader geometry.

**Close only when:** Independent weighted adjoints, mapped uniform fields, analytic/manufactured geometry and coupled wave/equilibrium refinements; face inflow/outflow/SAT/work, normal trace/positivity, late rollback and all physical/history restart; physical boundary behavior measured beyond ledgers.

<a id="p20"></a>
### P20 — Implement full curvilinear/toroidal scope

**Scope:** CAP-006A · **Acceptance dependencies:** [P19](#p19), [P05](#p05).

**Current state:** Solid axis/material/equilibria/import/general maps and backends remain. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/meshes/toroidal_annulus_mesh.jl`, `src/coupling/toroidal_annulus_hybrid_coupling.jl`, `src/initial_conditions/toroidal_annulus_hybrid.jl`.

**Deliver:** Add solid-axis topology/parity/crossings, curved material walls, finite-temperature/flow equilibria and imported/shaped/general topology, followed by owned/device implementations.

**Close only when:** Axis/orbit/metric refinement and physical force balance, invariant and work checks appropriate to symmetry, actual coupled evolution and full boundary/equilibrium persistence; independent geometry reference and error-matched resource evidence.

<a id="p21"></a>
### P21 — Accept self-consistent full-f slab gyrokinetics

**Scope:** CAP-006B · **Acceptance dependencies:** [P01](#p01).

**Current state:** Full-f constant-B slab family exists; method/science/cost gates open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/slab_gyrokinetic_pic.jl`, `src/io/slab_gyrokinetic_checkpoint.jl`.

**Deliver:** Accept the declared ring-TSC/full-k polarization/modified adiabatic electron model through actual nonuniform collective RK evolution.

**Close only when:** Independent mode/FLR/electron-field residuals, genuine nonzero truncation order, separate time/mesh/marker/velocity/ring refinement, number/momentum/energy, later-stage solve failure/reuse and exact-head full continuation/resource gates.

<a id="p22"></a>
### P22 — Implement broader gyrokinetic physics and composition

**Scope:** CAP-006B · **Acceptance dependencies:** [P21](#p21), [P20](#p20), [P35](#p35).

**Current state:** Nonuniform/toroidal/delta-f/collision/kinetic-EM/backend stages missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/gyrokinetics.jl`, `src/models/slab_gyrokinetic_pic.jl`, `src/particles/collisions.jl`.

**Deliver:** Derive and implement geometry/mirror/bounce, selected evolving-weight/equilibrium model, GK collision Jacobian/invariants/RNG, kinetic electron/EM/free-energy and valid source/backend/AMR compositions.

**Close only when:** Independent mirror/bounce and matched Cyclone or appropriate code cases, FLR/free-energy and weight/equilibrium/RNG work, full solve-failure/restart and transport uncertainty; slab results cannot accept absent geometry.

<a id="p23"></a>
### P23 — Accept adaptive ES and implement remaining boundaries

**Scope:** CAP-006C · **Acceptance dependencies:** [P01](#p01).

**Current state:** Static/dynamic 1D and grounded/periodic D2/3 ES exist. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/meshes/cartesian_amr_hierarchy.jl`, `src/fields/cartesian_composite_poisson.jl`, `src/models/adaptive_composite_electrostatic.jl`.

**Deliver:** Accept actual hierarchy/cloud/Q1/KDK/regrid histories, then add neutral-Neumann/material/electrode/mixed-axis/harmonic-source contracts to the canonical assembly.

**Close only when:** Independent constraints/dense quadrature/manufactured fields/all-row gauge residuals, arbitrary-width and interface/self-force work, separate fixed-width/time/space refinements, repeated no-step/regrid/changed-topology continuation and late solve/path failure.

<a id="p24"></a>
### P24 — Implement adaptive hybrid/EM, subcycling and backends

**Scope:** CAP-006C · **Acceptance dependencies:** [P23](#p23), [P17](#p17), [P18](#p18), [P35](#p35).

**Current state:** Compatible current/flux synchronization and broader hierarchy missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/adaptive_composite_electrostatic.jl`, `src/meshes/cartesian_composite_adaptation.jl`, `src/fields/cartesian_composite_poisson.jl`.

**Deliver:** Derive compatible face/edge currents, divergence/Gauss/reflux and work; implement all coupled stages/local time interpolation, pending subcycles, partition/rebalance and device hierarchy.

**Close only when:** Independent coarse-fine current/flux/field identities, uniform-fine reference and transmission/refinement, real subcycle/regrid histories, full phase/RNG/topology restart, actual hardware and whole-step cost versus equal-error uniform fine.

<a id="p25"></a>
### P25 — Complete pinned references and production wave matrix

**Scope:** CAP-007 · **Acceptance dependencies:** [P01](#p01).

**Current state:** Reference infrastructure and parallel/oblique fluid histories exist. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `validation/mode_history.jl`, `validation/27_nhds_dispersion_comparison/run.jl`, `validation/34_production_parallel_dispersion/run.jl`, `validation/35_production_oblique_hall/run.jl`.

**Deliver:** Replay pinned references; extend production oblique kinetic/thermal/mass-ratio/ES/EM/damping regimes with correct branch/normalization and uncertainty.

**Close only when:** Theory never supplies fitted omega/rate; frozen windows/cadence/seeds and alias diagnostics, independent primary references, time/space/particle refinement and held-out observables. Preserve historical fitter failure and its later correction distinctly.

<a id="p26"></a>
### P26 — Complete nonlinear profiles and datasets

**Scope:** CAP-007 · **Acceptance dependencies:** [P25](#p25), [P32](#p32).

**Current state:** Existing B/scalar summaries do not cover all full profiles. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `validation/28_hybridvpic_perp_shock/run.jl`, `validation/21_published_preisser2020_summary/run.jl`.

**Deliver:** Retrieve licensed/checksummed matched data and code builds, then compare full n/B/u/p and shock-frame/time statistics rather than summary self-roundtrips.

**Close only when:** Matched physical equations/geometry/units/boundaries, source/build/data hashes, alignment/uncertainty and convergence; missing data/build remains open and is never replaced with fabricated traces.

<a id="p27"></a>
### P27 — Complete independent cross-family physics validation

**Scope:** CAP-007 · **Acceptance dependencies:** [P04](#p04), [P09](#p09), [P19](#p19), [P21](#p21), [P23](#p23), [P35](#p35).

**Current state:** Broad matched independent campaigns remain open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `validation/`, `src/verification/instability_sweep.jl`, `validation/32_magnetic_reconnection/run.jl`.

**Deliver:** Validate bounded EM, collisions, gyro/geometry/AMR and reconnection against appropriate independent models; reconcile labels with actual eta/inertia/source configuration.

**Close only when:** Converged held-out observables and conservation/Ohm/flux budgets with uncertainty; explicit parameter scans distinguish models. A code label or a closing scalar ledger does not establish the physical mechanism.

<a id="p28"></a>
### P28 — Establish comparative error-matched effectiveness

**Scope:** CAP-007 · **Acceptance dependencies:** [P25](#p25), [P26](#p26), [P27](#p27), [P02](#p02).

**Current state:** No universal leadership/optimization evidence exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `validation/`, `docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md`.

**Deliver:** Freeze problems, versions, named hardware, physical error/uncertainty and cost budgets; execute fair serial/rank/device comparisons including failed cases.

**Close only when:** Time-to-solution, compile/setup/steady work, memory, communication/transfer and scaling at matched error. Any superiority statement is limited to measured regimes and machines.

<a id="p29"></a>
### P29 — Accept arbitrary-species schedules and IO

**Scope:** CAP-008 · **Acceptance dependencies:** [P01](#p01).

**Current state:** Ordinary/dividing/persistent collection source exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/empic_collection.jl`, `src/models/full_pic.jl`, `src/io/empic_collection_checkpoint.jl`.

**Deliver:** Accept generic and specialized collections, fixed/mobile/tracer rows, original q/m and all distinct half/full phases, pending deadlines, current/work/controller histories and schema1/2/3.

**Close only when:** Independent corrected-current/endpoints and every subinterval continuity, species permutations and unequal q/m production physics, variable macro calls/pending right edges, late whole-macro undo, stale/aliased/corrupt raw/archive continuation and full CI.

<a id="p30"></a>
### P30 — Review/adopt coupled error controller; close remaining domains

**Scope:** CAP-008 · **Acceptance dependencies:** [P29](#p29).

**Current state:** Nine-file controller patch frozen and applicable; NOT adopted. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/empic_collection.jl`, `src/models/full_pic.jl`, `src/io/empic_collection_checkpoint.jl`.

**Deliver:** Read complete portable final9 and contract, register model/diagnostic/IO/573-line essential, verify old inverses, accept CIC/TSC full-time pair; separately derive NGP jumps and any true multirate-cost policy.

**Close only when:** Actual stage charge differs between low/high branches; independent force/current/Maxwell/clock/work, knot/seam/C1/domain/budget and step0-current controls; fresh-grid/raw/archive replay, order/effectivity/tolerance sweeps and complete-operation resources. No FFT portability or universal true-error theorem inferred.

<a id="p31"></a>
### P31 — Implement schedule/source/boundary/backend intersections

**Scope:** CAP-008 · **Acceptance dependencies:** [P30](#p30), [P35](#p35), [P05](#p05), [P09](#p09), [P14](#p14).

**Current state:** Current collision wrappers intentionally reject many combinations. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/empic_collection.jl`, `src/integrators/collisional_empic.jl`, `src/integrators/collisional_electrostatic.jl`.

**Deliver:** Derive each supported phase/current/source/impulse/work and transaction across collisions, relativity, boundary events, nonmobile rows, partition and residency.

**Close only when:** Noncommuting coupled physical controls and order, continuity/work for actual events/subintervals, no partial macro/RNG commit, full pending state and actual-rank/device restart; removing a guard alone never passes.

<a id="p32"></a>
### P32 — Qualify current full-Hall open shock family

**Scope:** CAP-009 · **Acceptance dependencies:** [P01](#p01).

**Current state:** D2/3 CIC isothermal numerical-reservoir source exists. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/open_hall_hybrid.jl`, `src/fields/mixed_shock_operators.jl`, `src/diagnostics/open_hall_energy.jl`.

**Deliver:** Accept actual mixed/SBP/SAT/projected particle-field-reservoir model and declare a physically qualified numerical-boundary regime.

**Close only when:** Independent transpose/face work/manufactured induction/divergence, whistlers/upstream/1D controls, SAT/distance/perturbation/reflection and mesh/time/particle/birth sweeps; full histories/RNG undo/restart. Finite projection or hnu admission is not boundary stability.

<a id="p33"></a>
### P33 — Implement open pressure/transport/inertia and material models

**Scope:** CAP-009 · **Acceptance dependencies:** [P32](#p32), [P40](#p40), [P41](#p41).

**Current state:** Broader open closures/thermal/material reservoir problems missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/open_hall_hybrid.jl`, `src/boundaries/open_hall_reservoir.jl`, `src/io/open_hall_checkpoint.jl`.

**Deliver:** Add actual stage pressure/enthalpy and compatible incoming/outgoing heat, then derived dissipative/inertial/material/sheath models and complete external work.

**Close only when:** Independent thermal/charge/momentum/energy fluxes, inlet/outlet/newborn phase, positivity without clipping, reflection/profile refinement and all reservoir/source/RNG/field state continuation and late failure.

<a id="p34"></a>
### P34 — Compose full-Hall shocks with distributed/device ownership

**Scope:** CAP-009 · **Acceptance dependencies:** [P33](#p33), [P10](#p10), [P14](#p14).

**Current state:** Broader rank/device open models missing. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/models/open_hall_hybrid.jl`, `src/boundaries/open_hall_reservoir.jl`, `ext/pencil_hybrid.jl`.

**Deliver:** Integrate admitted open reservoirs, growth/IDs, mixed operators and closure physics with real ownership/residency.

**Close only when:** Actual rank/device injection/migration/stream separation, physical flux/profile equivalence, collective/asynchronous late failure, full restart and scaling at matched shock accuracy.

<a id="p35"></a>
### P35 — Accept existing composed sources, thermal and CGL

**Scope:** CAP-010/011 · **Acceptance dependencies:** [P01](#p01), [P02](#p02).

**Current state:** Serial/pencil wrappers, NGP/CIC/TSC thermal and CGL source exist. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/integrators/collisional_hybrid.jl`, `src/electrons/thermal_exchange.jl`, `src/io/collisional_checkpoint.jl`, `src/electrons/field_aligned_conduction.jl`.

**Deliver:** Accept each precise new source/transport family at the integrated head; preserve scalar schema1, NGP2, overlap3 and CGL4 plus normalization and original coefficient policy.

**Close only when:** Independent physical relaxation/incidence/tensor/force/rate, actual native energy/momentum/entropy, sparse/cold/mixed-precision limits, true zero reduction, late RNG undo and full raw/archive continuation; no bath/drag/inertia/tensor scope inferred.

<a id="p36"></a>
### P36 — Implement reaction growth, namespaces and advanced source state

**Scope:** CAP-010 · **Acceptance dependencies:** [P35](#p35), [P29](#p29).

**Current state:** Primitives exist; complete reaction/neutral/growth composition open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/particles/collisions.jl`, `src/integrators/collisional_hybrid.jl`, `src/io/checkpoint.jl`.

**Deliver:** Compose ionization/neutral reservoirs, newborn phase/current, growth/storage and monotonic global/per-species IDs with full source clocks/counters/RNG.

**Close only when:** Independent threshold/rate/recoil/charge/energy statistics, IDs unique after deletion/growth/migration/restart, exhaustion/failure before partial commit and exact advanced-stream/ID/source continuation.

<a id="p37"></a>
### P37 — Implement variable/higher-order source evolution

**Scope:** CAP-010 · **Acceptance dependencies:** [P35](#p35), [P36](#p36), [P17](#p17).

**Current state:** Existing post-field stochastic split is first order. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/integrators/collisional_hybrid.jl`, `src/integrators/extrapolated_leapfrog.jl`, `src/io/collisional_checkpoint.jl`.

**Deliver:** Define stronger split/variable-step phases, accepted/rejected draws and all field/current/pressure refresh; reuse canonical source and field owners.

**Close only when:** Independent noncommuting source-field and weak/statistical-order refinement, rejected-step retry and uninterrupted/restarted histories with all clocks/RNG. A symmetric name or two half calls is not an order proof.

<a id="p38"></a>
### P38 — Implement accounted RF/ray/radiation/external sources

**Scope:** CAP-010 · **Acceptance dependencies:** [P36](#p36), [P35](#p35).

**Current state:** Standalone antenna/ray/reduced-radiation components exist. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/fields/antenna.jl`, `src/raycon/`, `src/particles/radiation_reaction.jl`.

**Deliver:** Reconcile preserved radiation edits, derive receiving reservoir/coupling time/current/force and external work, and own callback/resource commit behavior.

**Close only when:** Analytic/manufactured power/impulse/charge transfer, zero-source reduction, actual growth/phase if relevant, combined energy and full ray/source/external/RNG restart; rejected callbacks cannot hide partial external effects.

<a id="p39"></a>
### P39 — Implement complete rank/device sources and RNG growth

**Scope:** CAP-010/011 · **Acceptance dependencies:** [P03](#p03), [P10](#p10), [P14](#p14), [P36](#p36), [P37](#p37), [P38](#p38).

**Current state:** Full broader parallel/device source systems incomplete. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `ext/pencil_collisional_pic.jl`, `ext/cuda_hybrid/`, `src/particles/collisions.jl`.

**Deliver:** Compose source schedules/reservoirs/reactions, device-local RNG/ID/growth and collective publication with accepted physics and ownership.

**Close only when:** Actual devices/ranks with growth, empty/imbalanced partitions, stream uniqueness, one-rank/async failure, full-generation restart; deterministic same-layout versus changed-layout statistics, transfer/memory/cost explicitly measured.

<a id="p40"></a>
### P40 — Implement limited/state-dependent/anisotropic heat conduction

**Scope:** CAP-011 · **Acceptance dependencies:** [P01](#p01), [P35](#p35).

**Current state:** Constant scalar/tensor source exists; frozen limiter proposal unapproved. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/electrons/closures.jl`, `src/electrons/energy_equation.jl`, `src/electrons/field_aligned_conduction.jl`.

**Deliver:** Review frozen physical limiter proposal first; settle original-factor saturation, full Nyquist/classical limit, variable divergence-form mobility, actual prepared-flux cap, entropy and cold/boundary domains before coding.

**Close only when:** Independent actual-stage cap and source/mean/work/entropy defects, constant/zero/classical limits, varying-B/coefficient and finite saturated range rescue, positivity/late undo/full policy IO. Mixed-zero fronts need an admitted positive conservative model, not clipping or undefined-rate acceptance.

<a id="p41"></a>
### P41 — Implement full drag, electron inertia and evolving tensors

**Scope:** CAP-011 · **Acceptance dependencies:** [P35](#p35), [P40](#p40), [P17](#p17).

**Current state:** Reduced thermal/algebraic CGL/transverse filter are not these models. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/electrons/thermal_exchange.jl`, `src/electrons/ohms_law.jl`, `src/fields/field_state.jl`, `src/models/hybrid_pic.jl`.

**Deliver:** Separate derived contracts for momentum drag/recoil/Joule/thermal sharing, full electron momentum/bulk energy, and spatial evolving gyrotropic/tensor pressure with selected heat/isotropization terms.

**Close only when:** Independent drifting multi-ion equilibration and full momentum/energy/entropy rates, dynamic tensor compression/isotropic/CGL limits, electron-scale dispersion/Ohm/work/solver residuals, zero-mass/source limits and full persistent state; avoid double-counted heating and guard-only admission.

<a id="p42"></a>
### P42 — Finish whole-code audit and measured optimization

**Scope:** CAP-GLOBAL · **Acceptance dependencies:** [P01](#p01), [P02](#p02), [P03](#p03), [P04](#p04), [P05](#p05), [P06](#p06), [P07](#p07), [P08](#p08), [P09](#p09), [P10](#p10), [P11](#p11), [P12](#p12), [P13](#p13), [P14](#p14), [P15](#p15), [P16](#p16), [P17](#p17), [P18](#p18), [P19](#p19), [P20](#p20), [P21](#p21), [P22](#p22), [P23](#p23), [P24](#p24), [P25](#p25), [P26](#p26), [P27](#p27), [P28](#p28), [P29](#p29), [P30](#p30), [P31](#p31), [P32](#p32), [P33](#p33), [P34](#p34), [P35](#p35), [P36](#p36), [P37](#p37), [P38](#p38), [P39](#p39), [P40](#p40), [P41](#p41).

**Current state:** Final independent audit/resources remain open. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `src/`, `ext/`, `test/`, `docs/DEEP_DEBUG_TRACKER.md`.

**Deliver:** Audit every exported/public lifecycle and consumer after the last implementation; iterate confirmed defect repair; profile full operations and retained exceptional paths under fixed baselines.

**Close only when:** No unresolved confirmed defect in declared audited scope, every mandatory risk settled, meaningful independent wrong-method challenges, preserved complete compiler trees/numerical gates, fixed resource targets and all final compatibility rows rechecked.

<a id="p43"></a>
### P43 — Close compatibility, documentation and accepted release

**Scope:** CAP-GLOBAL · **Acceptance dependencies:** [P42](#p42).

**Current state:** Project-wide acceptance remains incomplete. **Operator status:** open; owner not assigned.

**Read/change canonical owners:** `README.md`, `docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md`, `.github/workflows/HybridPlasmaPIC-CI.yml`.

**Deliver:** Close all79 milestones and all required model/geometry/closure/species/source/precision/integrator/backend/restart rows, update examples/guides/schema/migrations and publish accepted main/release metadata.

**Close only when:** Every applicable G1-G10 gate has exact-head evidence; no required open/unknown/blocked row, no unsupported silent fallback; full required main CI green or evidenced billing-only exception, real hardware/science/resources separately complete and all docs match final behavior.

## Frozen controller and transport decisions

P30 first reads the [controller handoff](../handoff/2026-11-09/controller/HANDOFF_OPENCODE.md), [approved V3 derivation](../handoff/2026-11-09/controller/AMENDMENT_V3.md), final receipt and the entire patch. Complete root review, parser-fix reconciliation and normal module/IO/full/focused registrations precede acceptance. Its first domain is full-time projected Euler/Heun CIC/TSC with separate branch charge and actual scaled nonrelativistic/relativistic force; NGP jump/error rules, multirate cost and broader backends/compositions remain work. Test real error/effectivity and work, not merely the estimator's self-consistency. Fresh destination FFT plans/processes must establish exact retained-primitive replay in the permitted environment. This portability issue is an **unexecuted hypothesis**, not a measured plan defect; settle it without silent reprojection or weakening saved-phase semantics.

P40 reviews the [conduction proposal](../handoff/2026-11-09/flux-proposal/PRE_EDIT_CONTRACT.md) with its primary sources and proof notes. Before implementation settle saturation constants/units, odd plus Nyquist flux channels, variable divergence/geometry, active-zero states, anchoring and post-restoration actual cap, projected conductive entropy versus rounding/full-step entropy, and persistence. The proposed rejection of mixed zero active cells cannot silently become a universal cold-front capability. Add the missing physical law/domain extension or explicitly leave the full milestone open. Coordinate with P41 full momentum drag, electron inertia and evolving tensors; reduced thermal exchange is not full drag, and algebraic global CGL is not a transported pressure tensor.

## Scientific and performance completion

Use P25–P28 for pinned wave growth/damping/dispersion, nonlinear shock/reconnection/instability profiles, published datasets and independent solver comparisons over the actually implemented domains. Trace case 32's actual supplied resistivity/inertia before calling its mechanism collisionless: the inspected helper defaults to `eta=0.005`. Its causal mechanism is unmeasured; parameter, spatial/temporal and particle scans plus independent references must settle it. Do not rename a case solely because a ledger closes.

Fair comparisons match equations, units, geometry, resolution or measured error, particle noise, tolerances, output and hardware. Record initialization/FFT/setup versus steady evolution, communication, IO, peak memory and uncertainty. Preserve input decks, resolved dependencies, scripts, data licenses/provenance and raw results. Report the supported regimes and measured tradeoffs; the goal's superlative wording is not a claim established by repository size or passing test counts.

## Closure and release checklist

All 79 milestones have a stable crosswalk; every required item is implemented and all 44 work packages pass applicable gates. Complete the support matrix and source/API/checkpoint compatibility review after the last implementation. Perform the independent final whole-code audit against that exact tree; fix every confirmed defect in canonical owners and re-check changed dependencies until the finite audited scope has no unresolved confirmed defects. Keep limitations and unexecuted hypotheses explicit.

Run the final remote minimum/latest full matrix plus required real hardware and physical/resource campaigns. Read raw logs, including every changed whole-JET chain. Confirm all required statuses at the exact final source; perform docs/link/hash checks and dependency/release compatibility. Publish the accepted revision to main, read back remote head and its CI, and update the tracker with immutable evidence. Never mark the full goal complete just because implementation source has landed or CI is queued.

## OpenCode startup prompt

> Work from `implement/cgl-integration-20261109` and read `docs/OPENCODE_IMPLEMENTATION_HANDOFF.md`, `docs/OPENCODE_ACCEPTANCE_TRACKER.md` and `handoff/2026-11-09/README.md`. Verify the packet with `python3 handoff/2026-11-09/verify_handoff.py`. Preserve all dirty worktrees, protected files, the compiler ledger and inherited processes. Keep Julia/parser/compiler/formatter/tests/MPI/profiling off the shared PC; use GitHub or approved remote hardware. Start P00/P01 with both verified EMPIC parser errors, then fully dispose of the new JET/format/MPI failures without weakening gates. Assign sizeable independent packages to exclusive path owners, with one integration lead. Complete every required CAP-001–011 implementation and every mapped milestone, including real failure/rollback/restart and physical/resource acceptance. The controller is frozen but unadopted; reconcile the parser repair before reviewing/adopting it. The flux proposal is unapproved and requires its design/domain review before code. Push coherent candidates to validation branches; promote the exact source to main only after required green CI, and verify publication. Keep the MD/JSON trackers and docs current with raw evidence, report progress while waiting, and never substitute source presence for tested/scientific acceptance.
