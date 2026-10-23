# Plasma simulation audit and implementation tracker

Updated: 2026-10-24 (Asia/Taipei). Goal remains **active**.

User objective: audit the full codebase until an independent pass finds no new
confirmed defects; optimize measured bottlenecks; plan and implement missing
capabilities toward a leading, comprehensive plasma simulation tool. Compatibility with the latest stable Julia is an explicit acceptance gate. Every
verified fix or implementation is committed and pushed separately to `main`. Each push must
also pass GitHub CI; a confirmed billing block is recorded separately.
Passing a subset of tests does not establish whole-project correctness or
completeness. Comparative leadership requires reproducible external benchmarks.

## Current progress (2026-10-24)

- [x] Rechecked the [official Julia downloads page](https://julialang.org/downloads/):
  current stable is 1.13.2, released 2026-10-14. Native and complete CI logs
  verify execution on 1.13.2 alongside minimum supported 1.12.7.

- [x] GROWTH-001 pushed as `16c9776`; [all six exact-SHA jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39465803412).
  Downloaded logs verify 507,988 assertions on Julia 1.12.7 and 507,989 on
  1.13.2. Complete JET trees retain 23/22 diagnostic leaves against `b3e2eb6`;
  the only changed frame is a generated Base type-variable name. No findings
  were suppressed. The earlier `c4da551` documentation push also passes all six jobs.
- [x] Recovered CI logs and exact Git source into persistent evidence storage.
- [x] Reproduced COLL-001 against `16c9776` on both runtimes: 3,834 passes,
  127 behavior failures, no fixture errors in the reconstructed public suite.
- [x] Reconstructed COLL-001 passes final local gates on both runtimes:
  91,372 focused assertions, 2,358 isolated endpoint assertions, 54,959
  neighboring assertions, 34 concrete compiler checks and 24 scalar checks.
  All 120 ordinary full outputs and allocation counts match the parent per
  runtime. Full compiler groups match (30/27 reports); only group ordering
  differs on minimum Julia. Independent review found and repaired the public
  compact-range BigFloat conversion gap, then found no further confirmed defect.
- [x] COLL-001 pushed separately as `b251cb7`; all 351 committed Julia blobs
  match the validated candidate, remote main matches, and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for b251cb7 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39885899763).
  Downloaded complete logs verify 599,360 / 599,361 passing assertions on Julia
  1.12.7 / 1.13.2. Both complete compiler trees retain the reviewed 23 / 22
  diagnostic messages; changed source and dependency frames remain recorded.
  The full local Julia 1.13.2 run finished with 599,361 passes, zero failed
  assertions and one error in the metadata subprocess harness. Its child
  selects a manifest-free checkout instead of the resolved active environment;
  `full1132.log` retains the failed command and dependency-load error. A separate
  harness repair was pushed as `2993f5a`. This local full-suite run is not green.
- [x] Reconstruct and locally accept GROWTH-002 in a separate persistent worktree.
  Fresh serial probes retain 30 passes/30 failed integrity checks on both
  runtimes. Both reproduce negative MPI byte displacements. Frozen V2 passes
  all six focused MPI runs (2/4/8 ranks on both runtimes); a subsequent safe
  singleton witness finds a position-read exception outside the collective
  failure stage. V3 closes that path and passes all six focused MPI runs,
  8,495 direct assertions and 99,191 combined assertions per runtime. Complete
  single-thread compiler groups match the parent; the retained threaded entry
  context differences are reviewed in the plan. Both final full two-rank
  integration commands pass. All 34 command logs, frozen sources and unchanged
  committed Julia files are checked. Resource measurements and complete output comparisons are recorded
  in [the repair plan](MIGRATION_STORAGE_PLAN.md).
- [x] GROWTH-002 pushed separately as `7c6e1e4`; all 353 committed Julia files
  match the validated candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for 7c6e1e4 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39892881923).
  Downloaded complete logs verify 600,993 / 600,994 passing assertions on Julia
  1.12.7 / 1.13.2. Complete compiler trees match the reviewed CI trees with
  23 / 22 diagnostic messages. All six job-log hashes are retained.
- [x] Reconstruct ENERGY-010's Float32 field-energy accuracy repair in a separate
  persistent worktree. Fresh exact-`b251cb7` probes reproduce all four large-grid
  errors on both runtimes; 120 characterization assertions pass per runtime.
  Float64 and smaller-grid controls are exact. Frozen V5 now passes 18,862
  focused assertions, all 164 measured owner profiles, 50 concrete compiler
  checks and 210 concurrency assertions per runtime. Earlier failed candidates
  remain retained. The frozen composition onto `7be1916` passes 21,313 focused
  assertions per runtime and the same resource, mutant, concrete compiler and
  concurrency gates. Full package compiler trees match (30/27 reports). All
  11 native command logs, 272 frozen sources and 348 unchanged committed Julia
  files are checked. [Acceptance details](MAGNETIC_ENERGY_ACCURACY_PLAN.md).
- [x] ENERGY-010 pushed separately as `5c2796b`; all 357 committed Julia files
  match the validated candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for 5c2796b pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39897058507).
  Downloaded complete logs verify 614,784 / 614,785 passing assertions on Julia
  1.12.7 / 1.13.2. Both full compiler trees retain 23 / 22 diagnostic messages
  and match the reviewed parent after checked source-line shifts. All six
  job-log hashes and the archive hash are retained.
- [x] Locally accept the separately repaired periodic migration seam inconsistency.
  Both runtimes retain 144 control passes and 72 failed invariants against
  `b251cb7`; the public boundary folds the rounded upper endpoint while rank
  lookup and migration do not consistently do so. The separate candidate passes
  the original 216 assertions, 27 concrete compiler checks, 54,480 neighboring
  assertions and real 2/4/8-rank routing on both runtimes. All 108 ordinary
  complete-output/allocation samples match per runtime; each resource command
  passes 667,584 assertions. Complete package compiler groups match (30/27).
  [Separate repair plan](PERIODIC_MIGRATION_PLAN.md).
- [x] Periodic seam repair pushed separately as `e8925e8`; all 355 committed
  Julia files match the validated candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for e8925e8 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39893248890).
  Downloaded complete logs verify 607,041 / 607,042 passing assertions on Julia
  1.12.7 / 1.13.2. Complete compiler trees match the reviewed CI trees with
  23 / 22 diagnostic messages. All six job-log hashes are retained.
- [x] Locally accept the measured interval-validation allocation improvement.
  Both runtimes pass 21,723 neighboring assertions and 40 concrete compiler
  checks. All 270 complete public-output/RNG samples match per runtime, with
  zero candidate allocation in all samples. Complete package compiler groups
  match the parent (30/27). [Plan and evidence](INTERVAL_VALIDATION_ALLOCATION_PLAN.md).
- [x] Interval validation pushed separately as `7be1916`; all 356 committed
  Julia files match the validated candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for 7be1916 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39895197004).
  Downloaded complete logs verify 609,492 / 609,493 passing assertions on Julia
  1.12.7 / 1.13.2. Complete compiler trees match the reviewed CI trees with
  23 / 22 diagnostic messages. All six job-log hashes are retained.
- [x] Repair the metadata subprocess environment; both runtimes reproduce the
  dependency-load error from a resolved external environment. Keep both child
  processes and all checksum assertions; source identity is checked explicitly.
  Both runtimes pass all 80 metadata assertions and six rejection controls in
  both external and checkout environments, then pass the same checks on the
  final `5c2796b` composition. [Plan](TEST_SUBPROCESS_ENVIRONMENT_PLAN.md).
- [x] Metadata harness repair pushed separately as `2993f5a`; all 357 committed
  Julia files match the validated candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for 2993f5a pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39897483287).
  Complete logs verify 614,784 / 614,785 passing assertions on Julia 1.12.7 /
  1.13.2. Both full compiler trees match the field-energy CI baseline (23 / 22
  reports); all six log hashes are retained.
- [x] Locally accept PIC-007, the separately reproduced legacy longitudinal-current
  modal quotient loss. Both runtimes pass 24 baseline characterization checks.
  V1/V2/V3/V4 retain failed cost matrices. V4 passes 433,206 guard comparisons
  and 1,438 numerical assertions on each runtime, but its minimum-version
  matrix has one latency failure. V5 replaces the minimum operations with
  equivalent direct comparisons and passes the same guard and numerical gates
  on both runtimes, then all 200 cost/state assertions per runtime across 40
  profiles with zero allocation. Both final lanes pass 37,916 complete EM
  neighboring assertions, six concrete compiler checks and 96 exceptional
  output/state assertions. Complete package compiler trees match the parent
  (30/27 reports). All 358 committed Julia blobs match the validated candidate;
  held changes are unchanged. [Bounded repair plan](EMPIC_LEGACY_CURRENT_MODE_PLAN.md).
- [x] PIC-007 pushed separately as `7b50b51`.
- [x] [All six exact-SHA CI jobs for 7b50b51 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39905400737).
  Complete logs verify 616,222 / 616,223 passing assertions on Julia 1.12.7 /
  1.13.2. Both full compiler trees match the metadata-repair CI baseline
  (23 / 22 reports); all six log hashes and the archive hash are retained.
- [ ] Reconstruct the remaining ENERGY-004 complete-budget composition against
  the accepted field-energy source; retain the earlier rejected cost matrices.
  The isolated implementation and numerical checks are in progress. Root has
  reviewed validation, canonical component ownership and exact-offset completion;
  both timing-harness controls now pass 12,673 assertions and 3,072 pairs per
  runtime. These validate the measurement protocol, not candidate speed.
  The expanded numerical tail passes 678 assertions on each runtime. The first
  zero-state resource gate retains 307 passes and 72 allocation failures per
  runtime. V2 passes all 379 zero-state and 65 signed-state resource assertions
  on each runtime, matching parent bytes in 126 zero-state and 16 signed
  profiles. Both versions pass 2,346 numerical and 5,760 lower-bound checks;
  the timing screen rejects the third minimum-version profile at median ratios
  1.683 / 1.651. The remaining 57 profiles and latest lane are unrun. V2 remains
  rejected; further optimization, full integration and CI remain open.
  V3 passes the expanded untimed gates but its fixed timing screen rejects
  Float64/4/unit/iso at ratio-of-medians 1.524693 (limit 1.5); paired median
  1.458859 passes. Thirty Float32 profiles pass; the remaining 29 minimum
  profiles and all 60 latest profiles are unrun. V4 passes all untimed gates
  and 66 concrete compiler checks per runtime, but its fixed screen rejects
  Float32/4/unit/polyhalf at median ratios 1.525357 / 1.535156. The remaining
  57 minimum profiles and all latest profiles are unrun. Complete-owner
  profiling informed V5 arithmetic changes. A subsequent public zero-sign
  probe retains 54 passes / 18 failures: exact quadratic cancellation can return
  negative zero with a zero polytropic reservoir. The canonical offset branch
  is repaired in V6; no V5 timing was accepted. Both V6 lanes pass 16,030 focused
  assertions, the unchanged 72-assertion sign probe, 5,760 endpoint checks,
  379 / 65 / 121 zero/signed/small resource assertions, 193 kinetic-underflow
  checks and 66 concrete compiler checks. The fixed small-profile timing screen
  passes all 60 profiles and 11,761 assertions per runtime, with unchanged
  source/protocol hashes. Maximum median ratios are 1.485094 / 1.459885
  (limit 1.5). Broader cost, threading and cold-path gates remain open.
- [x] Locally accept the NGP nearest-node boundary repair on `bf2a8bc`.
  Both CPU runtimes reproduce 14 baseline failures; actual Metal retains 13.
  The shared upward-tie owner passes 7,472 focused, 46 concrete compiler,
  1,410 actual-Metal integration and 204 additional scalar-type assertions
  per runtime. Prior complete CPU coupling/solver integration passes 49,590
  assertions per runtime; checked source composition records the sole added
  CPU difference (the accepted unused-helper removal). Complete final compiler
  trees match (30 / 27 reports). All 126 fixed cost profiles per version preserve
  outputs and 0-byte allocation; maximum median ratios are 1.210 / 1.195,
  within the fixed 1.5 gate. All controls pass. Fifteen final command logs,
  276 frozen Julia files and 356 unchanged committed files are checked.
  [Plan and evidence](NGP_NEAREST_NODE_PLAN.md).
- [x] NGP pushed separately as `df75656`; all 361 committed Julia files match
  the validated candidate, remote main matches and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for df75656 pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39912564815).
  Complete logs verify 623,694 / 623,695 assertions on Julia 1.12.7 / 1.13.2.
  Full compiler trees match the Metal uniform-push CI baseline (23 / 22 reports).
- [x] Locally accept the separate Metal uniform-push repair. Both supported
  runtimes reproduce the parent's recursive dispatch and stack overflow from an
  initialized state. The kernel reuses the canonical scalar Boris update and
  validates host parameters before launch. Each runtime passes 778 focused
  actual-device assertions, 894 full backend assertions, 57,583 CPU neighboring
  assertions, 14 concrete compiler checks and 60 state/resource assertions.
  Complete package compiler trees match (30 / 27 reports). All 274 frozen Julia
  files, 355 unchanged committed files and 14 command logs are checked; held
  changes are unchanged. [Plan and measured limits](METAL_UNIFORM_PUSH_PLAN.md).
- [x] Metal uniform-push repair pushed separately as `bf2a8bc`; all 359 committed
  Julia files match the candidate and held hashes are unchanged.
- [x] [All six exact-SHA CI jobs for bf2a8bc pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39909168466).
  Complete logs verify 616,222 / 616,223 passing assertions on Julia 1.12.7 /
  1.13.2. Both full compiler trees match the legacy-current CI baseline
  (23 / 22 reports); all six log hashes and the archive hash are retained.
- [ ] Repair Metal gathered-field preflight. Both actual-device runtimes retain
  114 passes and 51 failed assertions, with all 45 finite-input CPU/device
  controls passing and no fixture errors. Invalid host parameters and late
  nonfinite gathered samples bypass CPU rejection semantics and can mutate
  particles and midpoint output. **VERIFIED:** the Metal owner uses raw `q/m`
  and `T(dt)` and has no finite-sample preflight (`ext/HybridPlasmaPICMetalExt.jl:643-675`);
  `metal-gathered-validation/baseline*.log` retains the public witnesses.
  V1 passes the original 165 and full 6,394 integration assertions on each
  runtime. Both parent differential preflights pass all 120 alias
  checks and retain 234 passes / six resource failures: the 262,144-particle
  profiles add at least 9,344 host bytes, exceeding the fixed 8,192-byte limit.
  Full allocation profiling attributes all additional sampled bytes to the
  reduction. V2 uses a per-call atomic flag; its independent device probe passes
  344 assertions per runtime. The full V2 sequence retains a test-loading error
  (optional macro expanded before import). The corrected fixture passes the
  original 165 and full 6,426 integration assertions on each runtime. All 120
  alias checks and 240 untimed cost/state checks pass per runtime; maximum
  extra host allocation is 4,048 / 4,016 bytes, within the unchanged 8,192-byte
  limit. Full package compiler trees match (30 / 27 reports). Complete concrete
  device captures retain 24 matching parent diagnostic leaves plus three
  upstream status-transfer/kernel-launch reports; every full signature and
  caller path is retained and reviewed. The fixed minimum-version timing matrix
  passes all four identical/doubled controls but retains 254 passes and 46
  latency failures; all 30 profiles are recorded. Latest-version timing is
  unrun. V2 is rejected. A separate six-profile fence intervention preserves
  complete outputs (54 checks); its observations motivate further profiling,
  not a relaxed gate or acceptance claim. V3 queues conditional device commits
  after validation and reads status at the end. Both runtimes pass 6,459 full
  integration checks and 2,034 expanded alias/rejection checks; additional host
  bytes are at most 1,088 / 1,040. Full package compiler trees match, and all
  device diagnostics are retained/classified. V3 still fails 44 latency assertions
  across 23 minimum-version profiles; latest timing is unrun. V4 reuses the
  canonical scalar Boris update in one guarded kernel for disjoint arrays and
  retains the existing ordered path for aliases. Its gates are running; no
  gathered-push candidate is accepted or pushed.
  [Repair plan](METAL_GATHERED_VALIDATION_PLAN.md).
- [x] Locally accept the shock CN field-solve repair. Both runtimes retain the
  parent's 66 control passes and 50 failures across 28 direct/public profiles.
  **VERIFIED:** public resistive decay at `dt=4` gives maximum field 3.125
  versus direct CN 0.533709 (`shock-cn-solve/baseline-*.toml` and logs).
  The residual-checked matrix-free solve and full-step transaction pass 666
  broad and 24 independent boundary assertions per runtime, including
  solver-history aliases, stationary fields, absorption/growth and restart.
  All RK4 state/allocation comparisons match. Four first-matrix timing flags
  pass the one prescribed reverse-order confirmation within 1.10; both matrices
  are retained. Fresh matched compiler records match 28 minimum and all 27
  latest records, with two removed and one reviewed generic KdV record on
  minimum Julia. Matched concrete KdV/FFT checks pass 86 assertions and the
  existing 61-assertion suite in all four processes; complete outputs match.
  All 363 Julia files in the final composition match the validated candidate.
  [Contract, limitations and measured repair costs](SHOCK_CN_SOLVE_PLAN.md).
- [x] Shock CN repair pushed separately as `ce52836`; all 363 committed Julia
  files match the accepted candidate and held hashes are unchanged.
- [ ] [Exact-SHA CI for ce52836](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39999012813)
  has green formatting and 2/4-rank MPI jobs; both Julia suites and 8-rank MPI
  remain running. All six jobs are required.
- [ ] Resolve the exported memory-status constructor invariant. Both runtimes
  verify that keyword/full-positional construction accepts negative byte counts
  rejected by the checked wrapper (185 passes, ten failed rejection assertions).
  All current telemetry callers use the wrapper. The repair plan now explicitly
  applies that existing nonnegative-integer-or-missing contract to every public
  constructor, with the intentional raw-constructor compatibility tightening
  documented. Fresh parent characterization and implementation remain open.
- [x] Reconcile external-validation READMEs with the checked runners: missing
  reference output fails; NHDS compares the eigenvalue oracle; VPIC compares
  scalar metrics. Record current gates and unresolved source/cache provenance
  instead of asserting fresh historical measurements. This documentation
  revision `61356f1` and the earlier plan revision `ff4a973` have all six
  exact-SHA CI jobs green. Complete logs on each verify 623,694 / 623,695
  assertions and unchanged full compiler trees (23 / 22 reports).
- [ ] Complete remaining audit/implementation tracks listed below.

Recovery limitation: the earlier `/tmp` candidate and evidence directories are
no longer present. Old process handles are unavailable. The historical records
below describe prior observations; their local artifacts cannot currently be
rechecked and are not acceptance evidence for reconstructed code. New source,
logs and manifests live outside `/tmp` at
`../HybridPlasmaPIC-audit-evidence/2026-10-23`; the ionization worktree is
`../HybridPlasmaPIC-ionization-20261023`. Existing held radiation changes remain
excluded from this repair.

## Historical evidence (artifacts subject to the recovery limitation)

- **VERIFIED — PERF-006:** the isolated rank-lookup candidate passes 46,188
  focused, 16 concrete compiler and 137 MPI/quality integration checks per
  runtime. Complete migration output fingerprints match across all 35 samples
  per runtime while allocation falls by 59.68–84.19% on Julia 1.12.7 and
  47.04–76.16% on 1.13.2. No elapsed speedup is claimed. Complete package compiler
  trees match baseline (22/21 reports, no added or removed groups). Pushed as
  `b3e2eb6`; all 11 committed files and remote main matched the promotion manifest
  at that push.
  Its [exact-commit CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39432448426)
  passes all six jobs. The package logs verify actual Julia 1.12.7/1.13.2 and
  504,062/504,063 assertions. Complete CI compiler trees retain 23/22 reports
  with no added or removed groups. All six downloaded log hashes were checked.
  [Evidence and contract](RANK_CLASSIFICATION_ALLOCATION_PLAN.md).

- **VERIFIED — GROWTH-002:** serial migration can partially replace particle
  storage before a late view conversion fails, and periodic position wrapping
  can change aliased weights from 1.25 to 0.25. Both runtimes retain 36 control
  passes and 24 failed integrity checks. The append repair does not cover this
  separate path. The first private serial candidate passes 216 checks per
  runtime. A fresh repeated-owner probe retains six failures per runtime in
  both baseline and that candidate. V2 then passes 240 initial and 870 expanded
  serial checks per runtime, plus focused real 2/4/8-rank rejection/recovery
  runs on both runtimes. Byte-count/displacement preparation is a further
  confirmed gap. V4 passes 1,262 focused checks per runtime, but its first
  concrete compiler gate retains six failed assertions. Matched captures on
  both runtimes verify that the serializer diagnostics already exist on the
  baseline, with only the expected preparation-wrapper frames added.
  V4 also passes all six focused real-MPI runs across both runtimes; full
  fixture-loaded and unloaded package trees match baseline on both runtimes.
  Both full two-rank suites pass, with
  all eight MPI command/log hashes and frozen sources checked. V5 preserves
  the revised growth/rank source; its seven-path union is verified, and native
  checks on the combined snapshot are running.
  [Evidence and repair plan](MIGRATION_STORAGE_PLAN.md).

- **VERIFIED — GROWTH-001:** public append probes reproduce partial growth
  after a late nonresizable field and inconsistent lengths with shared storage
  on both Julia versions. The corrected ionization fixture also reproduces
  partially changed populations; V2 retains 172 passes and 48 failed invariants
  per runtime, with all RNG/success controls passing. The repair contract and
  preserved failures are in the [particle-growth plan](PARTICLE_APPEND_PLAN.md).
  Candidate V1 passes all 220 original checks on minimum Julia, but its
  allocation gate and two expanded-test fixtures require correction. V2 fixes
  the measured resize-check allocation; V3 repairs the fixtures without
  changing V2 production code. Both V2 runtimes pass 226 original/resource
  smoke checks and 3,060 retained-capacity append checks. V3 passes 1,784
  focused checks on each runtime. Tiny
  append costs increase and are recorded in the plan. A fresh finite-state
  cross-species backing-storage probe then reproduces 24 baseline and 12 V3
  failed assertions per runtime. V3 and its first rank-integration snapshot
  must not be promoted. The revised combined candidate preserves both rank
  and append coverage and passes 4,236 focused checks on both Julia versions;
  both also pass 47,838 neighboring integration and 52 concrete compiler
  checks. All 120 ordinary-birth output samples per runtime match V3 exactly,
  with identical allocation bytes and zero measured compilation/recompilation.
  Complete unloaded and fixture-loaded compiler trees match committed
  `b3e2eb6` on both runtimes (22/21 reports, no changed groups); all 92 reference
  source/extension blobs match that commit. The old-parent comparison retained
  a path-normalization failure; a same-file proof and corrected commands are
  recorded in the plan; all four corrected comparisons now pass. Both original
  neighboring suites pass 43,941 assertions. Restoring a missing historical
  baseline from its exact Git commit also completes all four original compiler
  comparisons, with no changed report groups. The revised seven-file repair
  passes its local gates and is pushed separately as `16c9776`;
  [exact-commit CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39465803412)
  remains pending. Original failures and rejected candidates remain retained.

- **VERIFIED — COLL-001:** ionization can miss eligible births, leave primary
  energy unchanged, or select an event with an erroneously saturated probability
  when finite velocity components have an unrepresentable norm. Both Julia
  versions reproduce 30 energy/threshold failures and 18 probability/state
  failures, with all independent arithmetic/RNG controls passing.
  [Evidence and repair plan](IONIZATION_ARITHMETIC_PLAN.md). The corrected
  near-threshold probe adds 19 reproduced failures per runtime, with 77
  arithmetic/RNG/control passes and no fixture errors. It also demonstrates
  a Float32 birth at exact threshold equality and missed Float64 births just
  above threshold. The original ambiguous RNG fixture remains retained. The
  V2 repair passes 53,519 focused and 5,879 integration checks per runtime;
  complete compiler trees remain unchanged, and ordinary full-call outputs and
  allocations match the parent. A separate BigFloat
  endpoint probe then reproduces 12 failures per runtime in the new primitive,
  with 78 other assertions passing. V1/V2 are rejected. V3 passes those
  endpoints but a wider-precision probe exposes six integer-intermediate
  failures per runtime; V3 is also rejected. V4 passes 53,521 focused checks
  per runtime plus 252 child-process endpoint checks. Both integration suites
  pass 5,879 assertions. Ordinary outputs and allocations match the parent,
  while exact-fallback costs are recorded separately. Both concrete compiler
  gates and all four complete package comparisons pass. Its 92-file comparison
  parent matches pushed `16c9776`. Parent CI and separate promotion/CI remain
  required; the arithmetic repair has not been promoted.

- Documentation checkpoint `0da4e09` is pushed to `main`; all six committed
  documentation blobs and the remote SHA match the promotion manifest.
  Source/test trees match `ffdb522`; held radiation files are excluded.
  Its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39401716377)
  passes all six jobs: 502,711/502,712 assertions on actual Julia
  1.12.7/1.13.2, formatting and all MPI ranks. Complete compiler trees match
  the accepted baseline (23/22 reports, no added or removed groups). No billing
  exception is used.

- Inertia and zero-state resource repair is pushed as `ffdb522`; all 16
  committed blobs and remote main were verified after promotion. Its
  [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39150381077)
  passes all six jobs: 502,711/502,712 assertions on actual Julia 1.12.7/1.13.2,
  formatting and all MPI ranks, with unchanged complete compiler trees (23/22).
  Complete-total composition is the next ENERGY-004 implementation, with its
  [pre-edit contract](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#complete-total-evaluator-pre-edit-brief)
  recorded before candidate edits.

- Complete-budget candidates remain isolated. V8 passes 13,636 permanent/Aqua,
  770 threaded, 12,689 shared-field, 3,681 solver, 453 shock and 38,055 checkpoint
  checks per runtime; fixture-loaded full compiler trees remain unchanged.
  V7's controlled cost gates retain all values and bytes, but one of 324 minimum
  Julia profiles reaches 1.536× and one of 16 signed latest-Julia profiles reaches
  1.515×, above the 1.5 target. V9 bound reuse and the wider reduction pass both profiled cases at
  0.94–1.05×. V10 passes permanent/compiler/threaded/field gates but its component
  regression catches a stale private test alias. V11 migrates that caller,
  preserves its limits and passes 49,321 component, 14,128 permanent/Aqua,
  770 threaded and 12,689 shared-field checks per runtime, with unchanged
  complete compiler trees. Its first timing matrix fails three profiles
  (1.572× minimum Julia; 2.218×/1.569× latest). V12 shares closure validation
  after profiling and passes formatting/migration, but its small timing matrix
  fails four/five latency checks on minimum/latest Julia with all values and
  bytes unchanged. Identical-function controls then fail eight/two timing
  profiles on minimum/latest Julia with no compilation or GC time. A longer
  calibrated protocol remains unaccepted: its first prototype has a slowdown
  wrapper allocation mismatch and an identical-call timing failure.
  The corrected 50 ms/48-pair prototype passes 16 numerical/allocation
  preflight assertions on each runtime, but also fails identical-call timing
  bounds on minimum Julia (9,342 passes, two failures). Latest Julia passes
  all 9,344 control assertions. The protocol remains rejected. Thread-clock
  diagnostics pass 436 assertions per runtime: across six selected cases,
  CPU-time ratios span 0.982–1.016 while elapsed ratios span 0.793–1.369.
  These diagnostic observations do not accept the candidate or establish an
  external process as the cause.
  No composition change has been promoted.
  [Evidence and contract](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#scalable-total-certificate-pre-edit-brief).

- **VERIFIED — ENERGY-010:** exact-`ffdb522` magnetic-energy probes on both
  Julia versions return 3,092,971.5 for a uniform Float32 160³ field whose
  rational integral is 3,200,000 (3.34464% error). All eight Float64 controls
  are exact. The ordinary scalar reduction reproduces each public result.
  [Evidence and implementation plan](MAGNETIC_ENERGY_ACCURACY_PLAN.md).

- **VERIFIED — remaining ENERGY-004:** 316 baseline witness assertions per runtime establish complete-budget and inertia range failures. Inertia V2 now passes all local gates: 2,846 focused, 1,987 budget, 12,689 shared-field and 666 cost assertions per runtime, unchanged complete compiler trees and formatting. The rejected 1.32–1.33 GB zero-state allocation is reduced to 3.15/6.29 MB at 64³. Its exact-SHA CI now passes all six jobs; complete-total composition remains required. [Evidence and phased plan](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#complete-budget-composition-and-inertia--current-energy-004-scope).

- Exact-node reuse passes canonical local acceptance on Julia 1.12.7/1.13.2: 25,091 focused, 462 shock/Aqua and 480 cost assertions per runtime, complete compiler trees (22/21 unchanged), migration and formatting. Measured large cold allocation is 50.4–50.8% of the pushed owner, with ordinary outputs/allocation unchanged. The earlier tuple-return regression is retained. Pushed as `99760f0`; its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39147850360) passes all six jobs, including 499,865/499,866 assertions on actual Julia 1.12.7/1.13.2 with unchanged complete compiler trees (23/22). [Contract and evidence](SHOCK_ENERGY_FLUX_RANGE_PLAN.md#exact-node-reuse-optimization-brief).

- Complete shock totals are pushed as `dd539da`. Its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39146029176) passes all six jobs, including 499,511/499,512 assertions on actual Julia 1.12.7/1.13.2 and unchanged complete compiler trees (23/22). Both native versions pass all focused, resource, shock, Aqua and complete compiler comparisons. All 12 committed blobs and the remote SHA match the final promotion manifest; the held radiation work is excluded.

- Shock component range repair is pushed as `daea194`. Its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39144368129) passes all six jobs: 491,943/491,944 assertions on Julia 1.12.7/1.13.2, formatting and all MPI ranks. Full compiler trees remain unchanged at 23/22. All 14 committed blobs and remote main matched the promotion manifest. ENERGY-009 separately tracks complete-total arithmetic.

- Complete polytropic energy is pushed as `7104c61`, with all six [exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39141208484) green: 474,774/474,775 assertions on actual Julia 1.12.7/1.13.2, formatting and all MPI ranks. Complete compiler trees remain unchanged at 23/22. Both native integration lanes also pass: 3,681 solver, 453 shock and 38,055 checkpoint assertions each. [Full evidence](CI_RELIABILITY_PLAN.md#complete-polytropic-energy-integration).

- Algebraic dispatch/compiler correction and CI test repairs are pushed as
  `15662e1`, with [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39137529999)
  fully green: 446,514/446,515 assertions on actual Julia 1.12.7/1.13.2, formatting and all MPI ranks. Complete CI compiler trees remove the two targeted reports with no added groups. The nine-file promotion includes six checked source/test paths and
  three trackers; held radiation files are excluded. Full native/compiler
  evidence and explicit remaining gates are recorded below. No pending or
  failed run is classified as green.

- Isothermal correction `6e7fb01` is pushed to `main` with its 15 checked Julia
  source/test paths and tracker documentation. Both native versions pass 15,427
  permanent assertions plus a smoke check and 72 concrete compiler checks.
  Formatter/documentation changes preserve all 205 source/test executable syntax
  trees. [Exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129297310)
  has finished with the fixture failure recorded below; corrected-source native broad checks pass and replacement CI acceptance stays open.
  Held radiation files are excluded. The traversal optimization is pushed as
  `5323c4d`, with its [own CI run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129564503)
  fully green. Pressure correction `0d8c8e8` is also pushed with
  [its own run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129830419).
  It passes 27,704 focused plus 1,314 independent assertions per native runtime,
  and 172/156 concrete compiler checks. Native broad acceptance continues on the final repaired source; the two failed
  historical revisions retain their exact outcomes.


- Two exact-SHA latest-Julia package jobs have failed: `0d8c8e8` has 14
  allocation-byte equality failures after 446,177 passes; `6e7fb01` has one
  uninitialized-position fixture error after 404,763 passes. The
  `6e7fb01` minimum-version job passes 404,768 assertions. The `0d8c8e8` minimum-version job has finished with 446,175 passes and 15 failures in the same two byte-equality assertions. Both historical runs remain failed. All six jobs for `5323c4d` are successful.
  TEST-002/003 track the corrections and retained logs; no billing exception
  applies. [Evidence](CI_RELIABILITY_PLAN.md#october-18-failed-package-checks-and-corrections).

- **VERIFIED — CLOSURE-002:** the pushed pressure owner rejects union-element
  input/output arrays in 15 of 27 storage/coefficient combinations where the
  legacy owner succeeds. Both runtimes pass 69 witness assertions. The isolated
  first dispatch correction passes 316 storage/broadcast assertions per runtime
  but retains two full-package diagnostics. V2 restores complete compiler trees
  to the 22/21 baseline and passes concrete, threaded, allocation and Aqua
  gates on both runtimes. Each native focused suite passes 57,165 assertions. The correction and both CI test
  repairs are implemented; exact-SHA GitHub acceptance now passes in `15662e1`.
  [Evidence and correction](POLYTROPIC_ENERGY_RANGE_PLAN.md#algebraic-storage-dispatch--closure-002).

- The isolated complete polytropic cold-integral component passes 273 checks on
  each native runtime, including the huge-gamma midpoint witnesses, plus 12
  compiler and 195 threaded assertions. The ordinary component passes 1,429
  numerical and 351 quality/resource/threaded checks per runtime, with 54
  measured calls allocating zero bytes. Gamma-one V4 passes 1,006 numerical,
  151 quality/resource/threaded and 28 general-resource checks per runtime.
  Canonical V3 passes 70,219 focused assertions per runtime, with unchanged
  complete compiler trees and passing concrete/concurrency gates. V4 exposes
  a 2.27 GB cumulative-allocation distinct-density case. The isolated fixed-log
  correction passes 15,086 component, 862 retained gamma-one, 120 public resource
  and 223 quality/concurrency assertions per runtime; all 24 measured public
  profiles allocate zero bytes. Canonical V5 passes 85,433 focused assertions per runtime, both complete
  compiler comparisons (22/21), 1,679 quality/Aqua and 1,608 threaded assertions
  per runtime, plus 223 direct fixed-log checks and 577 independent public
  range checks per runtime. The reviewed 20-file source/test overlay is integrated; both native broad lanes and exact-SHA GitHub acceptance now pass. [Scope and evidence](POLYTROPIC_ENERGY_RANGE_PLAN.md).

- VALID-001 charge conservation passes local acceptance on **1.12.7 and 1.13.2**:
  29,113 complete solver/restart assertions, 13,268 numerical/compiler/quality
  assertions, 2,273 concurrent checks and 504 resource/oracle checks per version.
  The unchanged 1,000-step wave passes at 2.88e-12 against 1e-8. The exact tested
  source is pushed as `e9c4fa7`; [its CI run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38767093527)
  is fully green: 368,455/368,456 assertions on actual Julia 1.12.7/1.13.2,
  plus formatting and all MPI ranks. [Contract and costs](CHARGE_DEPOSITION_CONSERVATION_PLAN.md).

- Latest-version gate advanced to **Julia 1.13.2**, listed as the current stable
  release (October 14) on the [official downloads page](https://julialang.org/downloads/manual-downloads/),
  checked again October 20. Existing 1.13.1 results retain their exact tested version;
  the latest fully green pushed revision is `ffdb522`. Its [exact CI run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39150381077)
  passes all six jobs: 502,711/502,712 assertions on actual Julia 1.12.7/1.13.2,
  formatting and all MPI ranks. No billing exception is used.
  The minimum-version gate remains 1.12.7.

- Pushed revisions `afc324f`, `b73abfe`, `2330379`, and `2bc7085`
  each pass both package jobs, formatting and
  2-/4-/8-rank MPI. That revision passes 369,711/369,712 assertions on actual
  Julia 1.12.7/1.13.2. No billing exception is used.
  [Exact revisions and logs](CI_RELIABILITY_PLAN.md#acceptance-run-at-afc324f).

- Pushed `240f77c` now has fully green [GitHub run 38741178196](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38741178196):
  actual Julia **1.12.7 passes 357,258 assertions** and **1.13.2 passes 357,259**;
  formatting and 2-/4-/8-rank MPI jobs also pass. Package-analysis counts are
  22/21, respectively. This establishes latest-stable acceptance for this SHA;
  later revisions retain their own pending gates.

- Pushed `4cd642c` has fully green [run 38750764319](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750764319):
  actual Julia 1.12.7 passes 358,978 assertions and 1.13.2 passes 358,979, with
  23/22 JET reports. Formatting and all MPI ranks pass. Pushed `6ad7276` also has
  a fully green [run 38754771205](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38754771205):
  359,870/359,871 assertions on actual 1.12.7/1.13.2, with 23/22 JET reports.
  Scheduling repair `4c03c87` also has a fully green
  [run 38750474902](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38750474902):
  357,258/357,259 assertions on actual 1.12.7/1.13.2, with 23/22 JET reports.
  The later `e9c4fa7` run also passes as recorded above.
  No billing block is evidenced. [Per-revision evidence](CI_RELIABILITY_PLAN.md).

- Docs revision `e339350` has fully green
  [run 39091128686](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39091128686):
  369,711/369,712 assertions on actual Julia 1.12.7/1.13.2, plus formatting and
  all MPI ranks. Revision `1150ed3` also has all six jobs green in
  [run 39093601976](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39093601976):
  369,711/369,712 assertions on actual Julia 1.12.7/1.13.2, with 23/22 JET reports.
  No billing exception is evidenced.

- The evolving-pressure implementation passes local acceptance on Julia
  **1.12.7 and 1.13.2**: 2,290 broad regression assertions, 336 resource checks,
  996 final energy/concrete/Aqua checks and 650 four-thread checks per runtime.
  Real two-rank MPI reruns pass on both. Complete final compiler report trees
  match the baseline (22/21 reports); concrete field-owner checks add no reports.
  Independent pressure physics, 3D transport, restart and rollback evidence is
  in the [implementation plan](HYBRID_ENERGY_MODEL_PLAN.md#evolving-pressure-implementation-tracker).
  The [user guide](EVOLVING_ELECTRONS.md) describes the ideal massless scope and
  positivity/timestep limits. Its example runs on both versions. Documentation
  and formatting preserve executable ASTs across all 186 candidate Julia files.
  Pushed as `013a70b`; [its exact-revision CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39097756477)
  has all six jobs green: 374,060 assertions on actual Julia 1.12.7 and 374,061
  on 1.13.2, plus formatting and all MPI ranks. Package JET counts remain 23/22.

- Completed [GitHub run 38722319490](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38722319490)
  at `55d996b` passes formatting, the 2-/4-rank MPI jobs and both package jobs:
  **352,519 assertions each**, in 66m46.4s on Julia 1.13.1 and 137m59.4s on
  1.12.7; JET reports 21/22. Its only failure is the 8-rank step's explicit
  20-minute timeout after passing checkpoint publication and moment-buffer
  validation. The bounded timeout repair and subsequent local repairs require
  acceptance in the next run; no billing block is evidenced. The combined native source passes JET
  (22/21 reports), Aqua and three concrete repair analyses on both versions.

- GitHub run [38710456393](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38710456393)
  at `fadc03f` is not green: Julia 1.13.1 cannot resolve the JET 0.9–0.11
  constraint against the current registry; the 8-rank checkpoint publication
  test reports missing rank files. Formatting and 2-/4-rank jobs pass; the
  Julia 1.12 package job ends with 352,282 passes and five failures: four
  RayCon normalization comparisons and one KdV error-ordering assertion. JET
  and Aqua pass. No billing block is evidenced.
  The local repairs pass both native-version gates and the delayed 8-rank
  replay. Fixes `fadfbf1` (MPI ordering) and `f0c9ce9` (JET compatibility/import)
  are included in that push; final GitHub acceptance remains pending. Logs: `/tmp/hybrid-ci-38710456393-julia-latest.log` and
  `/tmp/hybrid-ci-38710456393-mpi8-artifact/mpi-smoke.log`.

- The immutable `6341c4b` full `Pkg.test()` suites pass **352,287 assertions
  each** on native Julia 1.13.1 (132m00.6s) and 1.12.7 (160m03.1s), both exit 0.
  All 266 recorded Julia and environment files per version retain their snapshot
  hashes. These snapshots include archive inference, checkpoint serialization
  and continuity repairs; later KdV, RayCon, threaded-alias and held radiation
  work are excluded. Logs: `/tmp/hybrid-full-6341c4b-julia{112,113}.log`;
  source/environment hashes: `/tmp/hybrid-full-6341c4b-manifest.json`.

- Full package suites finished at immutable `fe519e9`, each with a separately
  resolved native environment. Julia 1.13.1 passes **338,415 assertions** in
  55m10.4s (exit 0). Julia 1.12.7 passes **338,414**, with one failed JET gate in
  68m53.2s (exit 1): 31 reports exceed the unchanged ceiling of 30. These include
  the electromagnetic field-energy and serial-publication fixes. Logs:
  `/tmp/hybrid-full-fe519e9-julia112.log` and
  `/tmp/hybrid-full-fe519e9-julia113.log`; snapshot hashes are in
  `/tmp/hybrid-full-fe519e9-manifest.json`. Later continuity and checkpoint
  serialization repairs, and uncommitted radiation work, are excluded.
  JET-002 tracks the failed integration gate; latest-version success does not
  establish that the minimum-version suite is green.
- The immutable `d2c3da7` full package suite passed **294,035 assertions** on
  each of Julia 1.13.1 (62m33.3s) and Julia 1.12.7 (82m43.6s); both processes
  exited successfully. The real two-rank MPI run also passed and exited successfully.
  Each version resolves its own environment. Logs:
  `/tmp/hybrid-full-d2c3da7-julia112.log`,
  `/tmp/hybrid-full-d2c3da7-julia113.log`,
  `/tmp/hybrid-full-d2c3da7-mpi2-julia112.log`. These runs include the clock,
  public-drift and moment-buffer repairs; later work is excluded.
- Starting revision: `905ea1b`; clean worktree; remote `main` matched that revision.
- Inspected package exports, manifests, CI, contribution rules, README, tests,
  validation catalog, hybrid integrators, coupling kernels, and checkpoint paths.
- At the initial audit, native gates used Julia 1.12.7 and 1.13.1. The latest
  gate has since advanced to 1.13.2 as recorded above. CI includes the minimum
  supported minor (1.12) and the moving latest stable version (`1`).
- Initial full suites were started in the live worktree on Julia 1.13.1 and
  1.12.7. Concurrent test/code edits mean these runs cannot prove the final
  revision; both were stopped through their own launch sessions after this was confirmed.
  Final full-suite checks use immutable checkouts. The initial 1.13.1 SHA
  compatibility finding is fixed; HDF5 and PencilFFT extensions have now
  precompiled successfully in both fresh version runs.
- Full package suite at fixed `1138317` on Julia 1.13.1: **106,662/106,662 tests
  passed**. Julia 1.12.7: **106,661 passed, 1 failed**; the sole failure was
  the JET report-count gate (31 findings against a ceiling of 30). The checkpoint
  inference repair passes focused JET checks on both versions without increasing
  the ceiling. Full-suite runs at immutable `c76e39c` passed **117,438/117,438
  tests on each of Julia 1.12.7 and 1.13.1** (48m56.7s and 37m37.6s respectively;
  both processes exited successfully). Later fixes have separate focused
  evidence and still require a final integrated run. Logs:
  `/tmp/hybrid-full-1138317-julia112.log` and `/tmp/hybrid-full-1138317-julia113.log`.
  New-run logs: `/tmp/hybrid-full-c76e39c-julia112.log` and
  `/tmp/hybrid-full-c76e39c-julia113.log`.
- Full package suite at immutable `67f01b7` on Julia 1.13.1 finished with
  **131,350 passed and 1 failed** in 81m34.8s. The collision rollback assertion
  compared uninitialized position storage containing `NaN` with `==`;
  `test/test_collisions.jl:653` and `src/particles/particle_set.jl:59` identify
  the fixture and allocation contract. The deterministic fixture repair passes all 512 collision assertions on
  each Julia version (`/tmp/hybrid-collision-fixture112.log` and
  `/tmp/hybrid-collision-fixture113.log`). The Julia 1.12.7 run passed all **131,351 assertions** in 95m08.0s and exited successfully. Logs:
  `/tmp/hybrid-full-67f01b7-julia112.log` and
  `/tmp/hybrid-full-67f01b7-julia113.log`. These snapshots include serial
  multispecies core/restart/diagnostics but exclude later AMR geometry,
  guiding-centre position, alias and MPI collection repairs.
- Julia 1.13.1 quick validation at fixed revision `b200c76`: **150/150 metrics
  passed across 20 cases**, no skipped results. Command:
  `julia +1.13.1 --project=validation validation/run_validation.jl --quick --no-plots`.
  Evidence: `/tmp/hybrid-validation-b200c76.log` and
  `/tmp/hybrid-validation-b200c76-results/artifacts/validation_summary.csv`.
- Julia 1.13.1 full electromagnetic suites at fixed `b62da7b`: **7,707/7,707
  assertions passed** across 40 test sets, including the new subcycling and
  Float32 cases. Log: `/tmp/empic-subcycling-b62da7b-julia113.log`.
- Julia 1.13.1 combined I/O, reproducibility and MPI checks at fixed `9856971`:
  **4,820/4,820 assertions passed** across 44 test sets. This includes the new
  two-species EM restart implementation. Log: `/tmp/hybrid-io-9856971-julia113.log`.
- CI run `38414172669` passed formatting and real MPI smoke at 2, 4, and 8 ranks.
  Its package suites were superseded by subsequent required pushes. Local full
  suites continue independently of push-triggered CI cancellation.
- Aqua passed all 9 quality checks on Julia 1.13.1 at `b62da7b`;
  the full JuliaFormatter 1.0.62 scan passed too.
- Baseline kernel harness: derivative 2.67 µs / 0 B; CIC deposit 305.29 µs for
  12,800 particles; hybrid step 1.54 ms for 3,200 particles. These are local
  measurements, not portable performance targets.
- Fused vector gather: all 304 deposition/gather assertions passed on
  Julia 1.12.7, including exact scalar equivalence and zero allocation. The repository benchmark
  measured 2.41–3.28× speedups against three scalar gathers with zero bytes allocated in
  every 1D–3D/NGP/CIC/TSC case. Hybrid and CAM-CL regression suites also passed.

## Audit coverage matrix

Each row requires domain/algorithm, implementation/error-path, resource, and
performance review. “Pending” means no current-run completion claim.

| Surface | Current work | Evidence needed to close |
| --- | --- | --- |
| Grids and spectral operators | Inventory complete; detailed pass pending | Manufactured solutions, boundary/Nyquist contracts, dependency tests |
| Particle storage/loading/pushers/boundaries | Coupling and step timing traced; remaining pass pending | Conservation, order, boundary maps, particle identity, invalid state |
| Deposition/gather/moments | Vector gather optimization verified; moment audit remains open | Adjointness, all shapes/dimensions/types, scalar equivalence, allocations |
| Electron closures and Ohm law | Integration traced; detailed pass pending | Tensor/force oracles, conservation budgets, dissipative terms |
| Hybrid/CAM-CL integration | Variable-step repair verified, including coupled-field convergence | Analytic timing, temporal order, subcycles, restart, MPI equivalence |
| Electrostatic/electromagnetic PIC | Timing, subcycling and restart verified; remaining audit open | Constant-force trajectory, Maxwell waves, charge continuity, mobile species |
| Hall-MHD, gyrokinetics, AMR | Hall-MHD domain/transaction/performance pass completed with no further confirmed defect; gyrokinetic position repair verified; AMR independent dyadic and finite-volume pass verified | Integrate latest Hall fixes into final full suites; reduce measured gyrokinetic costs and assess extreme BigFloat resource bounds |
| Collisions and radiation | Prior fixes exist; fresh pass pending | Conservation/rates, sampling, unequal weights, numerical robustness |
| Shocks/injection/campaigns | Inventory complete; detailed pass pending | RH laws, inflow flux, multidimensional/restart and published comparisons |
| Diagnostics and ray/Raycon solvers | Inventory complete; detailed pass pending | Independent references, units, transport/conversion invariants |
| Checkpoints/archives/field I/O | Independent audit active; defects confirmed | Configuration identity, exact bytes/checksums, atomic failure, restart |
| Threads/MPI/distributed FFT/GPU | Interfaces inspected; full pass pending | Real multirank tests, backend equivalence, ownership, hardware benchmarks |
| Packaging/quality/examples/validation | Full suites pass 117,438 assertions on each Julia version at `c76e39c`; later changes have focused checks | CI-version tests, Aqua/JET findings, formatting, examples, validation cases |

## Verified findings and fixes

| ID | Finding / acceptance criterion | State | Evidence / next action |
| --- | --- | --- | --- |
| PERF-001 | Vector gather repeats the same stencil three times | Pushed (`58fa769`) | 304 focused assertions plus hybrid/CAM-CL suites pass; reproducible speed/allocation comparison in `benchmark/benchmarks.jl` |
| IO-001 | `save_run` rejects its own valid dictionary output | Pushed (`d1db13d`) | Exact serialized bytes hashed under schema 4; 159 metadata/archive/field/async assertions pass, including cross-process and corruption cases |
| TIME-001 | ES/EM PIC accepts changed `dt` but uses the wrong half-step interval | ES pushed (`def6548`); EM pushed (`66ae6f3`) | ES suite passes including 180 new assertions; full EM suites pass 5,530 assertions; variable-step Maxwell/gyro convergence and charge/Gauss checks pass |
| IO-002 | Same-grid checkpoint silently accepts different model parameters | Pushed (`236d5d9`) | Restore saved model/shape parameters, reject incompatible concrete types; custom mutable configuration and alias tests pass |
| IO-003 | Truncated field file allocates declared payload before checking available bytes | Pushed (`c7b1c6f`) | 76 field/async/metadata assertions pass, including malformed extents, bounded allocation, empty and 0D arrays |
| IO-004 | HDF5 silently writes incorrect values from column-gapped views; range-backed arrays throw | Pushed (`b228b3f`) | 61 extension assertions pass, covering exact data/type/shape, contiguous-buffer identity and error recovery |
| TIME-002 | Hybrid/CAM-CL/MPI use incorrect kick interval after changing dt | Pushed (`236d5d9`) | 277 phase/convergence/MPI assertions; hybrid/CAM-CL/MPI regressions and real 2-rank suite pass; fixed-step comparison with original code is bitwise equal |
| IO-005 | Checkpoint restore fails on matching-length particle views and can partially mutate destination | Pushed (`236d5d9`) | 271 checkpoint/archive assertions pass, including 60 view/unsupported-resize assertions |
| PKG-001 | SHA compatibility excluded Julia 1.13’s bundled SHA 1.0 | Pushed (`4b5611d`) | Installed-version compatibility check plus 60 metadata assertions pass on each of Julia 1.12.7 and 1.13.1 |
| RESTART-001 | CAM-CL lacked checkpoint/archive methods and requires persistent magnetic state | Pushed (`1138317`) | 710 checkpoint/archive/MPI assertions pass, including 1D–3D, Float32/64, phase resets, edited fields, timestep/subcycle changes, malformed-state rejection, and exact continuation |
| CI-001 | Only Julia 1.12 was exercised; the test pipeline could mask failures | Pushed (`68b854d`) | YAML/matrix/artifact wiring checks pass; explicit CI Bash returns 1 for failed pipelines and 0 for successful ones; full version suites running separately |
| FORMAT-001 | Three files failed CI's JuliaFormatter 1 formatting rules | Pushed (`cf53bc4`) | Read-only scan with JuliaFormatter 1.0.62 identified exactly three Julia files; formatting-only repairs are idempotent and their diff was inspected. Contributor instructions now specify the CI formatter major |
| CFL-001 | Timestep recommendation lost representable results through intermediate overflow/underflow | Pushed (`27b4c8d`) | Valid 4-cell grid of length 1e-80 returned zero instead of the 512-bit oracle's 1.4184965709927286e-162. Logarithmic evaluation passes 3,802 helper/range/coupled-stability assertions on each of Julia 1.12.7 and 1.13.1 (`test/test_phase_helpers.jl`) |
| RESTART-002 | Electrostatic PIC lacked checkpoint/archive support | Pushed (`b200c76`) | Julia 1.13 combined IO/repro/MPI: 1,194 assertions; Julia 1.12 final focused suite: 540 assertions. Covers Float32/64, 1D–3D, exact variable-step continuation, changed fields, invalid destinations, views, empty populations, and archives |
| SUBCYCLE-001 | General EMPIC lacked electron subcycling | Pushed (`c5fef77`) | Julia 1.12.7 full EM suites: 7,257 assertions; added permanent Float32 coverage: 450 assertions in a separate final-file run. Includes independent constant-field oracles, 1D equivalence, 2D/3D charge/Gauss checks, changing dt/counts, mobile-ion coupled convergence. Default n_sub=1 arrays match prior code; allocations unchanged and scale linearly with subcycles |
| GYRO-001 | Guiding-centre drift loses finite results by reconstructing perpendicular speed squared from magnetic moment | Pushed (`2eeecd6`) | Direct moment-based grad-B drift; independent 2,048-bit equations cover extreme and ordinary Float32/64 inputs, charge signs, mass, zero moment, and failed-update rollback. Fusion suite: 256 assertions passed on each of Julia 1.12.7 and 1.13.1, including 76 new assertions; warmed push remains allocation-free |
| RESTART-003 | EM PIC lacked restart/archive support for two species and their staggered phases | Pushed (`9856971`) | Julia 1.12.7 combined I/O/reproducibility: 4,692 assertions, including 3,620 EM assertions. Covers legacy/general 1D–3D, Float32/64, all shapes, mobile/immobile and relativistic/nonrelativistic runs, exact continuation, diagnostics, changing dt/subcycles, invalid-state rollback, views and empty populations. Log: `/tmp/hybrid-io-empic-final112.log` |
| HALL-001 | Final Hall-MHD projection/Ohm evaluation can throw after changing public fields | Pushed (`f1f946e`) | Candidate fields are fully evaluated before copying to the public state. 114 Hall assertions pass on each of Julia 1.12.7 and 1.13.1; 108 independent before/after comparisons on 1.12 cover 1D–3D, Float32/64, three closures and projection on/off, with bitwise-equal successful steps and unchanged allocations. Failure tests compare every public field and preserve array identities. Logs: `/tmp/hall_transaction_final_tests.log`, `/tmp/hall_transaction_final_comparison.log`, `/tmp/hall-transaction-julia113.log` |
| TIME-003 | Finite solver steps overflow simulation time or wrap the integer step counter | Pushed (`f437600`) | Shared preflight covers hybrid, CAM-CL, Hall-MHD, both EM families and MPI. Clock and variable-step regressions pass 953 assertions on each of Julia 1.12.7 and 1.13.1; covers 1D–3D, Float32/64, moving particles, evolving fields, whole-state rollback, zero steps, Int boundary and EM substep rounding. Logs: `/tmp/hybrid-clock-final112.log`, `/tmp/hybrid-clock-final113.log` |
| JET-001 | Dynamic checkpoint field names widen inference to unrelated fields; Julia 1.12 full-suite JET gate fails | Pushed (`c76e39c`) | Typed field bundles and explicit ion-presence branches retain the existing schema. Frozen repair: JET reports 24 on Julia 1.12 and 22 on Julia 1.13, both below the unchanged ceiling of 30; no checkpoint reports. Concrete capture/validate/restore checks report zero across six families on both versions. I/O/repro: 4,692 assertions pass; eight serialized fixtures remain byte-identical. Logs: `/tmp/hybrid-io-typed-fields-*.log` |
| HALL-002 | Hall-MHD rejects stationary density exactly at the shared Ohm-law floor | Pushed (`d950dde`) | Inclusive bound matches the floor operation. Hall suite: 402 assertions pass on each of Julia 1.12.7 and 1.13.1; 288 new assertions cover exact/adjacent representable densities, three closures, 1D–3D, Float32/64, unchanged stationary fields and rejected-state rollback. Logs: `/tmp/hall_inclusive_floor_julia112.log`, `/tmp/hall_inclusive_floor_julia113.log` |
| HALL-003 | Final magnetic projection can produce an infinite electric field that Hall-MHD accepts | Pushed (`f96502c`) | Validate final E before committing candidate fields. Hall suite: 558 assertions pass on each of Julia 1.12.7 and 1.13.1; 156 new checks include unchanged rejected state, an unprojected retry and a finite projected neighbor. Independent comparison: 108 successful steps match prior fields/clock exactly, with unchanged allocations. Projected nonfinite B is already rejected by the spectral curl input guard (16 additional checks). Logs: `/tmp/hall_final_e_final_julia112.log`, `/tmp/hall_final_e_final_julia113.log`, `/tmp/hall_final_e_final_comparison.log`, `/tmp/hall_projected_b_guard.log` |
| RESTART-004 | Hall-MHD has no checkpoint/archive interface and carries inactive stage history | Pushed (`5c47d75`) | Add no-particle methods preserving all public fields and the retained stage `lapJ`. Focused Hall tests: 2,418 assertions on each Julia version; combined I/O/repro on 1.12: 7,110 assertions. Covers Float32/64, 1D–3D, three closures, initialization/live states, model changes, mutable closure ownership, array identity, rejection rollback and exact continuation. JET remains 24 reports on 1.12 and 22 on 1.13; Hall capture/validate/restore concrete analysis has zero reports. Logs: `/tmp/hybrid-io-hall-regressions112.log`, `/tmp/hybrid-io-hall-focused113.log` |
| GYRO-002 | Parallel velocity loses representable results by rounding field projection, force, acceleration and final addition separately | Pushed (`e9991e9`) | Compensated binary factors retain the complete expression; exact rational evaluation handles uncertain cancellation/range boundaries. Fusion: 550 assertions pass on each Julia version; 59,623 independent certificate/range/physical assertions pass on 1.12, including Float16 residual-underflow witnesses. Measured nonstationary and certified-equilibrium paths allocate 0 B; an uncertified mass=3 equilibrium still allocates 8,296 B. Local nonstationary latency increased about 20%; certified-equilibrium latency increased 36–44%. Logs: `/tmp/gyro_certified_final_112.log`, `/tmp/gyro_certified_final_latest.log`, `/tmp/gyro_certified_interleaved_performance.log` |
| PERF-002 | Hall-MHD repeats forward transforms for each velocity-gradient axis | Pushed (`870b544`) | Shared gradient in 2D/3D removes 12/24 forward FFTs per RK4 step. Hall tests: 627 assertions pass on each Julia version, including independent helical-wave, fourth-order convergence and nonlinear advection oracles. 108 normal-grid trajectory comparisons are exact; extreme anisotropic scaling can round differently. Benchmark `benchmark/hall_mhd_gradients.jl`: median baseline/current CPU ratios 1.085 at 128² and 1.133 at 64³; wall results vary, including mixed 2D repeats. Allocations increase 384 B/step in 2D/3D, unchanged in 1D. Logs: `/tmp/hall_gradient_final_julia112.log`, `/tmp/hall_gradient_final_julia113.log`, `/tmp/hall_gradient_final_comparison.log`, `/tmp/hall_gradient_repo_benchmark.log` |
| MSPEC-001 | Hybrid/CAM-CL accept one population despite existing multispecies moment deposition | Pushed (`b57d11e`) | Serial collection initialization/evolution preserves charge-weighted shared moments and independent q/m. Both Julia versions pass 2,286 integrated assertions, including 918 new collection checks. Independent uniform orbit, Boris force balance and wave-matrix convergence pass; 180 old/new scalar states agree exactly. Latest-Julia package JET remains 22 reports; all 24 concrete collection calls have zero reports; Aqua passes 9 checks. Example passes both versions with maximum velocity error 2.162e-8. Benchmark reports 18,992–35,456 B/two-ion step on Julia 1.13; alias checks account for 9,216–19,200 B and are a measured optimization target. Scalar extraction adds 16 B in 1D and 32 B in 2D/3D on Julia 1.12; latency samples vary. Logs: `/tmp/hybrid-multispecies-reviewed112.log`, `/tmp/hybrid-multispecies-reviewed113.log`, `/tmp/hybrid-multispecies-core-jet113.log`, `/tmp/hybrid-multispecies-benchmark113.log`, `/tmp/hybrid-multispecies-scalar-comparison.log` |
| AMR-001 | General scale division introduces nonzero curvature on exact affine profiles | Pushed (`d23dec5`) | Power-of-two normalization and a TwoSum rounding certificate preserve exact curvature; uncertain cases use the existing Rational fallback. BigFloat stored precision is retained, with directed rounding handled exactly. New tests: 7,223 assertions; existing AMR: 62, on each Julia version. Additional 166,013 independent oracle assertions pass on 1.12. Ordinary Float32/64 reductions allocate 0 B; measured local latency rises 19–26%. Logs: `/tmp/amr_curvature_final_julia112.log`, `/tmp/amr_curvature_final_julia113.log`, `/tmp/amr_curvature_final_resources.log`. Geometry and prolongation have separate confirmed findings; the AMR audit remains open. |
| RESTART-005 | Multispecies hybrid evolution requires ordered population restart and archives | Pushed (`ecb8c86`) | Shared capture/restore preserves all species and Hybrid/CAM-CL phase. Full I/O/repro on 1.12: 8,694 assertions; focused final-file collection suite on 1.13: 1,584 assertions, including 144 continuation configurations. Malformed later populations, shared empty storage and incompatible destinations reject before mutation; independent source/destination views and typed empties work. Final production JET remains 24/22 reports on 1.12/1.13, with zero reports in all 12 concrete scalar/collection I/O calls; eight scalar fixtures remain byte-identical. Logs: `/tmp/hybrid-io-collection-verified-regressions112.log`, `/tmp/hybrid-io-collection-handoff-focused113.log`, `/tmp/hybrid-io-collection-verified-jet112.log`, `/tmp/hybrid-io-collection-v4-jet113.log`. Hybrid clock and EM empty-alias repairs remain separate confirmed follow-ups. |
| MSPEC-002 | Collection energy/momentum budgets are absent; empty collection momentum fails | Pushed (`67f01b7`) | Extend canonical budgets and kinetic reduction; preserve heterogeneous nonempty momentum collections. Both Julia versions pass 1,098 diagnostics/regression assertions with independent 512-bit physical sums, all closures/inertia, unequal masses/weights, views and typed empties. Latest-Julia JET remains 22 package reports, all 24 concrete diagnostic calls have zero reports; Aqua passes 9 checks. Warmed energy/momentum particle reductions allocate 0 B for Float32/64 in 1D–3D. Logs: `/tmp/hybrid-multispecies-diagnostics112.log`, `/tmp/hybrid-multispecies-diagnostics113.log`, `/tmp/hybrid-multispecies-diagnostics-jet113.log`. |
| AMR-002 | Cell centers lose finite coordinates through intermediate overflow and rounded half-cell indices | Pushed (`06b637a`) | Complete-coordinate fused evaluation with exact half-index bounds and Rational fallback; Float16 widens to Float64 to avoid double rounding. Both Julia versions pass 9,388 AMR assertions, including 2,103 geometry checks. Additional 250,242 exact-oracle assertions pass on 1.12. Ordinary Float16/32/64 center calls allocate 0 B; local varying-index measurements show no regression. Logs: `/tmp/amr_center_final_julia112.log`, `/tmp/amr_center_final_julia113.log`, `/tmp/amr_center_resources.log`. Prolongation repair evidence follows. |
| AMR-003 | Premature slope rounding erases representable reconstructed children | Pushed (`c083043`) | Complete-child reconstruction shares the curvature rounding certificate; uncertified results use exact stored inputs. Julia 1.12.7 passes 710,102 assertions: 15,198 permanent/existing and 694,904 additional exact-oracle checks. Julia 1.13.1 passes the same 15,198 permanent/existing assertions. This includes both sides of the smallest-normal boundary. Ordinary measured profiles allocate 0 B; extreme cancellation allocates 4,096–8,664 B and boundary fallback 2,120–3,144 B. Logs: `/tmp/amr_prolong_boundary_final_julia112.log`, `/tmp/amr_prolong_boundary_final_julia113.log`, `/tmp/amr_prolong_boundary_final_cpu_resources.log`, `/tmp/amr_prolong_boundary_cold_resources.log`. Stored-precision restriction evidence follows. |
| AMR-004 | Restriction loses stored BigFloat bits before cancellation at lower working precision | Pushed (`ffb23ae`) | Complete dyadic averages retain both child operands. Each Julia version passes 15,468 permanent/existing AMR assertions; 1.12 additionally passes 186,072 independent exact averages and 101 public finite-volume/convergence/ownership/failure checks. Ordinary hardware restriction allocates 0 B; retained BigFloat samples allocate 3,064–12,008 B. Logs: `/tmp/amr_restrict_final_julia112.log`, `/tmp/amr_restrict_final_julia113.log`, `/tmp/amr_restrict_resources.log`, `/tmp/amr_final_public_audit.log`. Public finite-volume checks passed; a later dependency audit found AMR-005. Arbitrary BigFloat exponent-span resource use has not been validated. |
| AMR-005 | Exact-output conversion can round twice through ambient BigFloat precision | Pushed (`e5b082a`) | Shared integer quotient/remainder conversion rounds Float16/32/64 directly, independent of ambient BigFloat settings; BigFloat materialization preserves exact operands and caller rounding. AMR centers, curvature, prolongation and restriction use the shared owner. Corrected independent dyadic/neighbor oracles pass 1,196,587 assertions on Julia 1.12.7 and 65,268 on 1.13.1; 40 concrete JET calls per version have zero reports. The isolated 1.12 package analysis has 26 reports against the unchanged ceiling of 30 (`/tmp/shared_exact_jet_package112.log`); a count within the gate does not classify each report. Ordinary hardware AMR paths remain 0 B; extreme prolongation measures 4,096/4,768/10,272 B for Float16/32/64. Bounded stored-BigFloat samples measure 25,088–112,640 B; arbitrary exponent-span resources remain unverified. Logs: `/tmp/shared_exact_final_julia112.log`, `/tmp/shared_exact_final_julia113.log`, `/tmp/shared_exact_amr_jet112.log`, `/tmp/shared_exact_amr_jet113.log`. Gyro migration is a separate pending change. |
| IO-006 | Shared empty storage escapes alias checks and couples particle resize/write operations | Pushed (`9397d9b`) | Canonical alias checks trace data parents and mutable view indices for empty storage, retaining disjoint nonempty views. Final focused suites pass 1,752 assertions on each Julia version; broad 1.12 regressions pass 4,073 across 51 testsets. Original corruption repros pass 4/4. Package JET remains 24/22, with zero findings in all 48 concrete calls on each version. Two-population validation allocation falls from 9,216–19,200 B to 3,456–7,040 B across 1D–3D; cross-population helper remains 0 B. Logs: `/tmp/hybrid-empty-alias-v3-regressions112.log`, `/tmp/hybrid-empty-alias-handoff-focused113.log`, `/tmp/hybrid-empty-alias-v3-jet112.log`, `/tmp/hybrid-empty-alias-v3-jet113.log`, `/tmp/hybrid-empty-alias-v3-allocation.log`. Intra-population checkpoint aliasing and clock validation remain separate confirmed repairs. |
| MSPEC-003 | MPI hybrid APIs lack ordered ion collection evolution | Pushed (`ac88424`) | Canonical raw moment/push/predictor stages serve serial and MPI paths; collection preflight is collective. Final frozen v3: Julia 1.13.1 real 2-rank collection tests pass 48,049 assertions per rank; 1.12.7 real 8-rank suite passes, including 13,993 collection assertions per rank. Earlier v2 real 4-rank suite passed 25,329 collection assertions per rank; v2 integrated serial/MPI regressions pass 95,567. Scalar before/after comparisons preserve 180 serial and 72 MPI states exactly. Latest-Julia package JET at v2 has 20 reports and Aqua passes all 9 checks; final v3 concrete calls execute successfully with zero init/moment/nonmigrating reports and the same two Serialization reports as the scalar migration baseline. Density-floor narrowing and scratch alias regressions are covered. Two-rank benchmark smoke at 256 proton +128 alpha particles per rank measures 125,376 B maximum local allocation per two-ion step; this is not cluster scaling evidence. Logs: `/tmp/hybrid-multispecies-mpi-v3-2rank113.log`, `/tmp/hybrid-multispecies-mpi-v3-8rank112.log`, `/tmp/hybrid-multispecies-mpi-v2-integrated112.log`, `/tmp/hybrid-multispecies-mpi-v3-concrete-jet113.log`, `/tmp/hybrid-multispecies-mpi-v3-benchmark-2rank113.log`. Distributed collection restart is now verified in MPI-002. |
| MPI-001 | A destination validation failure on one rank leaves other checkpoint destinations restored | Pushed (`230fd58`) | VERIFIED: the baseline two-rank probe changes rank 0 while rank 1 rejects a non-resizable view (`/tmp/hybrid-mpi-restore-collective-probe.log`). Canonical restore preflight now finishes collectively before mutation. Final real two-rank tests pass on both Julia versions: 216 assertions on the failing rank and 258 on the healthy rank, including missing/corrupt files, incompatible state, unchanged destinations and successful continuation. Serial restore/repro/MPI regressions pass 2,342 assertions across 40 testsets on 1.13.1. Package JET has 24/22 reports on 1.12/1.13; all 30 concrete calls execute successfully on each. Latest Julia has zero concrete reports. Eight CAM-CL calls on 1.12 retain three reports each, exactly matching the unchanged baseline; no permanent gate was weakened. Logs: `/tmp/hybrid-mpi-restore-v2-2rank112.log`, `/tmp/hybrid-mpi-restore-v2-2rank113.log`, `/tmp/hybrid-mpi-restore-v2-serial113.log`, `/tmp/hybrid-mpi-restore-v2-jet112.log`, `/tmp/hybrid-mpi-restore-v2-jet113.log`, `/tmp/hybrid-mpi-restore-baseline-jet112.log`. The first test snapshot incorrectly compared copied Ref identities; final tests compare values. This preflight guarantee does not cover failures during subsequent allocation/copy. |
| GYRO-003 | Position loses finite updates by rounding drift, displacement and coordinate terms separately | Pushed (`109fe7a`) | Complete scaled displacement retains streaming, electric, magnetic-gradient and curvature terms through the final coordinate addition; exact fallbacks use the shared conversion owner. Each Julia version passes 2,422 fusion assertions; independent oracles add 157,888 on 1.12.7 and 245,045 on 1.13.1. Each version passes 60 concrete JET calls with zero reports and 76 runtime checks; matched package analysis retains 24 reports. All 16 measured ordinary hardware cases remain 0 B. Alternating local batches measured mixed-field costs 30–76% higher than `e9991e9`; equilibrium costs range from 21% lower to 17% higher. These shared-machine timings are not portable limits; further optimization remains open. Bounded retained-BigFloat pushes allocate 88,120–121,800 B; arbitrary exponent-span resources are unverified. Logs: `/tmp/gyro_shared_frozen_regressions112.log`, `/tmp/gyro_shared_frozen_regressions_latest.log`, `/tmp/gyro_shared_final_jet112.log`, `/tmp/gyro_shared_final_jet_latest.log`, `/tmp/gyro_shared_final_performance100k_settled.log`, `/tmp/gyro_shared_final_bigfloat_allocations.log`. Public drift and ring-averaging follow-ups remain open. |
| TEST-001 | Collision rollback test compares uninitialized positions containing NaN | Pushed (`44f749f`) | VERIFIED: the 1.13.1 full milestone log shows identical position vectors containing NaN failing `==`; the ParticleSet constructor allocates positions with `undef`. Initialize the fixture explicitly without changing the rejection or unchanged-state assertions. The full collision file passes 512 assertions across 32 testsets on each Julia version. Logs: `/tmp/hybrid-full-67f01b7-julia113.log`, `/tmp/hybrid-collision-fixture112.log`, `/tmp/hybrid-collision-fixture113.log`. A new integrated run is still required. |
| TEST-002 | Pressure resource gate equates allocator byte capacities with invariant allocation work | Pushed (`15662e1`); exact-SHA CI green | VERIFIED: `0d8c8e8` latest/minimum Julia CI records 14/15 failures only in two exact-byte equality assertions. Full allocation-call counters are unchanged in the displayed repeated samples. The replacement compares every allocation counter against compressed scalar work, retaining all numerical and zero-byte controls. [Evidence and contract](CI_RELIABILITY_PLAN.md#october-18-failed-package-checks-and-corrections). |
| TEST-003 | EMPIC subcycling fixture initializes weight but leaves position undefined | Pushed (`15662e1`); both full native files and exact-SHA CI pass | VERIFIED: `6e7fb01` latest-Julia CI records one setup error after 404,763 passes. Constructor and validator identify undefined position storage and the finite-position check. Initialize the one test particle at the origin; preserve all rejection/rollback assertions and test deliberately poisoned storage. [Evidence](CI_RELIABILITY_PLAN.md#october-18-failed-package-checks-and-corrections). |
| MOM-001 | Particle weights used as moment scratch are silently overwritten | Pushed (`d2c3da7`) | VERIFIED: baseline scratch aliasing changes weight 0.4 to 9.6 and bulk velocity (2, 3, 4) to (2, 6, 24). Shared output/scratch preflight now precedes mutation in serial and MPI moment APIs. Final focused tests pass 1,956 assertions on each Julia version; real two-rank tests pass 264 per rank on each. An earlier fixture error is recorded in the [repair plan](MOMENT_STORAGE_PLAN.md), alongside 2,145 passing neighboring assertions and the unchanged allocation gate. Package JET has 24/22 reports; all 36 concrete calls per version have zero reports and 72 analysis/runtime checks pass. Matched old/new results pass 288 checks exactly; all 18 warmed cases allocate 0 B, with paired local latency ratios 0.80–1.24. Logs and exact coverage are in the plan. |
| IO-007 | Overlapping components within one particle population corrupt checkpoint restoration | Pushed (`5af7a73`) | VERIFIED: the old EMPIC restore overwrites saved positions (0.2, 0.4) with velocities (0.6, 0.8); the final candidate rejects before mutation. A shared allocation-free component predicate now serves checkpoint and hybrid collection validation. Julia 1.12.7 broad IO/repro/collision/collection/diagnostic regressions pass 15,871 assertions across 101 testsets; Julia 1.13.1 final focused coverage passes 3,220, including every component pair, empty arrays, mutable view indices, later ion failures and exact view continuation. All 96 concrete JET calls per version have zero reports; package reports and leaf diagnostics match the 9397d9b baseline exactly (26/24). Warm validation measures 0 B across Float32/64, 1D–3D and concrete/erased storage types. Logs: `/tmp/hybrid-intra-alias-original-baseline113.log`, `/tmp/hybrid-intra-alias-original-fixed113.log`, `/tmp/hybrid-intra-alias-v3-regressions112.log`, `/tmp/hybrid-intra-alias-v3-focused113.log`. Shared clock validation and configuration-to-particle aliases remain separate follow-ups. |
| PIC-001 | The general EMPIC continuity diagnostic suppresses changes in total charge | Pushed (`43b6587`) | VERIFIED: 12 Float32/64, 1D–3D, odd/even baseline cases report zero for a uniform charge increment whose normalized residual is one (`/tmp/empic_dc_residual_probe112.log`). Retain true DC while excluding only derivative-null non-DC modes. Both Julia versions pass 3,923 assertions across all 10 complete affected caller testsets. An existing transverse-wave test was corrected: it retains the original non-DC relative bound and adds a density-scale absolute charge bound with an injected DC negative control; its frequency and Gauss checks are unchanged. Six concrete JET analyses have zero reports and 12 runtime assertions pass on 1.12. Warm allocations are unchanged in all six measured type/dimension cases. Logs: `/tmp/empic_dc_final_callers112.log`, `/tmp/empic_dc_final_callers113.log`, `/tmp/empic_dc_jet112.log`, `/tmp/empic_dc_resources_before.log`, `/tmp/empic_dc_resources_after.log`. |
| IO-008 | Hybrid checkpoint state accepts a nonfinite or negative clock that cannot advance | Pushed (`ba7da29`) | VERIFIED: 24 baseline clock cases restore successfully and then fail a positive step (`/tmp/hybrid-clock-baseline-probe.log`). The shared clock validator now rejects invalid time before save/restore. Each Julia version passes 1,590 focused assertions across nine family/layout variants, valid origins, exact continuation, archives and unchanged destinations. Julia 1.12 neighboring regressions pass 8,413 assertions across 46 testsets. Package JET reports 24/22 on 1.12/1.13 against the unchanged ceiling of 30; all 20 concrete clock checks per version have zero reports. Logs: `/tmp/hybrid-clock-v1-regressions112.log`, `/tmp/hybrid-clock-v1-jet112-native.log`. MPI save/load reach the shared validator; distributed invalid-clock save coordination is not covered by this slice. |
| GYRO-004 | Public drift helpers lose finite cross components through independent field normalization | Pushed (`3dd9f24`) | VERIFIED: tiny components disappear and a near-parallel cross differs by about 11.8% (`/tmp/gyro_public_drift_repro112.log`). Raw cross terms now retain their complete physical coefficient. Each Julia version passes 221,373 fusion, exact-range, physical and retained-precision assertions plus 165 mixed-input checks. All 88 concrete JET calls per version have zero reports and 104 runtime assertions pass; the fresh Julia 1.12 package gate has 24 reports against the unchanged ceiling of 30. Ordinary Float32/64 calls remain 0 B, but paired local normal-case latency is 1.90–4.25× the baseline; extreme fallbacks allocate 11,112–26,216 B. This correctness repair is not an optimization claim. Logs: `/tmp/gyro_public_drift_frozen_regressions112.log`, `/tmp/gyro_public_drift_frozen_regressions_latest.log`, `/tmp/gyro_public_drift_clean_jet_package112.log`, `/tmp/gyro_public_drift_clean_concrete112.log`, `/tmp/gyro_public_drift_compiler_resources_latest.log`, `/tmp/gyro_public_drift_performance112_run2.log`. Ring averaging and measured cost reduction remain open. |
| MPI-002 | Failed MPI save can destroy the previous rank file; collection restart is absent | Pushed (`aba03a9`) | VERIFIED: generation publication preserves the previous checkpoint on recoverable pre-publication failures; ordered collection restart shares serial validation and restore. Both-version real two-rank tests pass 360 continuation checks per rank and final 30 restore-failure +48 publication checks per rank. The complete two-rank Julia 1.12 suite passes. Schema 1 compatibility, external symlink ownership, cleanup error preservation and successful retries are covered. Warmed resource checks pass 120 assertions with no accumulating generations. Package JET is 26/25 against unchanged ceiling 30; all 24 concrete calls per version execute. [Plan and evidence](MPI_CHECKPOINT_PLAN.md) classify the reports and scope limitations. |
| PIC-002 | Finite light speed and finite curl lose representable Maxwell updates through intermediate range loss | Pushed (`0f53361`) | VERIFIED: the baseline zero-field vacuum acquires nonfinite fields when squared light speed overflows. The complete-expression repair passes 9,501 EM assertions on each Julia version; all 64 ordinary configurations retain exact state and allocation counts. Julia 1.12 package JET retains 24 reports; all 16 concrete helper/ParticleSet-ion calls have zero reports. Optional-Nothing findings exactly match the baseline. [Evidence, resource costs and follow-up plan](EMPIC_RANGE_PLAN.md). |
| PIC-003 | Raw spectral curl and half intervals lose representable complete Maxwell updates | Pushed (`8e7d5be`) | Scaled derivative bands and original timestep operands cover magnetic priming/centering/advancement and electric updates. Both Julia versions pass 19,765 complete EM assertions; all 54 concrete compiler signatures have zero reports, package count matches baseline 26. All 64 ordinary states and allocations remain exact; exceptional measured steps allocate 36,040–1,359,512 B. Shared-machine timings are mixed, including specialized 1D vacuum CPU increases of 11.6–18.2%. [Evidence and limits](EMPIC_RANGE_PLAN.md). |
| PIC-004 | Squared finite wavenumbers erase nonzero Poisson/Gauss projections | Pushed (`4f4fec4`) | Shared complete-expression projection preserves finite Fourier results and exact null classification. Both Julia versions pass 27,434 ES/EM assertions; 108 concrete JET signatures each are clean, package reports match the baseline, and 146 ordinary state/allocation comparisons pass. Local vacuum timing costs remain visible; FFT/density boundaries are separate PIC-005 work. [Scope, costs and evidence](EMPIC_RANGE_PLAN.md). |
| GYRO-005 | Ring phase, geometry and sample accumulation lose finite results | Pushed (`de6151d`) | VERIFIED: the Float16 constant-one mean becomes 0.5 at 8,192 samples and zero at 70,000; 12,000 points overflow the phase. Complete geometry and streaming certified accumulation pass 268,187 assertions on each Julia version. All 118 concrete JET calls per version have zero reports; the package gate retains 24. Ordinary calls remain 0 B but measured paired latency is 3.78–20× slower; optimization remains open. [Evidence and performance plan](GYRO_RING_PLAN.md). |
| PIC-005 | Density subtraction, FFT materialization and standalone residual arithmetic lose finite continuity results | Diagnostic range repair pushed (`1214234`); current correction planned | Both Julia versions pass 29,113 final EM/EM1D/IO assertions, 33 concrete compiler checks and bitwise ordinary/allocation comparisons. Package findings match baseline. Canonical input restrictions and severe arbitrary-grid cancellation accuracy remain explicit; current materialization is a separate repair. [Pipeline scope](EMPIC_RANGE_PLAN.md); [diagnostic contract, resources and remaining limits](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md). |
| PIC-006 | Legacy transverse Nyquist subtraction overflows or loses subnormal means | Pushed (`b73abfe`); GitHub green | Both runtimes reproduce 68 range witnesses. The final source adds complex compatibility controls and passes 1,090 focused/resource/concurrency assertions. Final combined quality passes 14,100/version, with complete concrete/package reports matching unchanged source. Final complete EM/EM1D/checkpoint integration passes 29,473/version. Exact-revision GitHub CI passes 369,711/369,712 assertions plus formatting and all MPI ranks. [Contract and costs](EMPIC_NYQUIST_RANGE_PLAN.md). |
| PIC-007 | Legacy longitudinal current narrows the first modal quotient before completing it | Pushed (`7b50b51`); all six exact-SHA CI jobs green | Both runtimes pass 37,916 complete EM assertions, six concrete compiler checks, unchanged complete compiler trees, 40 ordinary resource profiles and 16 exceptional profiles. [Scope and evidence](EMPIC_LEGACY_CURRENT_MODE_PLAN.md). |
| PERF-003 | Ordinary gyro ring points repeat complete expansion work | Pushed (`bc7e209`) | Cached lossless alignment and compensated phase products retain fallback and exact sample-mean behavior. Both Julia versions pass 269,015 regressions and 70,308 independent geometry checks; 124 concrete calls per version have zero reports, package count stays 24 on 1.12. Measured default 16-point time falls 18–37%, retaining 0 B; three-point/fallback results are mixed and compilation costs rise. [Evidence and limits](GYRO_RING_PLAN.md). |
| RR-001 | Radiation inputs can overflow during conversion and silently skip cooling | Pushed (`fdca823`) | VERIFIED: converted inputs now reject before mutation, including invalid no-op calls. Both Julia versions pass 1,810 forcing/particle/orbit assertions; all 27 concrete JET calls per version have zero reports, with package counts unchanged at 26/25. All 27 measured ordinary configurations preserve exact velocities and 0 B allocation (82 comparison checks); shared-validator error/tuple checks pass 54 on baseline and final code. [Evidence and remaining arithmetic plan](RADIATION_REACTION_PLAN.md). |
| RR-002 | Radiation arithmetic loses finite cooling through intermediate overflow/underflow | Serial gates pass; concurrent failure open | VERIFIED: 12 Float32/64 cases give unchanged or zero velocities instead of finite exponential-substep results; 48 independent high-precision witness checks pass on each version (`/tmp/hybrid-radiation-baseline112.log`, `/tmp/hybrid-radiation-baseline113.log`). [Plan](RADIATION_REACTION_PLAN.md). |
| RR-003 | Weak radiation drag can increase represented particle speed | Confirmed; candidate under verification | VERIFIED: Float16/32/64 stored-component squared norms increase under positive weak drag on both Julia versions (`/tmp/hybrid-radiation-weak-baseline112.log`, `/tmp/hybrid-radiation-weak-baseline113.log`). The arithmetic candidate uses a nonheating velocity factor; final gates remain open. [Plan](RADIATION_REACTION_PLAN.md). |
| REL-001 | Rounded velocity norms admit superluminal stored components | Pushed (`43881de`) | VERIFIED: shared preflight now compares exact stored components near the boundary. Both Julia versions pass 12,801 new boundary/rollback assertions and 1,810 neighboring checks; all 18 concrete compiler calls per version have zero reports, package counts remain 26/25. All 27 ordinary measured cases retain 0 B allocation; exact boundary cases allocate 3,168–3,248 B. [Plan and evidence](RELATIVISTIC_VALIDATION_PLAN.md). |
| NUM-001 | Concurrent MPFR shared-cache updates corrupt BigFloat transcendental results | Native cause verified; runtime integration open | Standalone C and Julia probes reproduce wrong pi with both installed runtimes. The final C sample fails in 1/8 and 6/8 stock processes; repaired-shared and pristine-unshared controls each pass 8/8. Both temporary configurations pass fresh concurrent radiation checks on both Julia versions. Installed runtimes remain unchanged. [Prepared upstream report](MPFR_UPSTREAM_REPORT.md), not submitted; [full evidence and remaining work](MPFR_CONCURRENCY_PLAN.md). |
| IO-009 | Configuration links to solver storage detach during restart | Pushed | VERIFIED: capture/restore reject retained model/shape links to evolving fields, particles, clocks, FFT plans, packed buffers and dictionary backing allocations. Both Julia versions pass 8,354 final ownership checks, 15,100 neighboring checks and real two-rank MPI integration (173/165 checks); all 16 concrete JET calls per version have zero reports. [Evidence and separate serialization gaps](CHECKPOINT_CONFIGURATION_PLAN.md). |
| IO-010 | Internal configuration views and native handles can detach during serialization | Pushed (`876baf7`) | Final source passes 12,950 focused assertions, 16 clean concrete compiler analyses and real two-rank MPI (296/340 checks) on each Julia version. Package reports match the earlier baseline. Preceding v4 passes 118,065 neighboring assertions/version; final shared-index extension is covered separately. All 18 standard configuration guards retain 0 B allocation. [Contract, exact snapshots, resources and limits](CHECKPOINT_CONFIGURATION_PLAN.md#io-010-final-verification). |
| IO-011 | Failed serial checkpoint serialization truncates the previous checkpoint | Pushed (`fe519e9`) | Atomic publication preserves prior files across serializer/write/flush/close/rename errors and reports cleanup failures. Both Julia versions pass 23,747 neighboring assertions and 293 focused checks; four concrete compiler paths are clean and package findings match baseline. [Publication contract and resource evidence](SERIAL_PUBLICATION_PLAN.md). |
| ENERGY-001 | EM field energy loses finite complete integrals through squared fields, c and cell volume | Pushed (`bef9802`) | Complete-integral cold evaluation passes 27,896 ES/EM assertions on each Julia version. All 24 concrete compiler calls per version are clean; package reports match baseline. Ordinary results stay exact and allocate 0 B; measured latency remains an optimization target. [Scope, resource investigation and evidence](ENERGY_RANGE_PLAN.md). |
| ENERGY-002 | Relativistic particle energy rejects valid extreme speeds or loses finite weighted results | Included in green `240f77c` | Each native version passes 2,138 permanent assertions, existing budgets, concrete/package compiler checks, Aqua and 398 four-thread checks. Ordinary calls retain exact baseline outputs and 0 B allocation; measured larger-case overhead remains 1.247–1.986× on 1.12 and 1.426–1.846× on 1.13. The three integrated source/test files match the verified snapshot. [Plan and remaining gates](ENERGY_RANGE_PLAN.md). |
| ENERGY-003 | Energy documentation equates Ohmic quadratic form with total loss and overstates reduced-inertia conservation | Documentation pushed (`2330379`); physical extensions planned | VERIFIED: source equations and 17 exact symbolic assertions give nonzero moving-ion balance terms. Existing tests cover fixed-ion inertia. Numerical APIs and forces are unchanged; a new pressure closure initially targets the ideal massless case. [Balance, witnesses and implementation boundaries](HYBRID_ENERGY_MODEL_PLAN.md). |
| ENERGY-004 | Nonrelativistic particle, hybrid magnetic and algebraic electron energies lose finite complete integrals | Kinetic/field/CGL/algebraic/inertia repairs green; complete-budget latency gates remain open | VERIFIED: 16 independent rational witness assertions and 20 additional collection/boundary assertions pass on each Julia 1.12.7 and 1.13.2 against source byte-identical to `e339350`. Includes false Inf, zero-field NaN, NaN from individually overflowing populations whose exact total is one, and finite-looking true overflow. The kinetic repair passes 12,943 numerical/resource, 2,235 integration/concrete-compiler/Aqua and 715 threaded assertions per runtime; full compiler trees match baseline. [Exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39102304110) is green in all six jobs: 386,955 assertions on actual Julia 1.12.7 and 386,956 on 1.13.2, formatting and all MPI ranks; measured certificate latency remains an optimization task. A further 72 witness assertions per runtime confirm polytropic power/ratio, isothermal pressure and CGL magnitude range losses before integration. The CGL integration now passes 15,743 numerical, 1,487 quality and 2,247 threaded assertions per runtime with unchanged complete compiler trees; its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39113547976) is green in all six jobs, with 389,341/389,342 assertions on actual Julia 1.12.7/1.13.2. The stored-pressure correction does not repair these legacy owners. [Scope and acceptance plan](HYBRID_DIAGNOSTIC_RANGE_PLAN.md). |
| ENERGY-005 | Finite rounded hybrid and electromagnetic field sums conceal true output overflow | Pushed (`5df5c38`); all six CI jobs green | VERIFIED: 16 independent rational boundary assertions pass on each Julia 1.12.7/1.13.2 against field source byte-identical to `013a70b`. Hybrid magnetic, EMPIC1D and EMPIC diagnostics return finite values below `floatmax` although the exact integral exceeds the overflow midpoint. The first shared repair fails two zero-allocation controls; the revised traversal passes 880 focused, 126 resource and 542 threaded checks per native runtime. Its solver/restart/concrete-compiler/Aqua suites pass 30,668 assertions per runtime and full compiler trees match baseline. Combined integration passes 13,907 numerical, 1,271 quality and 1,257 threaded assertions per runtime; complete compiler trees also match baseline. [Exact-SHA GitHub CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39110639668) is green in all six jobs: 387,505/387,506 assertions on actual Julia 1.12.7/1.13.2, package JET counts 23/22, formatting and all MPI ranks. [Evidence and scope](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#field-overflow-boundary--energy-005). |
| ENERGY-006 | Isothermal energy loses accuracy for large weak waves and uniform states | Pushed (`6e7fb01`); historical CI failed in TEST-003; traversal `5323c4d` green; corrected-source native gates pass | VERIFIED: on each Julia 1.12.7/1.13.2, three independent exact-moment/remainder assertions establish more than 16.08% error in the public Float32 diagnostic at 4,194,304 cells: 0.6159884 versus a complete-input enclosure near 0.7340845405043187. Source is byte-identical to `850fe46`. A separate corrected uniform-state witness passes 22 assertions per runtime and establishes 1.38%–3.97% Float32 errors outside the near-unit box; the range-only candidate certifies those inaccurate values. The bounded prototype passes 5,991 focused and 34 accuracy/resource assertions per native runtime, with zero allocations in all six measured cases; public integration remains open. [Evidence and implementation contract](ISOTHERMAL_ENERGY_ACCURACY_PLAN.md). |
| ENERGY-007 | Mixed-precision polytropic gamma selects the isothermal branch after narrowing | Pushed (`7104c61`); core, broad native and all six exact-SHA CI gates pass | VERIFIED: 40 witness assertions pass on each actual Julia 1.12.7 and 1.13.2 against the density-energy owner byte-identical to `63f2420`. Six unit-density 1D–3D cases return zero although the exact stored-gamma integrals are finite and nonzero; Float64-grid and exact-gamma controls pass. [Evidence and complete repair contract](POLYTROPIC_ENERGY_RANGE_PLAN.md). |
| ENERGY-008 | Shock bulk, particle and enthalpy fluxes lose finite complete values | Pushed (`daea194`); native/compiler and all six exact-SHA CI gates pass | VERIFIED: 52 baseline scalar/nodal witness assertions per runtime establish the range failures. Canonical V2 passes 17,169 permanent, 1,044 resource/cost and 462 shock/Aqua assertions per runtime, plus migration/format checks. Complete compiler trees match baseline (22/21). Earlier prototype and migration failures are retained. Ordinary latency overhead remains an optimization task. [Evidence and repair contract](SHOCK_ENERGY_FLUX_RANGE_PLAN.md). |
| ENERGY-009 | Complete shock total is lost by materializing components before combining | Pushed (`dd539da`); native/compiler and all six exact-SHA CI gates pass | VERIFIED: 36 witness assertions per runtime establish NaN for complete signed-one totals and zero for a two-subnormal total. The integrated repair passes 24,737 component/total assertions, 462 shock/Aqua and 414 cost assertions per runtime. Complete compiler trees match baseline (22/21); migration, formatting and four-thread checks pass. Ordinary latency and cold-storage costs are recorded. [Contract and evidence](SHOCK_ENERGY_FLUX_RANGE_PLAN.md#complete-total--energy-009). |
| ENERGY-010 | Ordinary Float32 magnetic sums lose percent-level accuracy on large grids | Pushed (`5c2796b`); all six exact-SHA CI jobs green | Both runtimes pass 21,313 combined assertions, all 164 cost profiles and complete compiler comparisons (30/27 reports). The final nine Julia files match the validated composition. [Contract and evidence](MAGNETIC_ENERGY_ACCURACY_PLAN.md). |
| GROWTH-001 | Append partially grows or corrupts particle storage on nonresizable/aliased inputs | Pushed (`16c9776`); all six exact-SHA jobs green | Both Julia versions pass revised focused/integration/compiler/resource gates; original baseline retains 48 append/ionization invariant failures. Cross-species backing-storage rejection is included. [Evidence and contract](PARTICLE_APPEND_PLAN.md). |
| GROWTH-002 | Migration changes aliased weights and partially publishes invalid replacement storage | Pushed (`7c6e1e4`); all six CI jobs green | V3 passes 8,495 focused and 99,191 combined assertions per runtime, all six real MPI 2/4/8-rank runs, complete compiler comparisons and resource/output gates. [Evidence and plan](MIGRATION_STORAGE_PLAN.md). |
| CLOSURE-001 | Polytropic pressure loses finite pointwise values and has large normal-output errors | Pushed (`0d8c8e8`); historical CI failed in TEST-002; corrected-source native gates and `15662e1` CI pass | VERIFIED: 44 witness assertions pass on each Julia 1.12.7/1.13.2 against `src/electrons/closures.jl` byte-identical to `63f2420`. Six cases return Inf or zero instead of finite nonzero pressures. Two finite cases return approximately 7.389 instead of a complete-input enclosure near 2.718. The public pressure owner feeds the Ohm-law pressure gradient. [Evidence and canonical repair plan](POLYTROPIC_ENERGY_RANGE_PLAN.md#pointwise-pressure--closure-001). |
| CLOSURE-002 | Concrete pressure path accepts union element types it cannot evaluate | Pushed (`15662e1`); native gates and exact-SHA CI pass | VERIFIED: 69 observations per Julia 1.12.7/1.13.2 establish 15 failures among 27 storage/coefficient combinations against `0d8c8e8`. The final correction passes 57,165 focused assertions per runtime and complete compiler trees match baseline (22/21). Solver, shock and checkpoint lanes also pass on both runtimes; all six exact-SHA GitHub jobs pass; complete CI report comparison removes only the two targeted diagnostics. [Plan and evidence](POLYTROPIC_ENERGY_RANGE_PLAN.md#algebraic-storage-dispatch--closure-002). |
| VALID-001 | Validation case 20 exceeds its charge-residual tolerance | Pushed (`e9c4fa7`); GitHub green | Final source passes the unchanged wave at 2.88e-12, 504 resource/oracle, 2,273 four-thread, 29,113 complete solver/restart and 13,268 numerical/compiler/quality assertions on each of 1.12.7/1.13.2. Exact GitHub CI passes 368,455/368,456 assertions plus formatting and all MPI ranks. Final capture intentionally changes ordinary bits; raw contracts remain tested. [Evidence, costs and production brief](CHARGE_DEPOSITION_CONSERVATION_PLAN.md). |
| PERF-004 | Legacy EM continuity diagnostic constructs FFT plans on every call | Pushed (`0f5b205`) | Reuses existing grid plans with 19,980 runtime assertions and 33 clean concrete compiler calls on each Julia version. Package findings exactly match baseline. All 24 measured configurations allocate 544 fewer bytes/call; timing variation is retained without a portable speed claim. [Evidence and compatibility](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md#legacy-diagnostic-plan-reuse--pre-edit-brief). |
| PERF-005 | Dyadic transform digit extraction constructs rationals per component/band | Pushed (`afc324f`); GitHub green | Each version passes 1,536 differential, 1,695 quality/state and 424 resource checks per snapshot. All outputs match; exceptional public diagnostics allocate 11.74–27.19% less across the two runtimes, with ordinary allocation totals unchanged. Full compiler report trees match baseline. [Contract, measurements and limits](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md#dyadic-digit-extraction-allocation--perf-005). |
| JET-002 | Julia 1.12 full-suite compiler-report gate exceeds its ceiling | Pushed (`6341c4b`); both native full suites pass | Reuse the validated archive byte-vector binding through hashing and decoding. Final focused compiler gates pass at 22/21 package reports on native 1.12/1.13, with zero concrete String/SubString loader reports. Exactly the nine loader findings disappear from the 1.12 baseline; all other leaves remain. Each version passes 234 neighboring assertions and 20 additional metadata checks; frozen full suites now pass 352,287 assertions per version. [Cause, compatibility and evidence](ARCHIVE_INFERENCE_PLAN.md). |
| CI-002 | Fresh latest-Julia resolution fails and MPI publication fixture races its readers | Pushed; both package jobs and 8-rank publication checks pass | Current-registry resolution requires JET 0.12 on Julia 1.13. The delayed 8-rank replay reproduces the missing-file failure; a pre-rename barrier passes 48 publication assertions on every rank on both native Julia versions. JET 0.12.3 passes all 3 analysis assertions (22/21 reports), 9 Aqua checks and 80 archive assertions per version; the ceiling remains 30. [Plan and GitHub status](CI_RELIABILITY_PLAN.md). |
| KDV-002 | Quadratic filtering allows excluded-input and cutoff aliases | Pushed (`c5d7e94`); both GitHub package jobs pass | VERIFIED: explicit cosine inputs change by 0.00150–0.00175 despite zero projected nonlinear RHS; Nyquist cross products also leak. Strict integer-mode operand/product projection passes 61 KdV checks, two concrete compiler analyses and eight downstream metrics on both versions. New mode tests fail 14 checks against the baseline; physical accuracy thresholds remain unchanged. [Evidence and plan](CI_RELIABILITY_PLAN.md#kdv-filtering-repair--pre-edit-brief). |
| RAY-001 | Magnetic finite-difference derivatives amplify unit-roundtrip perturbations | Pushed (`6a6b316`); both GitHub package jobs pass | VERIFIED: one-ulp launch changes perturb focusing by up to 4.87e-7 with the stencil versus 3.76e-13 with analytic derivatives. Both versions pass 961 complete-suite assertions and the final 10-assertion rounding block (two overlap). Existing normalized/SI and MATLAB tolerances remain unchanged; both Linux package jobs at `55d996b` now pass. [CI status](CI_RELIABILITY_PLAN.md). |
| THR-007 | Threaded deposition overwrites aliased values instead of rejecting them | Included in green `240f77c` | VERIFIED: one- and four-thread probes clear aliased inputs while serial deposition rejects them unchanged. Shared preflight passes 3,012 one-thread integration and 754 four-thread assertions per native version; 24 concrete analyses are clean, JET/Aqua pass. New tests fail 324 checks on the baseline. All 18 measured ordinary cases retain exact outputs and allocations. [Plan](THREADED_DEPOSITION_PLAN.md). |
| CI-003 | Eight-rank MPI job is interrupted by its 20-minute step limit | Pushed; complete GitHub run passes at `240f77c` | VERIFIED: the previous explicit deadline is replaced by a bounded 40-minute eight-rank window, retaining every test and the 20-minute 2/4-rank limits. Run `38741178196` passes all three MPI jobs; its eight-rank smoke step completes in 19m53s. Both package jobs also pass, with 357,258/357,259 assertions on actual 1.12.7/1.13.2. [Plan](CI_RELIABILITY_PLAN.md#eight-rank-execution-budget--ci-003). |
| MOMENT-001 | Per-species deposition, normalization and charge scaling lose finite first moments | Included in green `240f77c` | Final arithmetic passes 1,953 permanent assertions on both native versions; the final worker-state guard has separate preflight/threaded controls on 1.12 and the full permanent run on 1.13. Concrete and fixture-loaded package JET, Aqua and ordinary baseline/resource comparisons pass. Combined species/background expressions remain separate. [Contract, gates and limits](MOMENT_NORMALIZATION_PLAN.md). |
| MOMENT-002 | EMPIC species/background combination loses finite total charge | Pushed (`4cd642c`); complete GitHub run green | Both package jobs pass on actual Julia 1.12.7/1.13.2, alongside formatting and all MPI ranks. Native focused charge, solver, concurrent, neighboring, compiler and resource gates pass; all 16 ordinary cases retain bitwise results and 0 B. [Plan and acceptance gates](EMPIC_CHARGE_COMBINATION_PLAN.md). |
| MOMENT-003 | EMPIC current combination loses complete values; deposition errors retain temporary positions | Pushed (`6ad7276`); complete GitHub run green | Both native versions pass 4,805 checks on the final combined charge/current source, including 892 current assertions, public steps, 70 concrete analyses and unchanged package JET/Aqua. Each also passes 1,862 existing integration assertions. Ten ordinary cases preserve exact outputs/scratch and 0 B allocation. [Contract and gates](EMPIC_CURRENT_COMBINATION_PLAN.md). |
| CI-004 | Required Aqua import errors become skipped-success quality gates | Included in green `240f77c` | VERIFIED: an actual missing-Aqua baseline exits 0 with one broken test. The repaired gate passes all nine checks on both native versions; missing-package and injected import-failure probes now exit 1. Existing checks/options are unchanged. [Contract and gates](CI_RELIABILITY_PLAN.md#required-aqua-gate--ci-004). |
| CI-005 | New pushes cancel unfinished checks for earlier commits | Pushed (`4c03c87`); complete GitHub run green | The group now includes commit SHA; actionlint and parsed workflow/key checks pass. Older runs retained their package jobs across subsequent pushes, and `240f77c`/`4cd642c` have since finished green. All jobs, thresholds and matrices are unchanged. [Plan](CI_RELIABILITY_PLAN.md#preserve-acceptance-for-every-revision--ci-005). |
| PERF-006 | Rank classification allocates indexed labels and pencil orientation intermediates | Pushed (`b3e2eb6`); all six GitHub jobs green | V2 passes 46,188 focused checks per runtime and reduces measured complete-migration allocation by 59.68–84.19% / 47.04–76.16%, with identical outputs. Complete compiler trees match baseline locally and in CI; actual GitHub Julia versions are 1.12.7/1.13.2. [Contract](RANK_CLASSIFICATION_ALLOCATION_PLAN.md). |
| COLL-001 | Ionization loses threshold, cooling and probability accuracy at floating-point boundaries | Pushed (`b251cb7`); all six CI jobs green | Final native focused, endpoint, integration, compiler and ordinary-output gates pass on both runtimes. The local full run retains its separately repaired metadata harness error. [Contract and evidence](IONIZATION_ARITHMETIC_PLAN.md). |

## Completeness roadmap

These are capability gaps or assessment tasks, not all classified as bugs.
Inspect implementation contracts before choosing the design; add an acceptance
benchmark before implementation. The roadmap remains open to new findings.

| Priority | Capability / gap | Plan and completion evidence | State |
| --- | --- | --- | --- |
| 1 | Reliable variable-step evolution | Correct phase state in ES/EM, hybrid/CAM-CL and MPI; preserve it in checkpoints; verify analytic trajectories and convergence | Active |
| 1 | Restart across production solvers | Hybrid, CAM-CL, electrostatic/electromagnetic PIC and Hall-MHD share capture/validation/restore infrastructure; exact continuation includes EM species and Hall retained stage history | Existing solver families and serial hybrid collections verified; boundary validation follow-ups active |
| 1 | Reproducible archives and bounded field I/O | Fix checksum representation and malformed-file allocation; cross-process and corruption tests | Active |
| 2 | Multispecies hybrid evolution | Collection initialization/stepping absence reproduced; [implementation plan](MULTISPECIES_HYBRID_PLAN.md) defines shared charge moments, independent species push, restart and conservation/dispersion acceptance | Serial core, collection restart and conservation diagnostics verified; MPI collection evolution and distributed collection restart verified |
| 2 | Multispecies conservation diagnostics | Canonical reductions and budgets support typed collections and empty populations; independent particle/field/closure sums and scalar regressions pass on both Julia versions | Verified |
| 2 | Local moment accuracy under strong weighted cancellation | A 48-check Float32/64 witness retains raw non-DC error despite a correct reconciled total. Compare compensated/binned accumulation with exact spatial oracles; define ordinary rounding compatibility and require threaded, coupled and measured cost gates. [Scope](CHARGE_DEPOSITION_CONSERVATION_PLAN.md#remaining-local-accuracy-boundary) | Assessment |
| 2 | Exactly conservative coupled hybrid integration | Existing explicit predictor has second-order momentum error verified against a discrete Boris oracle. Assess a coupled midpoint field/current solve with convergence, conservation and nonlinear solve-residual acceptance | Assessment |
| 2 | Evolving electron pressure integrated into steppers | Ideal massless closure includes joint B/pressure stages, rollback, MPI, persistence and explicit-pressure budgets; both native Julia gates pass. [Contract and implementation checklist](HYBRID_ENERGY_MODEL_PLAN.md#evolving-pressure-implementation-tracker) | Implemented; all six exact-revision CI jobs green (`013a70b`) |
| 2 | Multidimensional EM electron subcycling | Common EM path extended; charge-continuity, Gauss-law, species phase, coupled convergence and Float32/64 checks pass | Verified |
| 2 | Nonperiodic EM-PIC boundaries | README declares periodic EM-PIC only; design compatible field/particle boundaries with reflection/transmission and charge tests | Planned |
| 3 | Production distributed stepping | MPI reference stepping replicates fields; design domain-owned fields/moments with migration and distributed operators; verify serial equivalence and scaling | Planned |
| 3 | Complete GPU-resident stepping | Metal particle kernels exist and uniform push now passes actual-device gates; hybrid/EM field storage remains CPU Array-based and CUDA currently supplies storage/telemetry. Generalize the canonical field/step owners, then verify complete backend trajectories and measured resources | Partially implemented |
| 2 | Spatially local, interspecies Coulomb collisions | Current public Coulomb operator takes one ParticleSet; isolated cell-local collection implementation is undergoing independent rate, normalization and rollback checks | Active; unaccepted |
| 2 | Integrated implicit field advancement | `run_whistler` evolves a linear reference packet; the shock CN field option uses four fixed Picard iterations without a residual check. Assess convergence-controlled field/current/particle coupling, stability, conservation and rejection benchmarks | Partial implementations; assessment |
| 3 | Model/geometry completeness | Assess Hall-MHD, gyrokinetic, AMR and fusion building blocks for missing end-to-end coupling before planning extensions | Assessment |
| 2 | Guiding-centre position, public drift and ring cost | Position, drift and ring correctness repairs verified. Ordinary hardware paths remain allocation-free; paired local latency increases are recorded in GYRO-003/004/005. Profile complete point construction, streaming accumulation and cross-expansion work; require unchanged exact-range certificates plus measured improvement | Correctness and ordinary-ring optimization verified; further costs remain open |
| 3 | Published validation and external comparisons | Run local validation; inspect full Preisser profile replay and external solver cases; add reproducible accuracy/cost comparisons | Planned |
| 4 | Independent final audit | Revisit every matrix row from another entry point, resolve all confirmed defects, rerun final acceptance gates | Pending |

## Capability implementation queue

The [expanded implementation plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)
records checked source, canonical owners, dependencies, model contracts and
acceptance gates for every slice below. This tracker owns accepted fixes and CI;
the expanded plan owns capability staging. Both remain active.

| Capability | Owner / current state | Next acceptance milestone |
| --- | --- | --- |
| CAP-001 Local multispecies collisions | Collision agent; isolated implementation, unaccepted | Independent rates, units, locality, rollback, legacy compatibility and resource gates |
| CAP-002 Nonperiodic EM-PIC | Unassigned; planned | Compatible bounded Maxwell/current operators and vacuum/wall oracles |
| CAP-003 Domain-owned distributed stepping | Unassigned; replicated reference exists | Local storage, halos and distributed operator equivalence on real ranks |
| CAP-004 Complete GPU-resident stepping | Root auditing primitives; coupled solver planned | Gathered-push validation, remaining kernels, then one resident solver with measured transfers |
| CAP-005 Coupled implicit advancement | Architecture agent auditing existing shock CN; general coupled solver planned | Independent trapezoidal solve, residual/stability/failure contract |
| CAP-006 Geometry, gyrokinetics and adaptive evolution | Unassigned; standalone primitives exist | Separate coupled model contracts and manufactured solutions for each branch |
| CAP-007 External validation and comparisons | Unassigned; partial reference comparisons exist | Pinned provenance, production time-series dispersion, full profiles and matched cost/accuracy |
| CAP-008 Arbitrary ES/EM species | Expanded scope; planned | Species/background/time-level contract and collection ES slice |
| CAP-009 Full-Hall multidimensional shocks | Expanded scope; planned | Full induction/boundary balance and transverse-mode oracle |
| CAP-010 Composed sources and RNG restart | Expanded scope; required for restartable collision integration | Split phases, complete RNG persistence and transaction semantics |
| CAP-011 Electron transport and energy | Expanded scope; planned | Term-specific heating/loss balance and transport convergence |

Every slice requires physics/conservation, convergence, failure/rollback,
restart, compiler/quality, applicable hardware/parallel and fixed performance
evidence, then an independent audit and exact-SHA CI. Unsupported combinations
remain explicit. Existing green CPU CI does not close missing physics, GPU
residency or cluster-scaling gates. Held radiation changes remain separately
owned and excluded from current promotions.

## Verification commands

- Full package/quality suite in immutable checkouts, on both Julia 1.12 and latest
  stable: `julia +1.12 --project=. -e 'using Pkg; Pkg.test()'` and
  `julia +1.13.2 --project=. -e 'using Pkg; Pkg.test()'`.
- Coupling: `julia +1.12 --project=. test/test_deposit.jl`.
- Kernel timing/allocation: `julia +1.12 --project=. benchmark/benchmarks.jl`.
- Validation: `julia +1.12 --project=validation validation/run_validation.jl --quick`
  followed by the complete case selection after inspecting runtime/dependencies.
- MPI: CI invokes `test/test_mpi_multirank.jl` at 2, 4 and 8 ranks using
  `MPI.mpiexec()`; record each observed result.
- Final diffs: `git diff --check`; format changed Julia files with the configured
  JuliaFormatter version; inspect the exact staged files before each commit.

No full-codebase, complete-physics, maximum-performance, GPU-production, or
cluster-scaling claim has been established by this run yet.
