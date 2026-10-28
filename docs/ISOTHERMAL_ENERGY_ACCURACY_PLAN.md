# Isothermal energy accuracy — ENERGY-006

Revision-scoped implementation record. The baseline failures and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

The numerical repair landed as `6e7fb01`. Its [exact-SHA run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129297310)
finished with one failed latest-Julia package job; the minimum package, formatter
and all MPI jobs passed. The separate traversal revision `5323c4d`
[passes all six jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129564503).
These outcomes were checked again 2026-10-29. The pre-promotion gates and
measurements below are historical; the earlier failure stays failed.

Recorded pre-promotion status: defects confirmed on Julia 1.12.7 and 1.13.2. Public v2 is frozen and
under verification. Public v3 separately integrates verified traversal optimizations;
no public isothermal repair has been pushed. Public v1
remains unaccepted because a fresh resource witness found excessive allocation.

## Confirmed defects and regression contract

| Case | VERIFIED observation | Evidence |
| --- | --- | --- |
| Large weak Float32 density wave | At 4,194,304 cells, amplitude 0.001, temperature `0.7f0` and unit spacing, the public result is `0.6159884f0`. The independent complete-integral enclosure is near `0.7340845405043187`, width `4.374e-17`: more than 16.08% relative error on both runtimes. | `/tmp/hybrid-isothermal-range-design-v1/large-accuracy-witness-results.json` |
| Large uniform Float32 states | At the same count/temperature/spacing, densities 0.5, 0.75 and 1.5 have errors of at least 2.6357%, 1.3758% and 3.9651%. The range-only prototype certifies all three inaccurate values. Each corrected native witness passes 22 assertions. | `/tmp/hybrid-isothermal-uniform-accuracy-20261018/corrected-evidence.json` |
| Intermediate range loss | Pressure, pressure-times-logarithm and spacing products can lose range before a normal complete result. The original range probe incorrectly accepts 19 cases, established by 71 assertions per runtime. | `/tmp/hybrid-isothermal-range-design-v1/rejection-and-accuracy-evidence.json` |
| General v1 resource regression | Four normal temperature-weighted integrals at 256 cells take the cold path. Complete calls allocate 388,632–16,429,056 B versus 1,120/2,112 B in the existing public routine. Both native witnesses pass 20 assertions against independent logarithm bounds. No large-grid cold allocation was measured or extrapolated. | `/tmp/hybrid-isothermal-logrange-resource-20261018/rejection-evidence.json` |

The weak-wave oracle sums six dyadic density-deviation moments exactly in Int128
and bounds the full series remainder. Its input/count bounds prove the integer
accumulators fit. A separate 56-assertion cross-check passes on 1.13.2 against
BigInt sums and independent alternating-logarithm bounds. Resource controls use
an eighth-moment enclosure, separately cross-checked in 144 assertions. These
oracles do not use floating logarithms or MPFR transcendental results as proof.
These witnesses use the public density-energy owner byte-identical to `63f2420`.

A normal-range certificate alone cannot repair these accuracy defects. Old bits
are not an acceptance criterion for finite isothermal reductions. Nonfinite and
unsupported-type behavior must retain the existing contract. Selected ordinary
Float64 controls allow two output ulps against independent complete integrals;
the bounded near-unit and cold paths require unambiguous final rounding.

## Current implementation and verification

The isolated implementation has three paths:

1. Finite states whose positive densities lie within `1/256` of one use exact
   linear sums, a bounded fixed-point entropy polynomial and its full analytic
   remainder. Bounded dyadic arithmetic applies temperature and each spacing;
   a result returns only when both endpoints round identically. Unresolved
   states go directly to complete logarithm evaluation.
2. Other finite states use Float64 terms and the existing compensated-addition
   primitive. Independent logarithm bounds and normal-operation checks cover
   pressure, terms, accumulation, volume and output range. This path does not
   promise universal single rounding.
3. Unresolved finite cases retain complete rational factors and bounded
   logarithms without MPFR state. UInt64 logarithm keys preserve the exact
   odd significands, which have at most 53 bits. Actual nonfinite inputs use
   the existing arithmetic.

| Component | Completed native evidence | Remaining status |
| --- | --- | --- |
| Near-unit v4 | 5,991 focused, 34 accuracy/resource, 24 clean concrete compiler and 1,226 parent-threaded assertions per runtime; all six measured calls allocate 0 B. | Tested component; public acceptance pending. |
| Combined UInt64-key v2 | Both seven-step lanes pass: 243 cold, 6,243 key-bound, 660 range, 552 near-unit, 34 resource, 24 concrete compiler and 1,226 parent-threaded checks. Key tests cover 1,037 densities. | Fixes the compiler representation issue; its old ordinary reduction does not repair uniform-state accuracy. |
| General reduction v1 | 1,416 focused, 72 uniform accuracy/resource and 1,627 quality/threaded assertions per runtime, including 24 clean concrete compiler and 17 independent cold rounding-cell checks. Eighteen ordinary measured calls allocate 0 B. | Unaccepted after the separate extreme-density allocation witness. |
| General reduction v2 | Both versions pass 620 enclosure, 1,416 focused, 72 ordinary-resource, 32 extreme-resource and 1,627 quality/threaded assertions. All 18 ordinary and eight extreme calls allocate 0 B; large controls use 4,194,304 cells. Other value arithmetic and guards are byte-identical to v1. | All isolated gates pass; public integration is running. |
| Public v1 | Both versions pass 8,532 permanent assertions plus one preflight check and 72 concrete compiler checks. Separate permanent key tests pass 6,243 assertions/version. All measured weak-wave and uniform public calls allocate 0 B. | Solver/restart, quality, threaded and complete compiler-tree lanes continue; resource regression prevents acceptance. |
| Public v2 | Both versions pass 15,427 permanent assertions plus one preflight check each. Both versions also pass 72 concrete compiler checks each. There are 27 helper-equivalence checks and executable-AST equality across 205 source/test files. All measured public calls allocate 0 B. The frozen 294-file snapshot has a checked 15-file Julia delta against `63f2420`; other integrator arithmetic is byte-identical. | Pushed as `6e7fb01`; broader native lanes and exact-SHA CI remain pending. |

The separate v5 components complete their isolated native gates on both versions:

- Near-unit v5 retains all 5,991 numerical, 34 resource, 24 concrete compiler and
  1,226 parent-threaded checks. Its extra 98 assertions per version cover eight
  certified and eight unresolved degree-eight first passes against independent
  integral bounds. Unresolved Float64 cases retain the degree-12 refinement.
- General v5 retains all 2,140 numerical/resource and 1,627 quality/threaded
  checks, including 24 concrete compiler calls. Each runtime also passes 15,360
  bit-for-bit complete-probe comparisons. Per-density reuse preserves every
  accumulation and input check; zero temperature skips the reuse-selection scan.

Their `all-native-gates-evidence.json` files are in
`/tmp/hybrid-isothermal-nearunit-v5` and
`/tmp/hybrid-isothermal-general-reduction-v5`. Public v3 combines these changes
in the existing isothermal helper, adds two permanent regression files and
updates their test includes. All other production files match public v2. Its
296-file Julia snapshot passes 29 parent-helper comparisons and executable-AST
comparison across all 207 source/test files. Both native public focused runs
pass 27,831 regression assertions plus one preflight check each, and 144
concrete compiler checks each. Broader lanes await their public-v2 parent lanes;
the optimization is pushed as `5323c4d` with its own GitHub run concurrently.

General v2 uses exact frexp range reduction outside `[1/2,2]`, combining bounds
for the mantissa logarithm and an integer multiple of `log(2)`. The real extreme-
density resource controls include 4,194,304 cells. The public owner preserves
closure validation order, B-adapter validation and generic/polytropic dispatch;
the bounded dyadic conversion joins the existing exact-conversion owner.

A 26-check public-v1 helper comparison verifies retained arithmetic bodies and
the statistics-entry split. Its documentation overlay preserves executable ASTs
across 203 source/test files. Public v2 passes its fresh helper and formatting checks; its complete compiler
comparisons remain pending. Its `6e7fb01` promotion excludes the held radiation
files and preserves all already-pushed kinetic, field and CGL changes.

Evidence directories:

- `/tmp/hybrid-isothermal-nearunit-v4`
- `/tmp/hybrid-isothermal-integrated-design-v2`
- `/tmp/hybrid-isothermal-general-reduction-v{1,2}`
- `/tmp/hybrid-isothermal-public-v{1,2}-verification`

## Retained failures and measured costs

| Earlier experiment | Outcome retained |
| --- | --- |
| Range probe v1 | Both native large-wave resource gates fail: 40 passes/one failure. The interval spans approximately -162.46 to 208.69 around the inaccurate old value. The large cold fallback was not executed. |
| Near-unit v1 | 1,902 focused and 28 resource/accuracy checks per runtime. Kernel allocation is 0 B; whole calls allocate 9,896–10,440 B. Large shared-host samples take 3.47–4.86 s. |
| Near-unit v2 | 5,823 focused checks per runtime. The minimum-version resource run records 0 B in all six cases but fails two inherited `0 < 0` comparisons; the latest lane does not run afterward. |
| Near-unit v3 | The corrected resource contract requires zero bytes directly. Both native focused/resource/compiler/threaded gates pass; the large Float32 paired median is 1.274 s versus 3.198 s for v2. |
| Near-unit v4 timing | Large Float32 paired median: 0.721 s versus 1.047 s for v3, ratio 0.689. The small Float32 ratio is 2.40; Float64 ratios are 1.04–1.06. No universal speedup or optimality claim follows from these shared-host samples. |
| Range guard v2 harness | The first run has 607 passes/five failures from requiring exact rounding outside its then-written ordinary-compatibility scope. Classification verifies unchanged public bits. A separate corrected harness passes 660 checks/version; the later uniform witness broadens the accuracy contract. |
| Combined isothermal v1 | Minimum Julia reports an undefined-variable possibility in Base BigInt hashing through the cold dictionary; latest Julia passes. UInt64 keys fix the representation without filtering reports or claiming a Julia runtime bug. |
| Uniform witness harness | The initial endpoint-distance expression only applies outside the oracle interval. A separate corrected witness uses general interval distance; only corrected runs support the finding above. |
| Helper-comparison harness | The first comparison has 24 passes/one failure from an `@assert` source-location node. Explicit source-location normalization yields 26 passes without changing production arithmetic. |

## Performance investigation

The paired public-v2 comparison uses the exact `63f2420` public owner, including
validation, with only its function name changed. Warmed calls alternate order on
the same data. All 16 cases pass 176 accuracy/stability/resource assertions on Julia 1.13.2;
all 48 measured corrected calls allocate 0 B. Paired median latency ratios are
4.16–16.72 relative to the old routine on this shared host. These costs remain
an optimization task.
Evidence: `/tmp/hybrid-isothermal-public-v2-timing`.

An isolated general-v3 experiment reuses the last active density's term and
term enclosure. Every cell still executes compensated accumulation, input
checks and outward sum bounds. Both native versions pass 9,216 bit-for-bit full-probe comparisons. A separate
latest-version paired experiment passes 132 output/resource assertions. Uniform
and repeated-block cases improve, while varied-input sample ratios span
1.002–1.304; the latter require follow-up before adoption. Public v2 remains frozen;
no v3 source is integrated or accepted. Its nine-sample large varied-input
follow-up passes 58 assertions and reports ratios 1.099 (Float32) and 1.030
(Float64). Scope and experiment: `/tmp/hybrid-isothermal-general-reduction-v3`.

General v4 selects reuse when every active density run has at least two
members, and otherwise uses the uncached specialization of the same loop.
Both native versions pass 15,360 full-state comparisons, including forced
checks of both specializations. Its 132-assertion latest paired run has large
uniform/block ratios 0.170–0.232 and large varied ratios 1.021/0.951. Small varied
samples have ratios 2.755/1.920; a 21-pair repeat passes 130 assertions and
reports ratios 1.420/1.233. The small-input cost remains unresolved. These are shared-host
measurements, and v4 remains unaccepted and outside the public candidate.
Evidence: `/tmp/hybrid-isothermal-general-reduction-v4/experiment-evidence.json`.

The v4 small-input batched repeat passes 70 assertions on latest Julia, with
100 calls per sample and 11 alternating pairs. Its ratios are 0.870/0.918; these
do not reproduce the slower single-call samples. Both sets remain recorded in
`/tmp/hybrid-isothermal-general-reduction-v4`, including
`small-varied-batch-evidence.json`.

Near-unit v5's five-pair latest timing run passes 44 assertions with zero measured
allocations. Its large Float64 ratio is 0.620; the large Float32 ratio is 1.269,
and small-input ratios are 1.011/0.244. The Float32 arithmetic is unchanged;
these shared-host timings establish no universal speedup. General v5's separate
latest run passes 132 paired output/resource assertions. The public-v3 paired
comparison is recorded below; complete public acceptance remains pending.

## Historical pre-promotion gates

Complete both public-v2 native lanes and retain the measured timing costs
against the exact previous public owner as an open optimization task. Native lanes include solver/restart, quality/Aqua,
threaded and full compiler-tree comparisons. Verify the real extreme-density
resource controls and preserve failed evidence. The correction is pushed as
`6e7fb01`; its original exact-SHA run later failed, while the separate traversal
run passed as recorded above. No current acceptance is inferred from this old
pre-promotion checklist.

The [polytropic plan](POLYTROPIC_ENERGY_RANGE_PLAN.md) and the
[complete budget/inertia scope](HYBRID_DIAGNOSTIC_RANGE_PLAN.md) remain separate
open work. Correcting isothermal energy does not close them.

Public v3 paired timing on Julia 1.13.2 passes 240 assertions over 16 cases at
4,096 and 4,194,304 cells. Both implementations return identical values against
independent bounds; all 48 measured calls per implementation allocate zero bytes.
Uniform-density v3/v2 median ratios are 0.131–0.209. Float64 weak-wave ratios are
0.468 and 0.664; Float32 weak-wave ratios are 1.023 and 1.124. The latter measured
cost is not reproduced consistently: an 11-pair repeat passes 92 assertions,
with ratios 1.005 (small) and 0.928 (large), identical values and zero bytes in
all measured calls (`repeat-evidence.json`). These shared-host measurements do
not establish a universal speedup. Exact source/commands and per-case data are
in `/tmp/hybrid-isothermal-public-v3-timing/evidence.json`.

Local broad acceptance now also runs the complete selected suites in separate
checkpoint and solver/diagnostic processes. The original combined runs remain
untouched and pending. All 13 files and their assertions are retained; production
and test hashes are unchanged. V2 split evidence is recorded under
`/tmp/hybrid-isothermal-public-v2-split-verification`, with the same scheme queued
for v3. Focused, concrete, quality, threaded and unfiltered compiler-tree gates
remain required, followed by exact-SHA full-package GitHub CI. A separate typed
checkpoint-test helper passed its original 276 assertions but showed no timing
improvement; it is not included in any repository or acceptance snapshot.

The checked 15-file isothermal-v2 correction is pushed as `6e7fb01` for
[exact-revision GitHub CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129297310)
while both native broad lanes continue. This overlaps the complete
Linux package/format/MPI run with native solver, restart, quality, threaded and
full compiler-tree checks; it does not accept or remove any pending gate.
The 15,427 permanent assertions plus one smoke check and 72 concrete compiler
checks already pass on each runtime. The performance-v3 and pressure candidates
remain separate. The reference to `d497436` in this campaign identifies its then-green baseline.

The separate public-v3 traversal change is pushed as `5323c4d`;
[its exact-revision CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39129564503)
finished with all six jobs successful.
Its four Julia paths add the checked Float64 degree-eight first pass and local
repeated-density reuse, with the original complete fallback and accumulation
order retained. Both native versions pass 27,831 permanent assertions plus one
smoke check and 144 concrete compiler checks. The 240-assertion paired run and
92-assertion Float32 repeat above pass. Broader native acceptance remains queued
and mandatory in that historical handoff; the completed GitHub outcome is
recorded above.
