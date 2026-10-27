# Shock boundary energy-flux range

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: ENERGY-008 is repaired in the isolated canonical V2 integration.
Native focused, concurrency, resource, shock-regression and full compiler
comparison gates pass on both runtimes. Exact-SHA CI acceptance remains open.

## Evidence and owning path

**VERIFIED:** both runtimes pass 34 witness assertions against the shock-flux,
shock-constructor and particle-storage sources byte-identical to `15662e1`.
`boundary_energy_flux` returns infinity or zero for finite complete bulk and
particle kinetic fluxes. Its electron enthalpy flux also returns infinity for
a finite complete product. Evidence, exact rational fixtures, source hashes,
commands and terminal results are in
`/tmp/hybrid-shock-energy-flux-range-20261018` (`source-proof.json`,
`witness.jl`, `results.json`, `witness112.log`, `witness1132.log`).

For example, Float64 density/weight `2^-1000` and x velocity `2^600` give
bulk flux `2^799`, while the current owner returns `Inf`. Velocity `2^-600`
and density/weight `2^1000` give `2^-801`, while the owner returns zero.
The particle witness includes the stored SBP boundary norm before rounding.
With gamma 1.5, pressure `floatmax(T)/2` and x velocity 0.25, complete enthalpy
is finite while the intermediate pressure product overflows.

The canonical owners are `boundary_energy_flux` and `_ion_energy_flux` in
`src/diagnostics/shock_surface.jl` (lines 42–103 at the checked revision).
The existing `test/test_shock_diag.jl` already checks constructed boundary
states, H-weighted kinetic-flux conservation and invalid particle velocities.
This is separate from volume-integrated `electron_internal_energy`.

A further 18 nodal witness assertions pass per runtime (`nodal-results.json`).
They establish lost subnormal contributions before accumulation, overflowing
opposite contributions whose complete wall flux is one, NaN in a zero-weight
neighbor, and finite flux lost before division by a large SBP norm. These cases
require retaining the complete nodal expression, not only repairing a scalar
per-particle value.

## Repair contract and plan

- Retain each finite stored factor through bulk kinetic and electron enthalpy
  flux evaluation when native intermediate arithmetic loses range. Reuse the
  shared exact conversion owner; preserve ordinary arithmetic where certified.
- For particle flux, include every deposited term, CIC weight and final SBP norm
  in the complete nodal expression. Correcting an individual particle's flux
  after materialization is insufficient when deposition, accumulation or
  normalization changes the representable range.
- Preserve generic/nonfinite behavior and the existing singular gamma-one
  enthalpy convention. Keep particle validation and exception order, inputs,
  deposited node locations, boundary ordering and the current mass convention.
  This repair changes arithmetic, not the shock model or its evolution.
- Audit scalar total-flux combination separately before extending the contract;
  do not silently claim a complete total when only its components are repaired.
- Keep ordinary allocation equal to the existing output-vector requirement.
  Measure complete cold-call cost, temporary storage and grid/particle scaling.

Acceptance must cover the original witnesses, signed and zero factors,
subnormal and overflow boundaries, exact/neighboring ties, zero CIC weights,
multiple particles whose individually rounded contributions lose the nodal
sum, views, empty populations, validation failures and unchanged inputs.
Reject patches that only check `isfinite`, only repair squared velocities,
or round per-particle flux before its stencil and H normalization.

Run independent rational nodal/flux oracles, all shock diagnostic and solver
regressions, ordinary bitwise/allocation comparisons, concrete and complete
compiler checks, Aqua, four-thread shared-input checks and formatting on both
supported runtimes. Commit and push the verified repair separately; require
green CI for that exact SHA. This item remains open until those gates pass.

## Implementation brief

- Mode: canonical arithmetic repair; prototype first outside the worktree.
- Acceptance: restore all 52 witnessed bulk/enthalpy/nodal cases per runtime;
  preserve ordinary bitwise results and output-only allocation; retain full
  validation, generic and nonfinite behavior; complete numerical, compiler,
  concurrency, solver and exact-SHA CI gates.
- Verified cause: the checked source materializes squared/cubic velocities,
  pressure products and per-particle flux before reducing factors, nodal
  cancellation or the stored SBP norm can act. The rational witnesses and
  byte-identity evidence above isolate those operations.
- Canonical owners: `boundary_energy_flux` and `_ion_energy_flux`. Use existing
  `_from_exact` and normal-product/error-bound arithmetic. No exported API,
  shock equations, mass convention, stencil location or validation changes.
- Strategy under test: retain ordinary scalar operation order behind
  conservative factor-range certificates. For deposited kinetic flux, track
  an outward global absolute-error bound while retaining the original sum;
  certify each normalized node or recompute its complete rational sum. Empty
  uncertain nodes must become zero without allocating an exact accumulator.
  Exact work is local to the call and only allocated for touched uncertain nodes.
- Compatibility: nonfinite/generic fields keep current arithmetic; gamma-one
  enthalpy remains singular. Particle input rejection order remains unchanged.
- Verification: independent nodal rational oracles, original and neighboring
  failures, existing shock suites, baseline bitwise/allocation comparison,
  measured cold resources, concrete/full compiler reports, Aqua, four-thread
  shared-input checks and formatting on 1.12.7/1.13.2. Review the full exact diff.
- Guardrail: an uncertain numerical certificate requests complete arithmetic;
  it never suppresses an input error or accepts a degraded answer. Certificate
  derivation and actual resource behavior must be verified before promotion.

## Prototype evidence and retained failures

The production shock owners remain unchanged at `7104c61`. Prototypes are
isolated under `/tmp/hybrid-shock-energy-flux-v1`, `-v2` and `-v3`.

- V1's first oracle inherited the directed MPFR rounding mode during destination
  conversion: 170 checks passed and six failed. Explicit nearest-even conversion
  corrects the oracle; the original logs and Julia MPFR source classification
  are retained. The unchanged algorithm then passes 316 numerical checks but
  fails one ordinary allocation check on 1.12.7 (6,896 B versus 320 B).
- V2 adds a checked dyadic-quantum certificate and passes 748 initial checks on
  both runtimes. The unchanged larger gate fails six assertions per runtime:
  three Float32 cases change ordinary rounding and allocate 6.2–45.1 MB.
  These failures remain in `large-resources-results.json` and their full logs.
- V3 accounts for Float32 product and addition errors with outward Float64
  intervals. Both runtimes pass the same 748 initial checks and all 24 larger
  checks. All 12 larger cases per runtime retain bitwise outputs and exactly the
  baseline output-vector allocation (320–4,160 B); no resource assertion was
  removed or relaxed. Single-call Float32 latency samples are slower than the
  baseline, so this is not a speedup claim or completed optimization.
- Independent rational checks pass 12,860 assertions per runtime: 5,006 outward
  interval/error checks, 2,400 exact minimal-quantum checks and 5,454 complete
  nodal, folded-stencil, cancellation, subnormal, midpoint and overflow checks.
  The nodal oracle constructs its stencil independently for dyadic fixtures;
  final values are bracketed with directed 8,192-bit MPFR before explicit
  nearest-even conversion. Input immutability is checked.

V3 commands, terminal exit codes, elapsed times and source/harness hashes are
retained in `results.json`, `large-resources-results.json`,
`properties-results.json` and their matching frozen-hash manifests. These are
component gates; the remaining acceptance steps above still apply.

## Canonical integration and measured cost

V4 replaces the Float32 interval calculation with a seven-operation normal
product error bound plus half-ULP addition bounds, accumulated outward in
Float64. Its independent rational checks pass 7,981 assertions on each runtime.
The original initial, larger and full-expression gates also pass. V3's now
unused interval helpers are excluded from production; their evidence remains
frozen. V4's complete compatibility/compiler/concurrency suite passes 536
assertions per runtime. The first contract harness had a metadata-tuple parse
error; its source and both failed logs are retained alongside the corrected
harness, with no algorithm or assertion changes.

The first canonical migration missed three shorthand `2factor_quantum` calls.
The independent AST comparison fails 1 of 20 checks, concrete JET fails eight
calls per runtime, and the full 1.12 report comparison fails. Those outcomes
remain in `/tmp/hybrid-shock-energy-flux-public-v1-verification`. Canonical V2
uses explicit multiplication with the renamed helper and adds a runtime test
of the affected lattice-certificate path; no diagnostics are suppressed.

Canonical V2 evidence is in
`/tmp/hybrid-shock-energy-flux-public-v2-verification`:

- Both native versions pass 17,169 permanent assertions across ten testsets,
  plus an explicit four-thread acceptance assertion. Coverage includes public
  boundary calls, complete nodal arithmetic, independent bounds, ordinary
  output/allocation equality, strided storage, empty populations, Float16 and
  BigFloat behavior, validation order, concrete JET and shared-input ownership.
- Both pass all nine existing shock files (453 assertions) and Aqua (nine).
- Both complete compiler trees match baseline: 22 reports on 1.12.7 and 21 on
  1.13.2, with no added or removed diagnostic groups.
- Twenty AST comparisons establish prototype migration, unchanged neighboring
  functions, the two intended public call-site changes and the exact legacy
  test comparator. Nine formatting comparisons preserve executable syntax;
  the whole-repository read-only formatting scan passes.
- Balanced ordinary measurements alternate three implementations over twelve
  repetitions; complete cold calls cover 16, 256 and 4,096 particles. All 1,044
  numerical/resource assertions pass per runtime. Ordinary allocation remains
  exactly the legacy output-vector allocation.

Measured median-call ratios against the legacy range-limited owner are
3.16–5.59 for Float32 and 1.84–2.29 for Float64 on 1.12.7; on 1.13.2 they are
3.26–5.88 and 2.18–2.78. Float32 ratios against V3 are 0.39–0.51 and 0.36–0.69.
The largest measured cold cumulative allocation is 27,403,472 B. Maximum
median cold time among these cases is 45.2 ms / 77.7 ms. These local samples
include normal runtime variability and establish neither peak memory nor
portable latency. Further ordinary-path and exact-accumulator optimization
remains open, with the numerical and output-only allocation gates retained.

## Complete total — ENERGY-009

The separate total-flux audit confirms that ordinary addition of already
rounded components loses representable complete totals. Both runtimes pass
36 witness assertions in `/tmp/hybrid-shock-total-flux-range-20261019` against
the frozen canonical V2 source. Finite positive density/pressure and finite
fields can produce opposite overflowing magnetic/kinetic components whose
complete total is signed one, while the reported total is NaN. Three exact
half-minimum-subnormal components round separately to zero while their joint
total rounds to two minimum subnormals. Both bulk and deposited-ion total
paths are covered by the cancellation witness.

This is a follow-up repair, not covered by ENERGY-008. Plan: retain the
complete magnetic, kinetic and enthalpy expression through the final total
rounding when range is uncertain; preserve the existing component tuple,
ordinary total order, singular gamma-one convention and validation behavior.
The particle path must retain boundary CIC contributions and H normalization
before combining with fields. Reuse canonical component owners, avoiding a
second public implementation. Require exact independent total oracles,
neighboring ties/cancellation, public callers, ordinary bitwise/allocation
controls, measured cold scaling, concurrency, full compiler comparisons,
shock regressions and green exact-SHA CI before marking this item complete.

The total-flux pre-edit brief uses the same canonical boundary owner and
retains the public return shape. The verified cause is component
materialization before the addition at `shock_surface.jl:156` in `daea194`.
The isolated experiment `/tmp/hybrid-shock-total-flux-v1` extracts the three
existing exact component expressions into shared private helpers. The ion
owner returns its flux and raw error bound internally; its existing vector
interface remains a wrapper over that owner. A normal-range certificate
combines field and ion bounds, while an uncertain total replays complete
boundary contributions and final normalization. This design must preserve
ordinary bitwise outputs and allocation, including empty and untouched
boundaries. Initial exact-total and allocation gates are running; the design
is not yet accepted or copied into production. No model equations, mass
convention, stencil, validation order, singular gamma-one or generic/nonfinite
behavior changes are authorized by this arithmetic repair.

Total prototype V1 passes 96 exact-total checks per runtime but fails 60
ordinary allocation assertions (84 pass). The retained 1.12 allocation
profile records a `Core.Box` and boxed Float64 values in the total-result
closure. Its optional lowered-code lookup subsequently fails because the
script guessed an absent generated method name; that failure is retained.
V2 replaces the closure with explicit calls for the two API boundaries.
With every original assertion unchanged, both runtimes pass all 240 initial
checks, including zero-byte bulk calls and output-vector-only particle calls.

The V2 independent numerical suite passes 6,480 assertions per runtime:
1,740 complete cancellation/midpoint/overflow checks, 2,000 complete-field
error-bound checks and 2,740 ion-state bound/immutability checks. Its adapted
compatibility/concrete-compiler/concurrency suite passes 546 assertions per
runtime, plus the explicit four-thread gate. Commands, hashes and original
failures remain under `/tmp/hybrid-shock-total-flux-v1` and `-v2`. Canonical
integration, full compiler comparison, cost measurements, regressions and
exact-SHA CI are still required for ENERGY-009.

The canonical total-flux integration is now verified in
`/tmp/hybrid-shock-total-flux-public-v1-verification`:

- Both native versions pass 24,737 focused assertions across 17 testsets,
  including the unchanged 17,169 component assertions and 7,568 total checks.
  The native driver explicitly requires four threads. The total concurrency
  suite contributes 812 parent-checked assertions per runtime.
- All nine existing shock files and Aqua pass 462 assertions per runtime.
  Complete compiler trees match baseline (22/21 reports), with no added or
  removed diagnostic groups. The final whole-repository formatting scan passes.
- Thirty-three AST checks preserve prototype arithmetic and all neighboring
  method overloads. Nine further checks establish that the measured prior
  owner is the pushed component implementation. Eight formatting checks
  preserve executable syntax. The first expanded verifier did not parse an
  existing parametric constructor; its corrected parser retains every overload,
  and the failed harness and log remain available.
- Both cost suites pass 414 assertions. Twelve ordinary cases retain bitwise
  results and allocation counts. Alternating paired median-call ratios against
  the pushed component owner are 1.25–1.72 for bulk and 0.99–1.56 for particle
  calls on 1.12.7; on 1.13.2 they are 1.13–2.00 and 0.87–1.75. Eighteen complete
  cold cases per runtime cover 16–4,096 particles, overflowing cancellation,
  subnormals and normalization. Their maximum cumulative allocation is
  47,985,808 B; the largest measured median is 94.0 ms / 219.1 ms. These are
  local samples, not peak-memory or portable-speed claims.

The internal ion result carries a raw reduction-error bound plus the flux
vector. A zero bound after global exact recomputation refers to reduction
error; the total certificate still includes final conversion/normalization
rounding. Generic or invalid-norm paths retain ordinary behavior. Exact bulk,
enthalpy and particle expressions have one shared private owner each. The
existing ion-vector interface delegates to the same accumulator. No public
API, component-output contract, validation order or physics equations change.
Further cold-storage and ordinary-latency optimization remains open. Exact-SHA
GitHub acceptance is required after this repair is pushed.

The component repair is pushed as `daea194035e85cd39c5a497e781f758844881ec4`. [Run 39144368129](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39144368129) passes all six jobs: 491,943/491,944 assertions on Julia 1.12.7/1.13.2, formatting and all MPI ranks. Full compiler trees remain unchanged at 23/22. The exact 14-path promotion and remote SHA were verified in `committed-evidence.json`; held radiation files retain their hashes. Exact-SHA logs and full-report comparisons are in `/tmp/hybrid-ci-daea194-20261019` and `/tmp/hybrid-ci-daea194-report-review-20261019`. No billing exception is used.

The complete-total repair is pushed as `dd539da0f267f5d18481d0dae7a68111adc0b5d8`. [Run 39146029176](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39146029176) passes all six jobs, including 499,511/499,512 assertions on actual Julia 1.12.7/1.13.2. Complete compiler trees match `15662e1` at 23/22 reports; exact-SHA logs and comparisons are retained in `/tmp/hybrid-ci-dd539da-20261019` and `/tmp/hybrid-ci-dd539da-report-review-20261019`. All 12 committed blobs and remote main match `committed-evidence.json` in the canonical total verification directory. Held radiation files remain unchanged; this completes ENERGY-009 acceptance without a billing exception.

## Exact-node reuse optimization brief

- Mode: measured resource optimization, with all numerical contracts retained.
- Baseline: `/tmp/hybrid-shock-exact-reuse-review-20261019` passes 48 assertions
  per runtime against the pushed total owner. At 4,096 particles, the two
  isolated boundary replays together allocate 19,432,880 B (Float32) or
  22,906,544 B (Float64). The component state has already accumulated those
  complete nodes in its local dictionary before the total traverses them again.
- Owner and contract: retain the same canonical ion accumulator and exact
  component helpers. Expose already-computed exact nodes only to the enclosing
  boundary calculation, with call-local ownership and no persistent cache.
  A missing dictionary entry must still trigger complete boundary evaluation;
  it never means zero. Preserve all public signatures and output semantics.
- Strategy: prototype a retained optional exact-node result. Measure ordinary
  return/storage overhead before choosing its internal representation. Reuse
  present entries before H normalization; retain the canonical evaluation for
  uncached nodes. No new workspaces on the certified ordinary path.
- Acceptance: unchanged numerical and bitwise outputs, ordinary allocation,
  validation and generic behavior; explicit partial-dictionary and zero-entry
  controls; four-thread ownership; complete compiler comparisons, all shock
  regressions and formatting. Target cumulative cold allocation at most 60%
  of the checked baseline for the 256-/4,096-particle cancellation cases, and
  record balanced latency rather than infer a speedup. Push only after these
  gates, then require exact-SHA green CI.

The retained V1 prototype passes its numerical checks but adds 32 B to ordinary
particle calls: 90 ordinary/resource assertions pass and 54 fail per runtime.
**VERIFIED:** the allocation trace identifies the three-element returned tuple at the state
owner; its log is `/tmp/hybrid-shock-exact-reuse-v1/allocation-trace112.log`.
V2 replaces that tuple with a concrete immutable internal result containing the
flux, raw error and optional exact-node dictionary. All 240 initial assertions
now pass on each runtime, including the original ordinary allocation checks.

V2 also passes 6,480 numerical/bound assertions, 812 public-total concurrency
assertions and 234 partial-cache/ownership/concrete-compiler assertions per
runtime. Both cost suites pass 480 assertions. The 256-/4,096-particle cold
cases allocate 50.4–50.8% of the prior owner on both runtimes; all ordinary
bitwise and allocation comparisons pass. Eighteen measured cold profiles per
runtime include normalization and subnormals. These are local cumulative
allocation measurements, not peak-memory or portable-latency claims. Logs,
commands and hashes are under `/tmp/hybrid-shock-exact-reuse-v2`.

The canonical integration is verified in
`/tmp/hybrid-shock-exact-reuse-public-v1-verification`. Only the private ion-state
and total-helper callers migrate; public signatures and the ion-vector wrapper
remain intact. Permanent regressions reject missing-key-as-zero, reusing a
normalized value, shared cache ownership and a repeat traversal of exact
boundary contributions. Existing numerical, generic/error-path and ordinary
resource assertions remain unchanged. The final local gates are recorded below.

Canonical local acceptance now passes on Julia 1.12.7 and 1.13.2:

- 25,091 permanent assertions across 19 testsets, retaining all 17 earlier
  component/total testsets and adding partial-cache/ownership and reuse-resource
  regressions. The native driver explicitly checks four threads.
- 462 assertions across the nine existing shock files and Aqua; 480 measured
  cost assertions against the AST-verified pushed `dd539da` owner.
- Complete compiler trees match baseline at 22/21 reports, with no added or
  removed groups. All 58 migration/prior-owner and six format-preservation
  checks pass; the whole-repository read-only formatting scan passes.
- All 12 ordinary profiles retain bitwise outputs and allocation. Paired median
  latency ratios are 1.00–1.20 (bulk) and 0.91–1.07 (particles) on 1.12.7;
  0.91–1.14 and 0.94–1.03 on 1.13.2. The 256-/4,096-particle cold allocation
  ratio is 0.5042–0.5085 on both runtimes, below the declared 0.60 target.
  A separate alternating-order cold measurement passes 246 assertions per
  runtime. Its median latency ratios span 0.22–1.87 / 0.27–1.49 across the
  sampled profiles. Shared-machine timings are mixed and do not establish a
  portable speedup. The largest measured new cumulative allocation is
  24,196,048 B. Earlier block-order measurements remain in the original cost logs.

The six tested source/test paths are copied to main with matching hashes;
all held radiation files remain unchanged. Public APIs, ordinary arithmetic,
complete-total rounding and validation are retained. There is one canonical
accumulator, one call-local optional dictionary and no compatibility overload
for the private helper signature. Uncached nodes still use the existing complete
evaluator. Commands, assertion inventories, prior-owner AST proof and logs are
in `/tmp/hybrid-shock-exact-reuse-public-v1-verification`. This local acceptance
still requires a separate exact-SHA GitHub run after push.

Exact-node reuse is pushed as `99760f0193a40c3a4bdb18ba9d1013e11a49d0ea`. [Run 39147850360](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39147850360) passes all six jobs: 499,865/499,866 assertions on actual Julia 1.12.7/1.13.2, formatting and all MPI ranks. Complete compiler trees match `15662e1` at 23/22 reports. All 11 committed blobs and remote main match `committed-evidence.json` in the reuse verification directory; held radiation hashes remain unchanged. Logs and exact-SHA evidence are in `/tmp/hybrid-ci-99760f0-20261019`, with full comparisons in `/tmp/hybrid-ci-99760f0-report-review-20261019/summary.json`. No billing exception is used.
