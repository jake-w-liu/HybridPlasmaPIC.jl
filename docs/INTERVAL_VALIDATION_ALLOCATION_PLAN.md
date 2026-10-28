# Interval validation allocation

Revision-scoped implementation record. The baseline failures and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

The repair landed as `7be1916`; [all six exact-SHA jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39895197004)
(checked 2026-10-29). The local acceptance results below retain their actual
Julia 1.12.7 / 1.13.2 versions.
Parent: `e8925e8663883832b6ffe453a559ac7be1903ae3`.

The retained `migration-seam/boundary-allocation-probe.jl` commands measure
640 bytes per axis on Julia 1.12.7 and 384 on 1.13.2, for both Float32 and
Float64. Empty `apply_periodic!` and `_validated_open_interval` allocate the same
amount in all three samples. The ordinary comparison confirms this cost for
nonempty populations too. This is an allocation measurement, not a latency
or whole-simulation bottleneck claim.

Pre-Edit Brief:

- Mode: bounded allocation improvement, preserving numerical behavior.
- Requirement: remove successful Float32/64 interval-validation allocation;
  retain public output bits, RNG progression, exception types/messages/order,
  conversions and rejection before particle mutation.
- Source: `particles/loading.jl:8` validates all lower bounds, then all upper
  bounds, then ordering. It constructs indexed labels on every finite input.
  `utils/validation.jl:8` owns conversion and finite rejection. The interval
  owner serves uniform/lattice loading and periodic/reflecting/absorbing
  boundaries. Velocity validation is outside this change.
- Mechanism to test: passing the label and index separately to the existing
  finite validator can defer label construction until rejection. The current
  allocation probes do not by themselves attribute every byte to strings.
- Canonical change: extend the existing private finite validator with an
  optional index; its existing three-argument callers retain the same owner.
  Change only the interval owner's two indexed calls. No new public API,
  cache, parallel validation implementation or dropped check.
- Compatibility: all valid operations retain their arithmetic and conversion
  order. Invalid lower/upper bounds retain indexed text and precedence. The
  existing non-indexed helper remains supported. BigFloat allocation is not
  required to be zero; its numerical semantics remain required.
- Verification: native minimum/latest public and error matrix, late-bound
  failures with full particle/RNG rollback, warmed allocation measurements,
  neighboring loaders/boundaries/deposition, concrete and complete package
  compiler comparisons, formatting and exact-source checks. Compare complete
  baseline/candidate outputs and error evidence before separate main promotion.
- Incomplete fixes to reject: skipping finite checks (NaN/Inf tests); validating
  after mutation (late-axis state/RNG tests); dropping axis labels or changing
  validation order (exact-message and multiple-invalid-bound tests).

- [x] Preserve measured parent allocation evidence.
- [x] Implement in an isolated worktree and preserve original source snapshots.
- [x] Pass numerical/error/resource and compiler gates on both runtimes.
- [x] Push separately as `7be1916`; all six exact-SHA GitHub CI jobs pass.

Final candidate: `../HybridPlasmaPIC-interval-validation-20261023`. Both native
neighboring suites pass 21,723 assertions. The same 2,405 indexed-message,
validation-order and whole-state/RNG rejection checks also pass on the unchanged
parent. All 40 concrete compiler calls have no reports; complete package report
groups and order match the parent (30/27 existing reports), with none suppressed.

The four ordinary comparison commands pass 1,350 assertions each. All 270
complete position/velocity/weight/id/tag/species/result/RNG samples match per
runtime across the five public operations. Every candidate sample allocates
zero bytes, compared with 640/1,280/1,920 parent bytes on minimum Julia and
384/768/1,152 on latest. Arithmetic, mutation and RNG calls are unchanged; only
the finite validator's indexed label construction moves to its error branch.
No elapsed-time claim is made.

Eight final formatter checks pass. The first formatter AST comparison retained
one failure because `Base.remove_linenums!` leaves macro line-number arguments;
both complete AST displays were checked and differ only in those source-line
annotations. The failure and corrected final idempotence result remain retained.
`interval-validation/final-verification.json` checks all accepted native logs,
both frozen source manifests, the four owned files and every unchanged
committed Julia file. No held radiation changes are included.
