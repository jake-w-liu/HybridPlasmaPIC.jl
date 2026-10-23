# Complete BigFloat rational conversion

Status: **NUM-002 V1 pushed as `245e73a`; all six exact-SHA CI jobs and complete compiler review pass.**
This shared arithmetic repair supports existing exact fallbacks. It does not
close the separate normalization, linear-CN, concurrency or capability work.

## Contract and checked failure

The canonical owner is `_from_exact(BigFloat, ::Rational{BigInt})` in
[exact_arithmetic.jl](../src/utils/exact_arithmetic.jl). Its callers in AMR,
gyrokinetics, moment deposition, field/electron diagnostics and shock diagnostics
require rounding the complete finite rational at the caller's BigFloat precision,
rounding mode and exponent limits.

**VERIFIED:** both Julia 1.12.7 and 1.13.2 reproduce 212 passes / 111 failures
across 160 baseline profiles. With MPFR exponent limits [-64,64], the old owner
converts numerator and denominator separately; `(2^200+1)/(2^200-1)` becomes
NaN although its complete rounded quotient is finite and near one. The witness
computes its reference from exact integers before restricting the exponent range.
All commands, source hashes, outputs and failed logs remain in
`HybridPlasmaPIC-audit-evidence/2026-10-23/exact-bigfloat-conversion/`.

## Canonical implementation

Retain the exact-integer division path when both integer operands fit the active
MPFR exponent range. Otherwise quantize the complete rational to a bounded integer
significand with guard and sticky bits; pass that integer and its binary exponent
together to `mpfr_set_z_2exp` for the final rounding. Production code must never
change global precision, rounding or exponent settings. Zero remains exact even
when one is outside a compact positive-minimum range.

The finite-rational caller contract remains unchanged. No public conversion API,
CN-specific materializer or normalization-only workaround is added. Integer work
is bounded by the supplied rational's bit lengths and requested precision; this
helper does not justify expanding arbitrary BigFloat exponents into huge rationals.

The MPFR manual cautions that [`mpfr_set_q`](https://www.mpfr.org/mpfr-current/mpfr.html#Assignment-Functions)
may fail when integer factors cannot individually be represented, so using it
alone is not accepted as a solution to this defect.

## Acceptance tracker

- [x] Inspect the canonical helper, all call sites and adjacent exact arithmetic.
- [x] Reproduce the unchanged baseline on both supported Julia versions.
- [x] Prepare V1 in an isolated worktree; preserve held radiation changes.
- [x] Pass all 323 original checks on both runtimes without changing the witness.
- [x] Add and pass 15,150 independent endpoint/public-owner checks per runtime.
  Exact rounding-cell inequalities cover precisions 2/3/24/53/128/256, all five
  public rounding modes, signs, tie sides, overflow/underflow and four compact
  exponent ranges, including ranges that exclude one. Public `cell_center` and
  `gyroaverage` exercise retained-precision inputs through the shared owner.
- [x] Pass the existing 49,080 exact-arithmetic checks per runtime; the child
  fixture uses the resolved active project and checks the loaded source path.
- [x] Reject lost-sticky, lost-guard and forced-nearest wrong-patch controls:
  they violate 176, 800 and 688 independent rounding-cell checks, respectively,
  on each runtime.
- [x] Pass 4,201 concurrent checks per runtime, with each static worker's TLS
  exponent limits explicitly configured and restored. The original fixture's
  assumption that main-thread bounds propagate to workers was disproved by
  native before/after probes; its 1,804 passes / 2,367 failures are retained as
  invalid fixture evidence. No production change or assertion relaxation followed.
- [x] Pass 1,081 resource assertions over 270 exceptional profiles per runtime,
  including 65,536-bit integer factors and precision up to 4,096 bits. Maximum
  allocation is 41,864 bytes; all first exceptional-input observations have
  zero compilation after ordinary-signature warmup.
- [x] Pass 148,054 neighboring assertions, five concrete compiler checks, nine
  Aqua checks and three JET checks per runtime. Both complete package trees
  match the parent (29 / 27 reports) after checked source-line shifts only.
- [x] Execute the frozen ordinary 75-profile cost matrix once per runtime.
  All complete outputs and controls pass, with zero extra warmed allocation.
  Maximum ratio of medians is 1.170083 / 1.183004; maximum paired-median ratio
  is 1.125468 / 1.087506, within the unchanged 1.25 limit. These are measured
  repair costs, not a speedup claim. No timing reruns or sample filtering.
- [x] Review the exact final source through the caller contracts and a separate
  guard/sticky rounding argument; find no additional confirmed defect in this
  bounded conversion pass. No full-simulator or general MPFR guarantee follows.
- [x] Push this repair separately as `245e73a`; all 366 committed Julia/Project
  files match the accepted source snapshot, and held radiation hashes are unchanged.
- [x] All six [exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40026749996) pass.
  Downloaded complete logs verify 625,494 / 625,495 assertions on Julia 1.12.7 /
  1.13.2. Full compiler trees match reviewed `ba3694c` (24 / 21 reports) after
  checked unchanged-source line mappings. Archive, log and tree hashes are retained.
  This closes this repair's CI gate, not the full simulator audit.
