# Threaded deposition storage checks

Revision-scoped implementation record. The baseline failures and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

## Confirmed baseline failure

**VERIFIED:** `deposit_scalar_threaded!(vals, ps, vals, g, NGP())` clears
`[1,2,3,4]` to zeros instead of rejecting the overlapping output and values.
The serial API raises `ArgumentError` and preserves the values. The independent
probe reproduces this with one thread on Julia 1.12.7 and four threads on Julia
1.13.1: `/tmp/hybrid-threaded-alias-baseline1.log` and
`/tmp/hybrid-threaded-alias-baseline4.log`, each with four passing witness
assertions. These are defect witnesses, not repair acceptance tests.

The threaded header checks sizes and positions but omits the serial alias guard
before `fill!`. Its public contract promises serial-equivalent deposition;
`density_threaded!` reaches the same unguarded path.

## Pre-edit brief

- **Mode:** repair threaded input validation.
- **Acceptance:** reject overlapping output/value or particle storage before
  any mutation, with the same error types as serial deposition.
- **Baseline:** the executable witnesses and missing header guard above.
- **Cause:** the threaded implementation duplicated only part of the serial
  preflight and then zeroed the output before reading aliased values.
- **Contract:** preserve accepted numerical inputs, stencil construction,
  partitioning, reduction order, offset-indexed values and empty-particle
  behavior. Validation failures preserve all caller-owned arrays.
- **Owner:** shared `_validate_deposition_inputs` in `src/coupling/deposit.jl`,
  called by serial and threaded scalar deposition. Both density APIs inherit
  the same check; existing particle-storage alias helpers remain authoritative.
- **Compatibility:** reject unsafe threaded aliases with `ArgumentError`; no
  public signature or valid-result change.
- **Strategy:** extract the existing serial preflight, invoke it before writes
  in both implementations, and document the threaded storage restriction.
- **Verification:** Float32/64, dimensions 1–3, all built-in shapes, direct and
  view aliases, particle weight/position/velocity aliases, empty value views,
  unchanged storage, serial/multithread parity, concrete/package compiler gates,
  allocation checks, formatting and GitHub CI.

## Checks against incomplete repairs

Checking only object identity misses views and particle aliases. Checking only
inside the multithread branch misses the one-thread path. Checking after output
initialization can throw while corrupting caller storage. Returning early for
empty inputs can bypass backing-storage checks. The new tests exercise each of
these cases and retain all existing numerical and error tests.

## Progress

- [x] Reproduce the discrepancy with one and four threads.
- [x] Share the existing validation owner and add rollback coverage.
- [x] Observe regression tests on both Julia versions with one and four threads.
- [x] Observe baseline failures of the new tests.
- [x] Check compiler analysis, allocations, formatting and exact diff.
- [x] Commit the verified repair.
- [x] Push the repair (`a2154cf`).
- [x] Observe the repair included in fully green `240f77c` acceptance.

This repair does not address the separate moment-normalization range defects;
those remain tracked in [the charge and moment plan](CHARGE_DEPOSITION_CONSERVATION_PLAN.md).

The four-thread suite passes **754 assertions** on each native Julia version.
The new storage testset contributes **648**; against the original source it
reports **324 passes and 324 failures**, including failure to reject aliases and
mutation of caller arrays. Logs: `/tmp/hybrid-threaded-alias-repair4-{112,113}.log`
and `/tmp/hybrid-threaded-alias-red113.log`. Final processes exit 0 / 0 / 1,
respectively. The later one-thread and compiler/allocation results are recorded below;
`/tmp/hybrid-threaded-alias-hashes.json` records the frozen source.

The warmed allocation comparison passes **36 assertions** for 18 configurations
on each of Julia 1.13.1 with one thread and Julia 1.12.7 with four threads.
All ordinary outputs match the baseline exactly, and every allocation count is
unchanged. The single-thread calls allocate **0 B**. These checks do not claim a
latency improvement. Logs: `/tmp/hybrid-threaded-alias-alloc-final1-113.log` and
`/tmp/hybrid-threaded-alias-alloc-final4-112.log`. The first temporary benchmark
module omitted `nparticles`/`width` imports; its driver failed before measurement.
The corrected driver uses the unchanged baseline file and exits 0 on both runs.

The three edited Julia files are formatter-idempotent at margin 100
(`/tmp/hybrid-threaded-alias-format-final.log`). The exact three-file production/test diff matches the frozen snapshots;
whitespace checks also pass.

Each native version passes **3,012 one-thread integration assertions**, covering
serial deposition, moment storage and the threaded suite. Logs:
`/tmp/hybrid-threaded-alias-repair1-{112,113}.log`, both exit 0. Together with the
four-thread runs, these exercise both threaded execution branches.

Each version also passes **three package JET assertions**, **nine Aqua checks**
and **24 concrete compiler analyses** with no reports. Package report counts
remain **22 / 21**, under the unchanged ceiling of 30. Evidence:
`/tmp/hybrid-threaded-alias-gates{112,113}.log`, both exit 0. The concrete gates
cover both deposit APIs with vectors/views and both density APIs, Float32/64 and
1D/3D. No numerical loop, thread partition or reduction was changed.

The `_validate_deposition_inputs` function and threaded regression file are
byte-identical between repair `a2154cf` and green `240f77c`; the surrounding
deposition kernel changed separately for MOMENT-001. [GitHub run 38741178196](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38741178196)
passes both package suites on actual Julia 1.12.7/1.13.2, formatting and all MPI
ranks. [Job-log evidence](CI_RELIABILITY_PLAN.md#acceptance-run-at-240f77c).
