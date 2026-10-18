# Isothermal energy accuracy — ENERGY-006

Status: confirmed on Julia 1.12.7 and 1.13.2; isolated reduction prototype under
verification. Public repair remains unimplemented.

**VERIFIED:** the public isothermal diagnostic returns `0.6159884f0` for a
4,194,304-cell Float32 density wave of amplitude 0.001 about density one,
temperature `0.7f0` and unit cell spacing. An independent enclosure of the
complete stored-input integral is approximately `0.7340845405043187`, with
width `4.374e-17`. The relative error exceeds 16.08% on both runtimes. The
diagnostic and closure source are byte-identical to `850fe46`.

The oracle accumulates the first six dyadic density-deviation moments exactly
in Int128, with checked input/count bounds proving the accumulators fit. It
uses the expansion of `(1+d)log(1+d)` and an explicit geometric bound on all
remaining terms. It uses no floating logarithm or MPFR oracle. Each native
witness passes three assertions; commands, terminal statuses, source hash,
moments and numeric enclosures are retained in
`/tmp/hybrid-isothermal-range-design-v1/large-accuracy-witness{112,1132}.log`
and `large-accuracy-witness-results.json`. A separate 56-assertion cross-check
passes on Julia 1.13.2, comparing the Int128 moments with BigInt sums and the
polynomial/remainder with independent alternating-logarithm bounds
(`accuracy-oracle-crosscheck1132.log`).

This changes the acceptance contract for the isothermal follow-up. Retaining
the original bits cannot be a success criterion for this case. A certificate
that only proves the exact result is normal and finite would not repair this
accuracy defect. The previously considered density-box positivity bound is
insufficient on its own. The isolated range prototype is also unaccepted:
both native resource gates fail at this grid size, and 71 witness assertions
per runtime establish 19 incorrectly accepted intermediate-underflow cases.

1. Preserve those failures and the independent exact moment oracle. Add
   explicit guards for pressure, pressure-times-logarithm and spacing products;
   a normal final value does not excuse a damaged intermediate.
2. Develop a bounded near-unit entropy reduction, using
   `n*log(n) = (n-1) + h(n)` with the linear part summed exactly and the
   nonlinear part enclosed without cancellation. Assess fixed-width integer
   bounds for a scaled polynomial and its full analytic remainder. A result
   may be returned only when the complete interval rounds unambiguously.
   Retain the complete rational/logarithm owner for unresolved inputs.
3. Verify arithmetic overflow bounds, signed reductions, exact zero, mixed
   precision, volume/temperature range and rounding boundaries independently.
   Test weak waves, offsets on either side of unit mean, alternating cells,
   small grids and large grids. Corrected values may intentionally differ from
   the old result; accuracy assertions must use independent complete integrals.
4. Measure allocations and latency, including scaling with grid size, before
   choosing the production implementation. Retain ordinary compatibility checks
   outside the corrected behavior. Complete native solver/compiler/threaded
   gates and exact-revision GitHub CI before claiming acceptance.

The existing [diagnostic range plan](HYBRID_DIAGNOSTIC_RANGE_PLAN.md) remains
open for public isothermal integration, polytropic evaluation and complete
budget/inertia composition.

The isolated `/tmp/hybrid-isothermal-nearunit-v1` prototype uses exact linear
sums, scaled fixed-point polynomial bounds and an explicit full remainder.
Both native versions pass 1,902 focused assertions, including independent
rational checks of fixed-point coefficients/products and complete-integral
comparisons against alternating logarithm bounds. Each runtime also passes
28 eighth-moment accuracy/resource assertions across Float32/64 and 256,
4,096 and 4,194,304 cells. Float64 large-state controls use exactly promoted
Float32 densities, so the independent integer-moment oracle remains applicable.
The large Float32 result is corrected to `0.73408455f0`.

The grid polynomial kernel allocates 0 B; complete calls allocate 9,896–10,440 B
independently of grid size in these measurements. That exceeds the old small-grid
allocations while avoiding the old 16/32 MiB pressure arrays at the largest size.
Shared-host large-grid samples take 3.47–4.86 s, versus 0.038–1.03 s for the
original calculation. No optimality or overall performance acceptance is
claimed. Allocation/latency optimization, operation-range guards, public
integration and broader gates remain open.

An isolated v2 replaces the rational final interval with bounded dyadic Int128
arithmetic and direct nearest-even integer rounding. Each native runtime passes
5,823 focused assertions, including 3,921 new exact-operation, rounding-boundary
and full-mantissa checks. The resource oracle separately passes 144 comparisons
against BigInt moments and independent logarithm bounds on Julia 1.13.2.

The v2 Julia 1.12.7 resource run records zero whole-call allocations in all six
cases, but its inherited strict ratio comparison fails twice because it becomes
`0 < 0`. The run is retained as 26 passes/two failures; the 1.13.2 resource lane
did not run after that failure. A separate v3 changes the resource contract to
require zero bytes in every case, retains the independent accuracy assertions,
and tests Int64 per-cell arithmetic only where the existing product bound leaves
three sign/carry bits. Whole-grid sums remain Int128. The v3 passes 5,991 focused assertions and 34 independent accuracy/resource
assertions on each native runtime. All six measured whole calls per runtime
allocate 0 B. Its 24 concrete compiler checks pass per runtime; an additional
threaded harness records every worker result in the parent testset on both versions.
That threaded harness passes 1,226 assertions per runtime across 36 fixtures,
with all four worker thread IDs observed.

A paired Julia 1.13.2 measurement passes 28 correctness/allocation assertions.
For 4,194,304 Float32 cells, three alternating-order samples give a candidate/prior
median runtime ratio of 0.398 (1.274 s versus 3.198 s). Float64 retains the Int128
kernel. These are shared-host samples, not an optimality result. The candidate
remains slower than the original inaccurate reduction in the separate resource
measurements. Public integration remains open. Evidence is retained in
`/tmp/hybrid-isothermal-nearunit-v{2,3}`.

The separate `/tmp/hybrid-isothermal-range-design-v2` probe adds explicit
promoted-pressure, grid-conversion, pressure-times-logarithm and contiguous
spacing-product guards. All 96 assertions for the original intermediate-underflow
family pass on Julia 1.12.7. Its first added test incorrectly required exact final
rounding for five ordinary cases outside the corrected domain; the classification
log verifies all five match the existing public routine bit for bit. That failed
run (607 passes/five failures overall) is retained. The separate corrected harness
preserves the prewritten ordinary-compatibility contract and requires complete
rounding for unsafe cases. It passes 660 assertions on each native runtime.
No prototype source was changed to repair that test contract.

The combined `/tmp/hybrid-isothermal-integrated-design-v1` snapshot now routes
finite near-unit states through the accurate bounded reduction, with unresolved
states sent directly to the complete logarithm evaluator. Other inputs use the
range-checked ordinary path and existing nonfinite propagation. Its 11 Julia
files are frozen; the range/cold sources and bounded arithmetic body match the
tested parents. The entry split avoids a duplicate density-statistics scan.
Combined numerical and resource gates pass on both runtimes, retaining zero
whole-call allocation in all six cases. The Julia 1.12 concrete compiler gate
reports an undefined-variable possibility in BigInt hashing reached through the
cold logarithm dictionary; its failure is retained. Julia 1.13.2 passes that gate.
This combined candidate is unaccepted. A bounded-key representation is under
investigation; no report filtering is permitted. Public dispatch is unchanged.
