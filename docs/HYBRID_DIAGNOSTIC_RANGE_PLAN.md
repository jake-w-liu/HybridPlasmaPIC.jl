# Hybrid diagnostic range follow-up

Status: confirmed defects. A kinetic correction is under isolated verification;
magnetic and algebraic-electron repairs remain planned. This work is separate
from the evolving-pressure feature and its stored-pressure correction.

## Evidence

**VERIFIED:** the public nonrelativistic `kinetic_energy`, `magnetic_energy`
and polytropic `electron_internal_energy` lose finite complete integrals on
Julia 1.12.7 and 1.13.2. A zero magnetic field can also produce NaN when the
cell-volume product overflows. The independent rational witness passes 16
assertions on each runtime in
`/tmp/hybrid-legacy-energy-range-20261018/witness{112,1132}.log` (both exit 0).

The tested integrator file is byte-identical to `main` at `e339350`;
`/tmp/hybrid-legacy-energy-range-20261018/evidence.json` records its SHA256,
command and revision comparison. The checked owners in
`src/integrators/extrapolated_leapfrog.jl` are lines 555, 839 and 868:

- Kinetic energy forms velocity squares before combining the mass. With
  Float64 mass `1e-300` and speed `1e200`, it returns Inf instead of the finite
  complete result `5e99`.
- Magnetic energy forms squared fields and their sum before applying volume.
  A Float64 field `1e200` over length `1e-200` returns Inf instead of `5e199`.
- Polytropic energy with two maximum finite pressures, gamma 3 and unit domain
  length returns Inf instead of half the maximum finite value.
- A zero field on a 2×2×2 grid with side lengths `2e103` returns NaN instead of
  zero. The grid spacings and field entries are finite.

Float32 witnesses reproduce each class. These results concern numerical
evaluation of diagnostics, independently of the continuum energy-model balance.

**VERIFIED collection and boundary cases:** 20 additional assertions per runtime
pass in `collection-boundary-witness{112,1132}.log` (both exit 0; commands and
hashes in `collection-boundary-results.json` under the same evidence directory).
For two finite populations with opposite weights, each separately evaluated
energy overflows; adding a third unit-energy population yields NaN although the
complete rational total is one. Scalar and collection sums can also retain
`floatmax(T)` when small positive terms bring the complete exact value to the
nearest-even overflow midpoint. The repair must therefore evaluate the requested
whole collection and check finite-looking overflow boundaries.

## Kinetic candidate acceptance

The first isolated kinetic candidate passes 2,646 exact-range, population/view
and nonfinite checks on Julia 1.13.2, but fails its ordinary-resource gate:
40 passes and eight failures. Four signed Float32 cases change their ordinary
result and allocate 64,936–16,433,960 bytes. It is not accepted. Frozen source,
commands and results are in `/tmp/hybrid-kinetic-range-v1-verification`.
A separate v2 retains the complete rational fallback, adds a forward-error
certificate for ordinary signed sums, and separates Boolean checks from the
ordered floating reduction. Its range, certificate and unchanged resource gates
pass 11,159 assertions on each runtime; all 24 measured ordinary configurations
retain their output bits and zero allocations. Measured latency remains higher
than the baseline. A concrete compiler gate exposes a possible `Missing`
condition for abstract-storage collection types. A separate v3 requires known
`true` for exact-fallback eligibility; its complete numerical, compiler, threaded
and package-report gates are running. No kinetic repair has been promoted.
Evidence: `/tmp/hybrid-kinetic-range-v{1,2,3}-verification`.

## Repair contract and sequence

1. Freeze ordinary and exceptional behavior at the public owners. Include
   scalar and collection particles, signed weights, cancellation, empty
   populations, supplied views, and existing nonfinite propagation. Assess
   collection accumulation separately from each population's returned energy.
2. Repair quadratic kinetic and magnetic integrals first. Retain ordinary
   order and allocations; detect range-sensitive products, sums and every cell
   spacing association. Eligible finite exceptional inputs should use a streamed
   complete rational expression and the existing final-rounding owner. Avoid
   full-sized field copies and independent rounding of population contributions
   when their complete sum is the requested result.
3. Define algebraic electron-energy semantics before changing its evaluation.
   A finite materialized pressure integral and pressure lost inside the closure
   are distinct cases. Assess polytropic powers, isothermal logarithms and CGL
   magnitudes/pressures with independent witnesses; their complete range behavior
   has not been established by the current 16 assertions. Reusing the stored
   pressure helper alone cannot establish correctness when pressure generation
   itself overflows or underflows.
4. Specify a verified bound and termination policy before adding any irrational
   cold evaluator. Retain true output overflow, underflow and exact cancellation;
   do not treat a larger precision, a finite retry count or a nonfinite temporary
   as proof. The separate MPFR concurrency finding rules out assuming that a
   process-wide precision lock makes concurrent MPFR calls reliable.
5. Verify complete energy-budget composition, including the existing inertia
   reservoir, after the individual components. The stored-pressure feature's
   component checks do not establish the range of every other budget reservoir.

Acceptance requires exact rational rounding-cell tests for quadratic expressions;
an independent certified oracle for any nonrational extension; both native Julia
versions; ordinary output/allocation comparisons; concurrency; concrete and
complete package compiler comparisons; affected physical/conservation tests;
and exact-revision green GitHub CI after each pushed repair. Do not alter physical
forces, closure laws, diagnostic conventions or validation tolerances to fix
arithmetic range loss.
