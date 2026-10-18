# Hybrid diagnostic range follow-up

Status: the kinetic correction passes local acceptance and is included in this
revision; exact-revision GitHub CI is pending. Field and CGL corrections remain
under isolated verification. Other algebraic-electron repairs remain planned.
This work is separate from the evolving-pressure feature and its stored-pressure
correction.

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
`true` for exact-fallback eligibility. Its frozen source passes 12,289
numerical/resource/concrete-compiler/Aqua assertions and 585 four-thread checks
on each native runtime. All 24 ordinary cases retain exact outputs and zero
allocations. Sampled latency ratios remain above baseline; the complete
measurements are retained in `resource-latency-summary.json` within the v3
verification directory. No performance-optimality claim is made.
Complete compiler-report trees match the baseline on both versions. The v3
pressure integration also passes its numerical and concrete compiler gates.
However, an expanded ordinary-input witness finds that Float32 weights `1e-10`
and `1e-20` change output bits and allocate 1,020,712/1,034,976 B for 256
particles, versus zero in the original reduction. This keeps v3 unaccepted.
The witness passes 12 assertions per runtime in `normalization-resources*.log`;
its assertions establish the regression rather than approve the candidate.

A separate v4 retains the complete evaluator and adds a checked normal-product
error enclosure, including signed population underflow. Its frozen 189 Julia
files pass 12,943 numerical/resource assertions on each runtime. These include
576 small-weight output/allocation checks and 1,208 independent rational
operation-enclosure checks. All 24 original resource cases still retain their
output bits and zero allocations. The broader pressure/budget, multispecies,
concrete compiler and Aqua gates pass 2,235 assertions per runtime; four-thread
checks pass 715 each. Complete package-report trees match the baseline:
22 reports/14 groups on 1.12.7 and 21 reports/13 groups on 1.13.2. The final
promotion's documentation preserves executable ASTs across all 189 source/test
files, and all five changed Julia files pass JuliaFormatter 1.0.62.

The kinetic repair is accepted locally and included in this revision; its
exact-revision GitHub run remains required. The 24 sampled ordinary configurations
retain output bits and zero allocations, as do the added small-weight checks.
Measured latency for 256–4,096 particles is 2.00–5.48× the old reduction on 1.12.7
and 1.81–3.45× on 1.13.2 in this shared-machine run. The correctness repair is
not a performance-optimality result; reducing certificate cost while preserving
these gates remains open. Full samples are in the v4 `resource-latency-summary.json`.
Evidence: `/tmp/hybrid-kinetic-range-v{1,2,3,4}-verification` and
`/tmp/hybrid-kinetic-range-integrated-v3-verification`.

## Field overflow boundary — ENERGY-005

**VERIFIED:** a second 16-assertion witness passes on each native runtime.
Hybrid magnetic, legacy electromagnetic and general electromagnetic diagnostics
all return a finite value below `floatmax(T)` when small positive squares lift
the exact complete integral past its nearest-even overflow midpoint. This
occurs for both Float32 and Float64. Commands, terminal statuses and the field
owner's byte identity with `013a70b` are recorded in
`/tmp/hybrid-legacy-energy-range-20261018/field-overflow-boundary-results.json`;
the neighboring logs retain each observed and exact result.

The field candidate reuses `_em_field_energy_exact/finish` for hybrid magnetic
energy and adds a conservative positive upper bound to the shared owner.
Its first version passes 438 existing range assertions but fails two ordinary
allocation controls (80 B instead of zero). Allocation profiling locates tuple
allocations in the grouped electric/magnetic component traversal; general
EMPIC remains allocation-free. A separate v2 replaces that traversal with two
typed loops. Its new test module omitted an import, causing 18 test errors;
v3 corrects the test and retains identical production code. The v3 focused
suite passes 880 assertions on each runtime. Its 42 ordinary configurations
pass 126 output/resource checks with zero allocations, and four-thread checks
pass 542 assertions per runtime. Measured latency is retained without an
optimality claim. The initial resource harness lacked a direct dependency;
corrected private environments resolve without changing the dependency graph.
Full electromagnetic solver/restart, pressure, concrete compiler and package
report comparisons remain in progress. All original tests and thresholds remain.
Pre-edit contracts, failed evidence and revised source hashes are retained in
`/tmp/hybrid-field-energy-v{1,2,3}-verification`. This repair does not cover
algebraic electron pressure generation or complete energy-budget composition.

## Joint quadratic integration

The isolated combined snapshot preserves both candidate implementations and the
pressure feature. Its 280 Julia files include 191 source/test files. A new
132-assertion joint-budget regression passes on Julia 1.13.2: individual velocity
squares, magnetic squares and pressure sums lose range while their complete
reservoirs and total remain finite. The full joint acceptance lanes wait for the
separate candidates' completed gates and report comparisons. The proposed public
docstrings and changelog are separate from the frozen test source; formatting and
executable-AST comparison pass across all 191 source/test files. No combined
repair has been promoted. Evidence:
`/tmp/hybrid-quadratic-energy-integrated-v1-verification`.

## Algebraic closure range — remaining ENERGY-004 scope

**VERIFIED:** 72 additional witness assertions pass on each Julia 1.12.7 and
1.13.2 against closure code and electron-diagnostic code byte-identical to
`013a70b`. The 18 Float32/64 cases return zero or infinity although their complete
integrals are finite, normal and nonzero. This confirms failures before the
pressure integral, beyond the previously observed finite-pressure sum overflow:

- Polytropic powers overflow or underflow before multiplication by `pe0`;
  the density/reference-density ratio can overflow before a square root whose
  exact result is finite. The witnesses include gamma 2 and gamma 1/2.
- Isothermal pressure overflows or underflows before the logarithm and cell
  volume can produce a finite complete free energy, including negative energy
  below unit density.
- CGL field-component squares overflow or underflow before the magnitude,
  corrupting both perpendicular and parallel internal energies.

The polytropic and CGL controls use exact rational expressions, with a single
nonzero field component making the CGL norm exact. Isothermal controls bound
`ln(2)` using a positive rational `atanh(1/3)` series and its geometric remainder;
no MPFR transcendental result is used as proof. Assertions require the observed
value to lose range, so ordinary floating rounding cannot satisfy the defect
witness. Commands, source comparisons, observations and terminal statuses are
in `/tmp/hybrid-algebraic-diagnostic-range-20261018` (`verified-results.json`).
These are confirmed defects, not accepted repair tests.

The next algebraic implementation should first address CGL's complete rational
plus square-root integral, assessing reuse of the existing exact radical
finisher. Preserve the current zero-field singularity convention and ordinary
arithmetic, handle signed/nonfinite inputs explicitly, and certify the fallback
selection as well as its final rounding. Polytropic powers and isothermal
logarithms need their own bounded evaluators and termination arguments; their
repair must not be described as complete merely by fixing materialized-pressure
sums. Mixed closure/grid precision is also part of that contract. No algebraic repair has been promoted. An isolated CGL candidate now reuses the
existing radical finisher and is being checked against independent integrals,
ordinary outputs and an operation-error enclosure; it is not yet accepted. Its
pre-edit contract and source are in `/tmp/hybrid-cgl-energy-v1-verification`.

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
