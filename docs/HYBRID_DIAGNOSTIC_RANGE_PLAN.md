# Hybrid diagnostic range follow-up

Status: the kinetic correction is pushed as `850fe46`; its exact-revision
GitHub CI is green in all six jobs on actual Julia 1.12.7 and 1.13.2.
The shared field correction has completed native integration gates and is
pushed as `5df5c38`; all six exact-revision GitHub jobs pass on actual Julia
1.12.7 and 1.13.2. The CGL correction
is pushed as `63f2420`; all six exact-revision GitHub jobs pass on actual
Julia 1.12.7 and 1.13.2. Isothermal integration is under
isolated verification; the polytropic repair remains planned.
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

The kinetic repair is pushed as `850fe462e7bb50eacda919dd6a5427411b4d151a`;
[run 39102304110](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39102304110)
must finish green before GitHub acceptance is recorded. The 24 sampled ordinary configurations
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
The full electromagnetic solver/restart, pressure, concrete compiler and Aqua
suites pass 30,668 assertions per runtime across 88 testsets. Complete package
compiler trees match baseline (22 reports/14 groups on 1.12.7; 21/13 on 1.13.2),
with evidence in `broader-complete-evidence.json`. All original tests and thresholds
remain. Combined integration with the pushed kinetic repair has completed all
native gates; this revision includes the shared field repair, with exact-SHA
GitHub verification still required.
Pre-edit contracts, failed evidence and revised source hashes are retained in
`/tmp/hybrid-field-energy-v{1,2,3}-verification`. This repair does not cover
algebraic electron pressure generation or complete energy-budget composition.

## Joint quadratic integration

The isolated combined snapshot preserves both candidate implementations and the
pressure feature. Its 280 Julia files include 191 source/test files. A new
132-assertion joint-budget regression passes on Julia 1.13.2: individual velocity
squares, magnetic squares and pressure sums lose range while their complete
reservoirs and total remain finite. Both parents have completed their native
gates and full compiler comparisons. Both joint numerical lanes pass 13,907
assertions across 18 testsets. Both also pass 1,271 quality/compiler assertions
and 1,257 four-thread assertions. Complete package compiler trees match baseline
on both versions (22 reports/14 groups and 21/13, respectively). The proposed public
docstrings and changelog are separate from the frozen test source; formatting and
executable-AST comparison pass across all 191 source/test files. The promotion
changes only the magnetic block of the current integrator, preserving the pushed
kinetic implementation byte for byte outside that block. Exact-SHA GitHub CI
remains required after the push. Evidence:
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

A separate [mixed-precision gamma witness](POLYTROPIC_ENERGY_RANGE_PLAN.md)
passes 40 assertions per native runtime and confirms that narrowing a stored
polytropic gamma into exactly one selects the wrong diagnostic branch. The
complete polytropic repair must address branch selection as well as power,
ratio and integral range.

The CGL repair below reuses the existing exact radical finisher for the complete
rational-plus-square-root integral. It preserves the zero-field singularity
convention and ordinary arithmetic, handles signed/nonfinite inputs explicitly,
and checks both fallback selection and final rounding. Polytropic powers and isothermal
logarithms need their own bounded evaluators and termination arguments; their
repair must not be described as complete merely by fixing materialized-pressure
sums. Mixed closure/grid precision is also part of that contract. The first
isolated CGL candidate reuses
the existing radical finisher and passes 944 focused assertions on each runtime,
including independent integrals, mixed precision, ordinary output/allocation and
operation-error enclosures. Broader resource controls keep it unaccepted:
16 witness assertions per runtime establish 1,676,808/1,684,736 B allocations
for ordinary 256-cell perpendicular field nulls / zero-density cells, versus
zero in the original owner. On Julia 1.13.2, its real 4,194,304-cell positive
grid also returns the correct unchanged value but allocates 27,682,407,880 B
and takes 714.96 s, versus 0 B and 0.304 s in the original reduction. The two
output/immutability assertions pass; the resource acceptance criterion fails.
The full measurement and source-hash check are retained in the v1
`large-grid-resource1132.log` and `resource-rejection.json`.

A separate v2 preserves the ordinary reduction, exact owner and signed error
certificate. It treats identically zero cells separately and adds a positive-term
range bound without the signed certificate's count limit. New checks include
zero-cell allocations, complete overflow midpoints, subnormal tails, mixed-precision
cancellation and a real 4,194,304-cell grid. Both Julia 1.12.7 and 1.13.2 pass
1,434 focused assertions. Separate measurements preserve the original large-grid
value and allocate 0 B on both versions: candidate/original elapsed samples are
0.768/0.133 s on 1.12.7 and 0.745/0.292 s on 1.13.2. These shared-host samples
show remaining certificate overhead, not an optimization limit. Source hashes
confirm the signed certificate and exact owner are unchanged from v1. Both
runtimes also pass 1,244 solver/multispecies/concrete-compiler/Aqua assertions
and 990 four-thread checks. Fresh full compiler-report trees match baseline
on both versions (22 reports/14 groups on 1.12.7, 21/13 on 1.13.2). Combined
integration has completed all native gates after the joint quadratic gates passed.
Pre-edit contracts, source hashes, `focused-acceptance-evidence.json`,
`large-resource-results.json` and retained failed measurements are in
`/tmp/hybrid-cgl-energy-v{1,2}-verification`.

The isolated CGL/quadratic integration snapshot contains 284 Julia files,
including 195 source/test files. Removing only its added CGL block recovers the
quadratic parent's integrator byte for byte; the CGL block and three CGL
regression files are byte-identical to v2. A new joint-budget regression covers
perpendicular, parallel and combined CGL reservoirs alongside particle and
magnetic energies, with scalar, collection and collection-view callers. Its
focused execution passes all 402 assertions on Julia 1.13.2. The combined native
suites pass 15,743 numerical, 1,487 quality and 2,247 threaded assertions per
runtime. Complete compiler trees match baseline on both versions (22 reports/14
groups on 1.12.7 and 21/13 on 1.13.2). Assembly scope and source hashes are
in `/tmp/hybrid-cgl-quadratic-integrated-v1-verification`. A separate public
documentation overlay passes formatting and executable-AST comparison across
all 195 source/test files. The CGL delta is pushed as `63f2420`;
all six exact-revision GitHub jobs pass. The promotion scope is limited
to the CGL block, public documentation and four regression files, preserving
the pushed kinetic and field implementations. Acceptance and current-main scope
are recorded in `acceptance-complete.json` and `promotion-current-main-scope.json`.


A separate isothermal cold-evaluation prototype is being checked in
`/tmp/hybrid-isothermal-range-design-v1`. It keeps complete rational factors,
encloses logarithms without MPFR, and handles exact zero using the dyadic
structure of stored densities. Its independent alternating-series and range tests pass 243 assertions on
each Julia 1.12.7 and 1.13.2. This verifies the isolated cold evaluator's tested
cases; 336 independent-bound and immutability assertions also pass under four
threads on each runtime. It does not establish a complete public repair. An
ordinary-path prototype uses outward bounds derived from the logarithm's
series. Both runtimes pass 468 numerical/ordinary-compatibility assertions,
and eight ordinary configurations retain their bits with 0 B allocations.
Both resource gates nevertheless fail (40 passes, one failure): a
4,194,304-cell Float32 weak wave has an uncertified interval spanning
-162.46265 to 208.69566 around the ordinary result 0.6159884. The large cold
fallback was not executed; no allocation or latency claim is made for it.
A further 71 witness assertions per runtime establish 19 incorrectly accepted
intermediate-underflow cases despite normal complete integrals; independently
scaled ordinary controls round correctly in those cases.

The [ENERGY-006 accuracy finding](ISOTHERMAL_ENERGY_ACCURACY_PLAN.md) also proves
that the existing public large-wave result differs from its independent
complete-integral enclosure by more than 16.08%. Consequently a positivity-only
certificate retaining that value would be insufficient. The corrected uniform-
density witness additionally establishes 1.38%–3.97%
Float32 errors outside the near-unit box. Both the bounded near-unit reduction
and a general accuracy repair are required alongside the range guards. The original
cold evaluator remains separately tested, but the public repair, broader
gates and polytropic implementation remain open.

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
