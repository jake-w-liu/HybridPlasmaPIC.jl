# Hybrid diagnostic range follow-up

Status: the kinetic correction is pushed as `850fe46`; its exact-revision
GitHub CI is green in all six jobs on actual Julia 1.12.7 and 1.13.2.
The shared field correction has completed native integration gates and is
pushed as `5df5c38`; all six exact-revision GitHub jobs pass on actual Julia
1.12.7 and 1.13.2. The CGL correction
is pushed as `63f2420`; all six exact-revision GitHub jobs pass on actual
Julia 1.12.7 and 1.13.2. Isothermal and polytropic repairs are also pushed;
their current exact-SHA results are tracked in [CI acceptance](CI_RELIABILITY_PLAN.md).
The inertia repair now has all six exact-SHA CI jobs green; complete budget
composition remains open.
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

## Complete budget composition and inertia — current ENERGY-004 scope

**VERIFIED:** both Julia 1.12.7 and 1.13.2 pass 316 defect-witness assertions
against immutable `99760f0193a40c3a4bdb18ba9d1013e11a49d0ea` source:

- 282 assertions cover 1D–3D, Float32/64, scalar/collection/view populations and
  polytropic, evolving-pressure and CGL closures. Opposite overflowing kinetic
  and magnetic reservoirs cancel exactly, leaving unit electron energy, yet
  `energy_budget.total` is NaN.
- Sixteen assertions show three positive half-minimum-subnormal reservoirs
  separately round to zero. Their complete total rounds to two minimum
  subnormals, while the public total is zero.
- Eighteen assertions cover electron-inertia square overflow, square underflow
  and zero current with overflowing cell volume. The public inertia result is
  infinity, zero or NaN although an independent rational expression from the
  same finite spectral-current output is finite and defined.

Source hashes, immutable package/environments, commands, observed values and log
hashes are in `/tmp/hybrid-budget-composition-review-20261019/verified-results.json`
and `source-hashes.json`. At this exact revision, `conservation.jl:77–81`
materializes the squared-current reduction before the coefficient/volume, and
`:88` adds already-rounded reservoirs. These are distinct arithmetic failures;
the witness does not claim new continuum conservation physics or exact spectral
differentiation.

### Inertia repair pre-edit brief

- Mode: repair of the existing electron-inertia diagnostic reservoir.
- Requirement: retain finite computed-current values, the validated grid-precision
  `de2` and each cell spacing through range-sensitive integration. Preserve
  normal output bits, allocation, validation order, zero-`de2` bypass and
  generic/nonfinite behavior. Complete-total composition remains the next step.
- Baseline/cause: the 18 inertia witnesses above and the checked reduction at
  `conservation.jl:77–81` establish premature squaring/summing/volume range loss.
- Owner: `energy_budget` in `src/diagnostics/conservation.jl`; existing field
  exact arithmetic in `src/models/full_pic.jl` supplies the common square sum.
  Keep `curl!` and its three output arrays unchanged.
- Strategy: extract the unchanged ordered inertia reduction into one private
  owner, add normal-product/volume and positive-upper-bound checks, and retain
  the complete rational expression when those checks cannot certify range.
  Share the streamed exact square sum with electromagnetic field energy.
- Compatibility: no public signature or returned-field changes; no model-law
  changes, arbitrary precision cap, external cache or alternative diagnostic.
- Verification: original witnesses, subnormal/overflow neighbors and ties,
  independent rational rounding cells, ordinary bitwise/zero-allocation
  controls, public budgets and shared field-energy regressions, invalid and
  nonfinite cases, four-thread ownership, complete compiler comparisons,
  formatting, then a separate pushed exact-SHA GitHub run.
- Wrong-patch controls: scaling a rounded magnetic energy fails the coefficient
  recovery cases; checking only nonfinite output misses underflow and finite
  overflow neighbors; changing curl/workspace ownership fails public bitwise,
  allocation and concurrent-input checks. None is an acceptable shortcut.

The isolated inertia component passes 456 assertions per runtime in
`/tmp/hybrid-inertia-energy-v1`, preserving ordinary bits and zero allocation.
Canonical V1 then passes 2,738 focused assertions and 1,987 complete public
budget/inertia/multispecies/stored-pressure assertions per runtime. Complete
compiler trees match baseline (22/21). Its first field lane omits the prerequisite
CGL oracle test; the corrected lane retains every selected test and adds that
prerequisite. Its first migration harness loses indentation in triple-quoted
match strings; the corrected harness asserts each exact match and passes all
77 checks. Original failures are retained in the canonical verification directory.

V1 remains **unaccepted on resources**. Both native runtimes measure
1,324,357,672 B (Float32) and 1,327,504,552 B (Float64) cumulative allocation
for one zero-field/current 64³ budget call with overflowing cell volume. Its
numerical result is zero. These are full-call allocation measurements, not
peak-memory measurements. The six zero-state samples per runtime are retained
in `costs{112,1132}.log` under
`/tmp/hybrid-inertia-energy-public-v1-verification`.

The follow-up resource contract certifies genuinely zero square inputs using
both the unchanged reduction and the existing nonzero-square minimum sentinel.
It must not classify an underflowed nonzero square, a nonfinite field or an
uncertified coefficient as zero. Preserve ordinary signed zero and generic
behavior. Add the certificate to the shared field finisher and the inertia
owner before exact accumulation. For 16³ and 64³ zero states, require public
allocation at most 110% of the corresponding ordinary zero-state call, retain
all V1 numerical/ordinary checks, and repeat complete compiler, field, budget
and concurrency gates. This is a measured resource correction within the same
unpromoted repair, not a waiver of V1's costs.

**VERIFIED allocation attribution:** the sampled public-call profile records
7,282 allocation samples in `_field_square_sum_exact`, at V1 `full_pic.jl:1484`,
with `BigInt` and GMP rational allocation types. The initial reporting script
has a Julia global-counter scope error; `zero_profile_fixed.jl` corrects that
harness, retaining the failed script/log. Profile command and output are in
`/tmp/hybrid-inertia-energy-public-v1-verification/zero-profile-fixed112.log`.
All V1 numerical gates eventually pass: 2,738 focused, 1,987 budget, 12,689
shared-field and 666 measured-cost assertions per runtime. The high allocation
keeps V1 unaccepted; a numerical pass does not waive that resource result.

Canonical V2 retains every earlier assertion and adds zero-state allocation and
false-zero controls. Its measured 64³ zero-state allocation is 3,145,968 B
(Float32) / 6,291,696 B (Float64) on both runtimes, compared with V1's
1.32–1.33 GB. Both measured cost suites pass 666 assertions; paired ordinary
latency ratios for `de2 > 0` are 0.87–1.18 / 0.89–1.19. All ordinary output and
allocation comparisons pass. Complete V2 local acceptance is still being
collected in `/tmp/hybrid-inertia-energy-public-v2-verification` before promotion.
The complete-total composition repair remains a separate required step.

V2 local acceptance is complete on both Julia 1.12.7 and 1.13.2:

- 2,846 focused assertions across 11 testsets, including all V1 checks, 108 new
  zero-state assertions, explicit four-thread arithmetic ownership, independent
  public grids and concrete public-budget compiler calls.
- 1,987 assertions across 22 public budget/inertia/multispecies/stored-pressure
  testsets, plus 12,689 across 19 shared electromagnetic/hybrid field,
  kinetic, CGL, joint-budget and Aqua testsets.
- 666 measured-cost assertions per runtime. All 24 ordinary profiles preserve
  results and allocation. The declared 110% ordinary-zero allocation gate
  passes for both 16³ and 64³ states. Nonzero exceptional profiles through
  4,096 cells allocate at most 26,764,128 B; largest sampled median calls are
  37.6 ms / 76.2 ms. Larger nonzero cold-grid costs and peak memory are unmeasured.
- Complete compiler trees match baseline (22/21 reports, no added/removed
  groups), all 77 final migration checks pass, and 11 formatting AST checks
  plus the whole-repository read-only scan pass.

The 11 tested source/test paths are copied to main with matching hashes; held
radiation work is unchanged. The shared exact square sum has one owner. Inertia
keeps its original floating reduction/order and validation, uses the same
spectral current and three output arrays, and checks the complete coefficient
and volume only when range cannot be certified. Zero certificates preserve
signed ordinary zero, reject nonfinite factors, and distinguish nonzero squares
whose floating intermediates vanish or overflow. Public signatures and the
model equations are unchanged. Commands, final assertion inventories, complete
compiler trees, prototype failures and cost samples remain in the V1/V2
verification directories. Exact-SHA GitHub acceptance is required after push;
complete budget composition remains open.

The verified inertia change is pushed as `ffdb522aa858622a645c9ace104726ab54782493`.
All 16 committed blobs, remote main and held radiation hashes match
`/tmp/hybrid-inertia-energy-public-v2-verification/committed-evidence.json`.
Its [exact-SHA run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39150381077)
passes all six jobs: 502,711/502,712 assertions on actual Julia 1.12.7/1.13.2,
formatting and all MPI ranks. Full compiler trees remain unchanged at 23/22.
The exact-SHA evidence is in `/tmp/hybrid-ci-ffdb522-20261019`; no billing
exception is used.

### Complete-total evaluator pre-edit brief

- Mode: repair the total of the existing budget; first validate the complete
  expression independently of the ordinary-path certificate.
- Acceptance: combine kinetic, magnetic, closure and computed-current inertia
  terms before final rounding when their separate rounding loses range.
  Keep the four reported reservoirs, public signatures, validation order,
  spectral current, nonfinite/singular behavior and ordinary output bits.
- Baseline: the retained 298 composition assertions above establish false NaN
  and lost subnormal totals. The canonical `energy_budget` still adds four
  rounded values. The inertia correction does not repair that operation.
- Owners: `conservation.jl` composes reservoirs; `extrapolated_leapfrog.jl`
  owns exact kinetic/stored-pressure/CGL evaluation; the existing isothermal
  and polytropic energy modules own logarithmic/power integration. Share their
  mathematical reductions instead of adding parallel closure formulas.
- Strategy: expose the existing exact kinetic and stored-pressure rational
  reductions; let the canonical exact CGL, logarithm and polytropic finishers
  accept a rational offset. The offset is the complete quadratic contribution.
  Extend scaled interval summation to signed terms, retaining bounded binary
  exponents, strict endpoints and exact midpoint/zero decisions. Keep the
  current public component calls delegating to these same owners with zero
  offset. No global precision state, approximation cap or external cache.
- Integration boundary: the public budget keeps its ordinary evaluation and
  returned components. A separately verified total-range certificate must
  establish when the ordinary sum can be retained. Prototype success alone
  does not authorize replacing every ordinary call by expensive exact work.
- Compatibility: no new public diagnostic, model law, closure family or
  conservation claim. CGL's zero-field parallel singularity and unsupported
  generic/nonfinite input arithmetic remain visible. Evolving pressure keeps
  its validated grid-precision gamma; algebraic gamma keeps stored precision.
- Verification: original witnesses plus cancellation to both signed neighbors,
  exact zero, half-subnormal and overflow midpoints, huge stored gamma,
  irrational powers/logarithms, mixed closure/grid precision, density floors,
  collections/views/typed empties, input ownership and concrete compiler calls.
  Retain earlier component suites, independent rational/interval oracles and
  four-thread checks. Measure ordinary allocations and paired latency plus
  exceptional allocations before canonical promotion; compare complete
  compiler trees and formatting, then push and require exact-SHA CI.
- Wrong-patch controls: summing rounded components fails the original witnesses;
  dropping small scaled terms fails midpoint-offset neighbors; a precision
  cap or tolerance-based zero fails huge-gamma exact cancellation and ties;
  unconditional rational replay fails ordinary/zero-state resource gates.

For logarithmic completion, clearing the rational log coefficients makes
`exp(d*L)` a positive rational. If nonzero `L` were rational, this would
contradict the Hermite–Lindemann consequence of the
[Lindemann–Weierstrass theorem](https://arxiv.org/abs/2306.14352).
This is a derivation from the theorem, checked October 19. The existing
odd-significand sign argument distinguishes zero log sums. Therefore a
nonzero log sum plus a rational offset cannot equal a rational rounding
midpoint or zero. For powers, the previously recorded positive-power
rationality derivation still applies after moving the rational offset to the
other side; all-rational sums require the existing exact equality decision.

The first isolated offset candidate is retained at
`/tmp/hybrid-budget-composition-v1`. It passes 912 signed-interval assertions,
then fails to terminate in the half-subnormal tests. Its selecting expression
chooses zero for the adjacent intervals `[0,qmin]` and `[-qmin,-0]`, so it
does not ask the exact equality owner about the actual midpoint. Sixteen
bounded witness assertions per runtime verify that selected-zero equality is
false and adjacent-midpoint equality is true. The current-session test tool
was interrupted; its source, frozen hashes, log and interruption record remain.
V2 prefers the adjacent boundary and uses the zero decision only for a strict
sign-straddling interval without an adjacent boundary. This is candidate work,
not a defect promoted to main.

The ordinary-total certificate will enclose complete reservoirs before allowing
the existing floating sum. Kinetic bounds reuse the absolute-operation probe
and its error owner. Positive quadratic/stored-pressure reservoirs use their
existing normal-operation or complete-rounding contracts, with an outward
rounding-factor allowance and one minimum subnormal for final underflow.
Positive CGL terms use the same allowance at the smaller input precision;
mixed-sign terms reuse the existing absolute CGL probe. Polytropic nonunit
gamma reuses the complete exponent bounds. Gamma one has a correctly rounded
component; other isothermal cases use the existing near-unit classification or
general logarithm enclosure. No certificate duplicates a closure's exact
formula or assumes the rounded components are exact.

For at most `k` normal rounding/inverse-rounding factors, a conservative
relative allowance is `4ku/(1-4ku)`, restricted to `ku <= 1/16` and rounded
outward, where `u=2^-p` and `p` is the smaller relevant precision. The exact to
computed ratio is between `(1-u)^k/(1+u)^k` and its reciprocal; bounding the
logarithm and then `exp(x) <= 1/(1-x)` gives this allowance. Signed reductions
need their absolute-sum bound instead. Nonfinite or unresolved intervals must
use complete evaluation when eligible; unsupported inputs retain their native
arithmetic. Only a structurally proved all-zero state may bypass a zero total's
subnormal uncertainty.

Public resource targets are unchanged ordinary allocations and a paired median
latency ratio at most 1.5 on the declared ordinary profiles, including all
closures and nonuniform density. This target allows the additional total-range
certificate; it is not a speedup claim. Zero states must retain the prior
64³ allocation bound. Record exceptional costs through 4,096 cells, preserve
all prior numerical tests and reject a numerically passing candidate that
violates these gates. Public acceptance is pending.

V2 completes on both runtimes with 3,968 offset/interval assertions, 328 complete
expression witnesses, 48,973 retained component assertions across 29 testsets,
and 5,580 initial certificate/numerical/allocation checks. The four reported
components remain unchanged. These results are recorded in
`/tmp/hybrid-budget-composition-v2/completed-counts.json`.

V2 is **unaccepted on latency**. Of 324 paired ordinary profiles per runtime,
95/91 exceed the 1.5 ratio target; the largest sampled ratios are 3.26/3.80.
Both zero-state suites retain the prior allocation exactly, including
3,145,968/6,291,696 B at 64³. The passing numerical/resource assertions do not
override the latency measurements in `ordinary-costs{112,1132}-summary.json`.

**VERIFIED repeated work:** instrumentation inserted only at the existing
general-logarithm probe entry counts one traversal in the baseline and two in
the V2 budget for Float32/64 arrays of 64, 512 and 4,096 cells. All 18 count and
output assertions pass on both runtimes (`probe-count{112,1132}.log`). This
attributes repeated logarithm work; it does not explain every slow profile.

V3 retains the canonical isothermal component result and its existing bounds
in an immutable state. The public scalar component reads its value; the budget
consumes its bounds. Algebraic `pe` rejection has one shared owner, and the
budget's existing grid validation remains before calculation. Other closures
still call the existing public diagnostic, preserving custom-closure dispatch.
Same-sign kinetic contributions use the proved rounding-factor allowance;
mixed-sign contributions retain the existing absolute-operation probe. The
minimum-runtime instrumented count is back to one for all six profiles. The
first V3 edit script completed these changes but failed an unnecessary final
instrumentation-marker replacement; the original script and precise correction
are retained, with no test assertion removed.

V3's initial comparison uses the changed component implementation on both
sides, so it cannot establish the whole change's latency against pushed main.
Canonical public V1 therefore pins the baseline scalar/B validation, closure
selection and both isothermal reductions to `ffdb522`, with source hashes and
exact extraction recorded in `legacy-source-evidence.json`. It retains all
324 profiles and compares equal-sized alternating batches to reduce timer
granularity on submicrosecond calls. The original single-call measurements
remain evidence. The 1.5 target, numerical equality and allocation checks are
unchanged; the new public candidate is still unaccepted.

Canonical public V1 fails 37/28 assertions in its 8,748-assertion paired suite.
Public V2's full compiler comparison finds three new mixed-endpoint dispatch
groups. V3 defines adjacency on the left endpoint's floating lattice while
accepting either stored endpoint precision; both complete compiler trees return
to the unchanged 22/21 baseline. Its 64 endpoint/unsupported-precision checks
pass per runtime. The first comparison used a resolved `/private/tmp` root
instead of the actual manifest `/tmp` root; both original and corrected
comparisons are retained, with every report frame preserved.

Public V3 still fails the resource gate: 23/24 latency profiles exceed 1.5,
and 1/2 alternating-batch allocation equalities fail on Julia 1.12.7/1.13.2.
The full suites pass 8,724/8,722 assertions and fail 24/26. These are not all
latency failures; the batch allocation discrepancies remain unexplained.
Evidence: `/tmp/hybrid-budget-composition-public-v3/paired_costs{112,1132}.log`.

**VERIFIED large-grid replay:** with Float32 80³ unit-spacing cells, unit Bz,
zero temperature and no particles, the pushed owner allocates zero bytes;
public V3 allocates 1,331,206,192 B on both runtimes. Its magnetic operation
count is 1,536,010, exceeding the certificate's 1,048,576 limit and routing the
normal state into the rational replay. The corresponding warmed single-call
times are 3.56/6.03 ms versus 3.75/2.96 s; these are not median or peak-memory
measurements. Both 10-assertion witnesses pass; 32³ remains allocation-free.
Evidence: `large_grid_witness{112,1132}.log` in the same candidate directory.

Both corrected threaded suites pass 770 assertions with separate spectral
workspaces and shared physical inputs. The original fixture fails before the
concurrent calls because two snapshots destructure three entries as four;
the corrected driver changes those two destructuring patterns only. Original
failure logs, correction note, frozen hashes and corrected logs are retained.

### Scalable total certificate pre-edit brief

- Preserve the complete-expression evaluator, old component output bits,
  validation, supported storage/precision and all existing acceptance gates.
  The current candidate is isolated and must not be promoted on numerical
  checks alone.
- Replace the count cutoff in the same-sign rounding enclosure with a wider
  power-of-two bound when the existing tight bound is unavailable. For
  `u=2^-p <= 1/4`, `log((1+u)/(1-u)) <= 2u/(1-u) <= 8u/3 <= 4u*log(2)`.
  Thus any count `k` of the already certified normal rounding/inverse factors
  gives a ratio in `[2^-r,2^r]`, with `r=ceil(4*k*u)`. Compute `r` by integer
  ceiling division, retaining outward rounding and the final minimum-subnormal
  allowance. Unbounded intervals still require complete evaluation. This
  extends the same certificate; it introduces no second energy formula.
- Retain the kinetic owner's already computed contribution signs and count
  in an immutable internal result. Public kinetic energy reads its value;
  the budget consumes the same state. Keep the ordinary reduction and all
  finite/nonfinite decisions in one owner. Preserve generic-precision and
  scalar/collection behavior; do not infer sign from the rounded total.
- Validate the new ratio bound independently with exact rationals across the
  old cutoff, signed values, subnormals and overflow edges. Preserve all
  offset, composition, component, threaded, compiler and ownership assertions.
  Extend ordinary profiles to 80³/128³ grids and large particle counts, all
  closures, nonuniform densities, mixed coefficient precision and signed
  reservoirs. Same allocations and paired median ratio <=1.5 remain required.
- If large signed cases still trigger expensive replay, measure those cases
  before deciding on a wider native refinement. Do not add another formula,
  weaken a gate or attribute failures to hardware without evidence.

Public V4 passes the 5,580 initial certificate/value/allocation checks on both
runtimes and 631 independent scalable-bound/large-particle checks per runtime.
The directed oracle was corrected to round its `log(2)` denominator downward;
the original driver and correction are retained. Its paired suite already
records small-grid latency failures, so V4 is not accepted. The full paired
and large-grid lanes continue on its frozen source.

The next isolated optimization simplifies the tight rounding allowance:
when `4ku <= 1/4`, `4ku/(1-4ku) <= 8ku`. The latter is an exactly represented
dyadic factor under the same count guard. One outward product bounds the
relative error; its successor is at least one minimum subnormal, covering
correctly rounded subnormal component results. Outward endpoint additions
remain. This removes the division and redundant outward operations without
loosening the 1.5 acceptance target. Inlining the total certificate and its
kinetic-state consumer is a performance hypothesis to test, not an accepted
speedup claim. Retain the count-extension proof, old regression suites, pinned
component baseline and all large-grid gates.

**VERIFIED candidate compatibility gap:** the original `ffdb522` budget body
accepts custom `electron_internal_energy` results with a different numeric
type from the grid. V4 ties the internal electron-state parameter to the grid
type and throws `MethodError`. The corrected 36-assertion witness on each runtime
verifies both old acceptance and candidate failure. Its first fixture
used a renamed generic component that bypassed custom dispatch; the corrected
fixture renames only the original budget binding. Both drivers and extraction
hashes remain in `/tmp/hybrid-budget-composition-public-v4`.

The V5 compatibility repair leaves public dispatch and summation order intact:
accept the internal state regardless of its value type, form the original sum,
and restrict the new certificate to matching Float32/64 component data.
Other custom results keep that sum. The regression must cover integer,
rational, other floating precision, nonfinite and keyword-bearing custom
results, plus the original built-in numerical/ownership tests. This is a
repair to an isolated candidate, not a shipped regression.

V4's 132 large-grid profiles on Julia 1.12.7 meet the 1.5 latency target but
retain two byte-equality failures (1,714 passes, two failures). Brief custom
closure and signed-range witnesses also ran while the unaccepted V4 cost
lanes were active; those timings are exploratory, not final controlled
acceptance. The later candidate's cost lanes are scheduled after V4 finishes.

**VERIFIED remaining signed replay:** an 80³ Float32 unit-spacing state with
`Bz=0.5`, `n=0.8`, `Te=0.5` and no particles has finite total `18300.203f0`.
V4's magnetic enclosure is approximately `[32000,128000]`; adding the electron
enclosure crosses zero. The complete evaluator is selected and allocates
1,328,191,872 B, versus zero for the pushed budget, on both runtimes. Six witness
assertions per runtime verify the selection, unchanged components and an
independent complete-integral enclosure. Single warmed calls take 3.95/3.31 s
versus 13.2/8.87 ms; these are not paired medians or peak-memory measurements.
Evidence: `/tmp/hybrid-budget-composition-public-v4/signed_large_witness{112,1132}.log`.

### Wider field refinement pre-edit brief

- Keep the fast certificate and all public component results. If finite
  Float32 reservoirs remain unresolved, refine only magnetic/inertia bounds
  with Float64 arithmetic before considering the complete rational replay.
  Reuse the existing electron and kinetic enclosures. Unresolved cases still
  take complete evaluation; unsupported arithmetic retains its current sum.
- Generalize the canonical square-sum owner to the two required accumulator
  types, Float64 and `Rational{BigInt}`; migrate every existing exact caller
  and direct quality test. Exact accumulation must continue to call
  `_exact_rational`, including its stored-precision BigFloat conversion.
  No compatibility shim or second quadratic formula is needed. Zero squares
  may be skipped because the accumulator starts at mathematical zero.
- Every finite Float32 square is a normal Float64 value. With at most three
  positive Float32 spacings, a Float32 inertia coefficient and an Int-sized
  grid, complete nonzero products lie above `2^-896` and below `2^837`.
  Thus the wider reduction's normal-rounding enclosure applies without
  intermediate range loss. Use its actual sequential operation count and
  outward conversion to Float32; retain the original computed spectral J.
- Share one total-range predicate across both stages. Require a finite
  ordinary total whose sign agrees with the certified complete enclosure;
  refinement must not accept a known wrong sign. This strengthens the range
  check without changing the declared accuracy of ordinary component sums.
- Verify independent rational bounds, all zero fields, extreme finite Float32
  exponents and anisotropic spacings, scalar/collection/view ownership and
  unchanged stored-BigFloat conversion. Preserve original field/inertia tests
  after the private-symbol migration. Add signed large-grid costs, retained
  small/large profiles, full compiler comparisons and formatting. Neither
  the latency target nor allocation acceptance may be relaxed to ship this.

The V5 small-profile screen passes all 1,620 assertions per runtime, including
the unchanged 1.5 latency target and byte comparisons. Its full acceptance is
still pending. A separate same-owner allocation probe fails 2/1 of 4,800
assertions on Julia 1.12.7/1.13.2: identical pinned-budget calls sometimes record
extra malloc and pool allocations as well as extra bytes. Every observed
mismatch includes a GC pause. Therefore simply replacing bytes by allocation
counts would not resolve this measurement problem. Logs are
`/tmp/hybrid-budget-composition-public-v5/allocation_self_compare{112,1132}.log`.

**Hypothesis (untested):** grids from prior fixtures become collectible during
later samples; plan finalizers or their compilation add unrelated work to the
sample. Test this with the same owner and unchanged assertions, recording
compilation and GC statistics, then retaining all fixture grids through the
measurements. The grid-retention control must still execute GC pauses and must
not disable collection or exclude failing samples. If it succeeds, use that
fixture ownership in subsequent paired gates while retaining all numerical,
allocation and latency assertions. No production allocator change is proposed.

The retention control passes all 7,200 assertions on each runtime, with identical
bytes and all allocation counters in every sample, active GC pauses and zero
measured compilation time. Instrumented unretained controls again fail and
record compilation during the mismatching GC samples. This verifies the fixture
retention result; it does not by itself identify every allocating finalizer
frame. Commands, raw statistics and results are retained at
`/tmp/hybrid-budget-allocation-accounting-20261019`. Subsequent paired gates
retain their fixture grids and every existing assertion.

V6 passes 2,204 direct wide-field/sign/stored-BigFloat assertions and 344
endpoint/custom-closure assertions on each runtime. Its signed 80³ case now
retains the original result with 80/96 B, so the zero-allocation gate still
fails. Full allocation profiles locate Float32 and interval-tuple boxes at the
total's range-predicate calls; the direct wider-field helper itself allocates
zero bytes. Evidence: `/tmp/hybrid-budget-composition-public-v6/profile_wide_allocations{112,1132}.log`.
The next candidate checks all optional intervals before each predicate call,
so that call receives concrete intervals. Keep the same shared predicate and
arithmetic, preserve every zero-byte assertion, and verify the allocation
change rather than assuming the narrowing removes boxing.

V5 also retains six failed zero-bound assertions per runtime: its old direct
driver expected `nothing` above the count guard after the zero enclosure had
been moved ahead of that guard. V6's migrated driver checks the explicit
`[-qmin,qmin]` bound and both exact half-subnormal endpoints for every zero/count
combination. All nonzero count, sign, range and large-particle checks remain;
both migrated runs pass. The production source is unchanged by this test
migration; original failures and the migration note are retained.

Public V7 checks optional intervals before the shared range predicate. Its
signed 80³ witness passes all seven assertions on each runtime with zero bytes
allocated. Full allocation profiles also contain no samples. The first profile
driver failed while summing that empty vector; the corrected driver supplies
an additive identity and explicitly checks the empty profile and both zero-byte
measurements. Original logs and `profile-zero-correction.md` remain available.
Both full compiler comparisons match baseline: 22/21 reports and no added or
removed groups. Wider-field checks pass 2,204 assertions and endpoint/custom
contracts pass 344 per runtime. Evidence is in
`/tmp/hybrid-budget-composition-public-v7`; final numerical, paired cost,
permanent-test integration and exact-SHA GitHub gates remain required.

### Permanent composition test integration

The V7 source remains frozen during numerical and performance gates. A separate
integration snapshot will add the accepted offset, total, endpoint, wider-field
and scalable-bound cases to `test/budget_composition`, loaded by `runtests.jl`.
Use the existing independent logarithm oracle; retain integer square-root and
rational boundary oracles. Custom closure expected values come from the known
zero-particle, unit-field fixture, without installing a copied budget method.
Keep the four-thread native driver and permit the permanent suite to exercise
its available CI threads. Preserve ownership and per-call spectral workspaces.

Permanent resource checks compare ordinary component composition and current
scratch allocation, including the signed large-grid zero-byte witness. The
pinned previous-revision, alternating-order latency suites remain external
acceptance evidence; machine timing is not added as a CI correctness assertion.
Run the integrated suite, fixture-loaded compiler analysis, Aqua, migration
checks and final formatting on both Julia versions. Test integration must not
change production executable syntax or relax the frozen candidate's gates.

V6's complete budget regression lane also fails 54 allocation equalities per
runtime (2,625 passes before the failing testset stops the file). V7 passes the
same lane in full: 4,833 assertions per runtime, including those allocations.
V7 also passes 3,968 offset, 328 total, 5,580 initial certificate and 769
scalable-bound assertions per runtime. Full table inventories and log hashes
are retained in `checked-completed-log-tables-20261019.json` in each snapshot.
The older failures remain failures; the correction is accepted only on the
new frozen source.

The separate V8 integration snapshot adds ten composition test files and a
suite loader. Its 358 formatting assertions preserve executable syntax across
all 357 Julia files and pass the complete formatter scan. Migration passes
89 production AST comparisons against frozen V7 and 20 suite-wiring checks.
The 21 changed source/test paths are recorded in `changed-paths.json` under
`/tmp/hybrid-budget-composition-public-v8`. Its numerical, loaded compiler,
Aqua, shared-field, solver, shock and checkpoint gates are queued after V7's
performance work. No integrated result or performance acceptance is assumed.

V8's integrated custom-closure and wider-field smoke passes 2,548 assertions
on each runtime. Its expected custom values are independent of a copied
production budget, and all 280 type/dispatch/ownership assertions are retained.
The controlled cost queue also includes 16 signed 80³/128³ profiles and 72
exceptional closure/current profiles through 4,096 cells. Ordinary signed
profiles require unchanged values/bytes, finite correct signs and the existing
1.5 paired latency target; they do not impose correctly rounded totals where
the public contract deliberately retains ordinary arithmetic. Exceptional
profiles require independently constructed complete zero/unit totals and record
allocation ranges and median elapsed times. Peak memory remains unmeasured.

The cost queue additionally covers 88 public-budget profiles with 1.1/2.1 million
particles, both weight signs and all closure/current combinations allowed by
the model. These retain the value, byte and 1.5 latency assertions, supplementing
the direct scalable-bound tests. V8's private environment audit found a stale
V7 project source path while its manifest already loaded V8. Both project paths
are corrected before full integration; parent and copied-environment child
probes pass on both Julia versions. V5–V7 project/manifest paths already agree.
Evidence: V8 `environment-source-audit.json` and `environment_probe-results.json`.


### October 20 completed gates and remaining cost failures

Main `ffdb522` passes all six exact-SHA GitHub jobs, with 502,711/502,712
assertions on actual Julia 1.12.7/1.13.2 and unchanged complete compiler trees
(23/22). The isolated V8 permanent/Aqua suite passes 13,636 assertions per
runtime; its loaded compiler trees also match baseline (22/21). Each runtime
passes 770 four-thread, 12,689 shared-field, 3,681 solver, 453 shock and 38,055
checkpoint assertions. All seven commands per runtime exit zero in
`/tmp/hybrid-budget-composition-public-v8/integrated-results.json`.

**VERIFIED performance rejection:** V7's 60-profile small suites both pass;
its 324-profile suite fails one minimum-Julia latency assertion (1.535892×,
Float32 8³, unit density, polytropic gamma 1.5, zero inertia). Latest Julia
passes all 324. Both 132-profile large-grid suites pass all 1,716 assertions.
All ordinary value and byte equalities pass. The signed/large-particle/cold
suite passes 4,336 assertions on minimum Julia, but latest Julia stops after
the signed testset with 975 passes and one latency failure (1.514708×,
Float32 80³, temperature 0.75, zero inertia); its particle and cold sets remain
unexecuted. No gate is relaxed or failure rerolled. Complete profiles and log
hashes are in V7 `controlled-cost-summary.json`. V8 is not promoted because
its production ASTs match this rejected source. Profiles of both failing
shapes are being collected before the next canonical optimization.

### Bound reuse and wider-reduction pre-edit brief

- Mode: optimize the two measured costs without changing the range contract.
- Evidence: `/tmp/hybrid-budget-composition-cost-profile-20261020` records both
  native profiles. The polytropic budget interval contributes 311/359 samples;
  the wider field refinement contributes 518/1,633 samples, respectively.
  Source inspection shows the budget repeats the component's finite-density,
  count, maximum-density and exponent-bound work.
- Canonical change: retain the existing polytropic probe's exponent bounds in
  the same internal electron-energy state used by isothermal composition.
  Share the finite-input eligibility predicate with the public component.
  The public scalar still returns the state's original value; the budget uses
  its bounds without a second density traversal. Correctly rounded cold and
  gamma-one values use adjacent floating bounds. Remove the unused scalar
  private probe instead of retaining a compatibility wrapper.
- Wider refinement: expose the inner array reduction to SIMD in the shared
  field-square owner. Only the Float64 certificate sum may reassociate;
  preserve the exact rational loop and stored-BigFloat conversion. Finite
  Float32 squares fit exactly in Float64, and the existing conservative bound
  still covers at most three arrays' additions in any reduction tree. Public
  ordinary field arithmetic remains unchanged.
- Compatibility: retain validation order, unsupported/nonfinite arithmetic,
  custom closure dispatch, zero/sign behavior, computed current and ownership.
  No tolerance, report suppression or performance threshold changes.
- Verification: independent interval enclosures and all permanent/native gates;
  exact component equality; full loaded compiler comparisons; all original
  ordinary/signed/particle/cold profiles. Pin the polytropic baseline probe and
  owner to `ffdb522` as well, so metadata work cannot change both benchmark sides.
  Measure before claiming either optimization succeeds; retain failed V7 logs.


V9's original state-bound driver fails 12 assertions per runtime because its
768-bit independent logarithm remainder cannot resolve the Float64 zero bin
for an exact log identity. The corrected driver uses the same independent
series with a 1,280-bit remainder for all gamma-one fixtures; no production
change or assertion removal accompanies that correction. Both corrected
runs pass 3,809 assertions. The two profiled cost cases then pass 100 checks
per runtime, with exact original values and zero allocation. Paired ratios
are 0.9383/1.0457 on Julia 1.12.7 and 1.0339/1.0323 on 1.13.2 (small polytropic /
large signed). These measurements do not replace the complete cost matrices.
Evidence: `/tmp/hybrid-budget-composition-public-v9`, including the original
failed logs, `oracle-precision-correction.md` and `targeted_costs-results.json`.

V10 integrates the 492 new state-bound assertions into the permanent suite.
Formatting passes 359 checks across all 358 Julia files. Migration passes
89 unchanged production ASTs, 22 suite-wiring checks and two exact pinned
polytropic-owner comparisons. Its full acceptance sequence ran in
`/tmp/hybrid-budget-composition-public-v10/run_acceptance.py`: core numerical,
loaded compiler and four-thread gates on both versions, then serial controlled
cost matrices, then solver/shock/checkpoint integration on both versions.
All gates remain required before promotion. Main remains `ffdb522`.

V10's permanent/Aqua suites pass 14,128 assertions per runtime; both complete
fixture-loaded compiler trees match baseline (22/21 reports, no added or
removed groups). Both four-thread suites pass 770 assertions. Project and
manifest graphs still match V9 after substituting only the intended source
path (`environment-graph-after-permanent.json`). Both shared-field suites pass
12,689 assertions. **VERIFIED rejection:** both component-regression commands
then exit 1 after 10,606 passing assertions: the constant alias at
`test/polytropic_energy/ordinary.jl:5` still names the removed private scalar
probe. The earlier call-pattern search missed the alias. Cost and broad gates
did not start; `core-results.json` and both failed logs are retained.

V11 migrates that test directly to the state-returning probe, retaining every
existing numerical assertion and both four-ULP limits. It adds independent
rational containment assertions for the returned bounds and searches all
source/test files for the old symbol, including aliases. Production arithmetic
is unchanged from V10. README documents why the complete total can differ from
the sum of separately rounded components. Formatting passes 359 checks.
`/tmp/hybrid-budget-composition-public-v11/run_acceptance.py` runs component
regressions first, then the permanent/compiler/threaded/field gates, serial
controlled costs and solver/shock/checkpoint integration on both Julia versions.
No candidate has been promoted.

V11's component regressions now pass 49,321 assertions on each runtime,
including 28,606 polytropic checks. Migration passes 384 assertions: production
AST equality, suite wiring, pinned benchmark syntax, removal of every old
probe reference and preservation of the original test limits. All 358 frozen
Julia files and 16 external Julia drivers retain their recorded hashes.
Permanent/Aqua (14,128), four-thread (770) and shared-field (12,689) checks
also pass on each runtime. Complete fixture-loaded compiler trees match
baseline at 22/21 reports with no added or removed groups.
Evidence: `regression-results.json`, `preflight-results.json` and
`checked-completed-log-tables-20261019.json` in the V11 directory.

**VERIFIED V11 performance rejection:** the first 60-profile cost matrix
retains all public values and allocation counts, but fails one minimum-Julia
and two latest-Julia latency assertions. The failed near-unit-density cases
are Float32 1D/eight-cell evolving pressure (1.57204×, minimum Julia), Float32
1D/four-cell CGL (2.21840×, latest) and Float64 1D/four-cell evolving pressure
(1.56877×, latest). The full ordinary, large-grid, signed/particle/cold and
broad integration stages did not run. Original logs and
`controlled-cost-summary.json` remain unchanged. V11 is not promoted.
CPU profiles of these exact fixtures and their pinned baseline are being
collected in `/tmp/hybrid-budget-small-cost-profile-20261020`. No timing cause
is inferred before inspecting those profiles; no failed gate is rerolled.

### Shared closure validation — pre-edit brief

- Mode: remove repeated validation in the public budget's CGL and evolving
  pressure paths while retaining the complete-total contract.
- Verified source path: `energy_budget` validates B/n before evaluating
  components. Its generic `_budget_electron_state` then calls the public
  electron diagnostic, repeating the same B/n grid checks. The minimum-Julia
  CPU profiles include 1,124 samples in that CGL electron-state path and
  712/856 in the two evolving-pressure paths. These counts identify inspected
  work; they do not establish the cause of every observed timing ratio.
- Canonical change: extract the public evolving-pressure argument, storage,
  gamma and energy block into one internal owner shared by the public
  component and budget. The budget supplies already-validated B/n, while that
  owner retains every pe/gamma check. For supported CGL storage, the budget
  checks the algebraic pe argument and calls the existing checked CGL owner.
  Other precisions and custom closures keep their existing dispatch.
- Compatibility: preserve validation ordering and exception text; no new
  closure formulas, result tolerances, unchecked pressure input or public API.
  Keep scalar component values, range completion and exact coefficient rules.
  Benchmark baseline owners remain pinned to the pre-composition source.
- **Hypothesis (untested):** avoiding the repeated public wrappers will bring
  the rejected small profiles within the existing 1.5× gate. Settle this with
  the same full 60-profile matrix on both runtimes, followed by every remaining
  numerical, compiler, concurrency, cost and integration gate. Test the small
  matrix first so a failed optimization cannot be mistaken for acceptance.
- Migration: prove that all production syntax outside these declared blocks
  is unchanged; retain invalid-shape, missing/alternate pressure, wrong storage,
  invalid gamma, nonfinite and custom-closure tests. Keep all V11 failed logs.

V12 is prepared at `/tmp/hybrid-budget-composition-public-v12`. Formatting
passes 359 checks; the initial migration passes 387, including exact equality
of the original evolving-pressure block and new helper body. Normalization
removes only the three declared additions and restores only the factored
conditional before comparing all 89 production ASTs. Its 24 changed paths
remain isolated. Both profiling commands completed successfully; the full
acceptance controller tested the small cost matrix before broader gates.

**VERIFIED V12 rejection:** all value and byte comparisons pass, but the
60-profile small matrix fails four minimum-Julia and five latest-Julia
latency checks. Ratios span 0.6935–2.5158 and 0.3359–3.2461, respectively.
Full profile identities, commands and log hashes are retained in V12
`controlled-cost-summary.json`. Component, permanent/compiler, threaded,
field, remaining cost and broad integration stages did not run.

The identical-function timing control is complete at
`/tmp/hybrid-budget-timing-accounting-20261020`. Both measured sides call the
same pinned function object (asserted explicitly), retaining the original
fixtures, pairing, sample counts, value/byte comparisons and 1.5× threshold.
Each timed sample records compilation, recompilation and GC time. This tests
whether the current measurement can reject identical work; it cannot accept
V12 or erase either rejected matrix.

**VERIFIED minimum-Julia control:** the identical-function matrix finishes
with 1,614 passes and eight timing failures. Its 60 ratios span
0.3158–6.3817 even though both sides call the same function object.
All 720 sample pairs record zero compilation, recompilation and GC time.
`control-summary.json` retains every profile/sample and the log hash.
This proves false regression classifications in this procedure; it does not
identify the external cause of the elapsed-time variation or accept V12.
Latest Julia also fails: 1,620 passes and two timing failures, with ratios
0.5810–2.5838. All 720 of its sample pairs likewise have zero compilation,
recompilation and GC time.

### Timing calibration repair — pre-edit brief

- Keep the original fixtures, exact values, allocations and 1.5× elapsed-time
  ratio limit. Preserve all rejected measurements.
- Calibrate using warmed batches instead of one scalar call. Increase the
  count until five measured batches each take at least 20 ms, using checked
  integer growth. Remove the arbitrary 4,096-iteration ceiling that can leave
  fast cases with short samples.
- Use 48 alternating pairs and retain every sample. Keep the original ratio
  of medians gate; add the median of within-pair ratios at the same limit.
  Record compilation, GC and elapsed data outside the timed calls.
- Validate the procedure first on all 60 identical-function profiles, including
  a reciprocal bound so it cannot favor either label. Deliberate doubled-work
  controls must exceed the unchanged regression limit. If these controls fail,
  the procedure remains unaccepted; do not promote a candidate on its results.
- After control validation, apply the same procedure to the existing full
  matrices and rerun all required numerical/compiler/integration gates on a
  new frozen acceptance snapshot. This is a measurement repair supported by
  identical-function failures, not a reroll of the original gate.

The first calibrated prototype is running at
`/tmp/hybrid-budget-timing-calibration-v1-20261020`. It is already rejected:
its deliberate doubled-work wrapper allocates 112 bytes versus the baseline
wrapper's 16, and one identical-call profile reaches a ratio of 1.54928.
Both assertions remain enforced and their logs are retained. A subsequent
prototype must correct the wrapper, first pass a short numerical/allocation
preflight, and demonstrate timing reliability and slowdown sensitivity.
No replacement protocol or composition implementation is accepted.

Both first-prototype controls are complete. Minimum/latest Julia record
192 slowdown-wrapper byte mismatches each, plus three/four identical-call
timing assertion failures. The complete profiles and assertion classifications
are in `completed-summary.json`. The second prototype at
`/tmp/hybrid-budget-timing-calibration-v2-20261020` specializes the sampler
arguments and uses direct, observable doubled batch results. It increases the
minimum warmed batch to 50 ms while retaining all 48 pairs, both 1.5× limits
and reciprocal identical-call controls. Each runtime passes all 16 short
preflight assertions, including zero allocations in both wrappers. Full
identical-owner and slowdown-sensitivity controls are running; their results
remain required before any candidate timing matrix.

The second prototype is also rejected before completion. On minimum Julia,
Float32/N=8/general/polytropic-gamma=1.5 identical calls yield ratios
0.640656 and 0.609563, below the reciprocal limit 2/3. Both failures remain in
`control112.log`; all runtimes continue to completion for the full record.
The prepared V13 creator requires successful complete controls and therefore
cannot create an acceptance snapshot from this failed protocol.

Minimum Julia has now completed: 9,342 passes and two failures, both on that
profile. All four deliberate-slowdown controls pass, and all 3,072 recorded
sample pairs show zero GC, compilation and recompilation time. Identical-call
ratios span 0.640656–1.446462; paired medians span 0.609563–1.228870.
`minimum-completed-summary.json` retains every sample and the completed log
hash. Latest Julia subsequently completes all 9,344 assertions successfully:
identical-call ratios span 0.762858–1.332155 and paired medians
0.794188–1.182381. All 3,072 sample pairs also have zero GC/compilation time.
Both runtimes' four deliberate-slowdown profiles exceed 1.5. The complete
`completed-summary.json` recomputes every reported ratio from the samples and
records `accepted=false`; minimum Julia's failures prevent protocol acceptance.

A separate diagnostic is complete at
`/tmp/hybrid-budget-thread-timing-20261020`. It records elapsed and
calling-thread execution times for the three most
deviant identical-call profiles on each runtime, preserving their fixtures,
batch counts and all sample metadata. It uses the local SDK clock API and
checks thread identity and a sleeping control. Its ratios cannot accept a
candidate or replace the elapsed-time gate. The separate Float32 accuracy
prototype runs afterward.

Each runtime passes 436 checks, including the sleeping/monotonicity control
and 144 identical-call pairs. `summary.json` retains all samples and hashes.
Across the six profiles, elapsed ratios are 0.792513–1.369321 and
thread-execution ratios are 0.982471–1.015660. Median execution fractions of
the enclosing wall interval are 16.1–21.2%; median thread execution times are
14.4–24.8 ms for batch counts selected by the earlier 50 ms wall-time
calibration. All samples have zero GC, compilation and recompilation time.

These measurements establish that the selected batches did not ensure
50 ms of thread execution. They do not identify a competing process or
operating-system policy, establish a hardware limit, or accept the rejected
composition candidate. A subsequent protocol can use measured thread execution
to calibrate batch work while retaining the elapsed-time ratio gates; it still
requires complete identical-owner and slowdown controls before use.

The separate accuracy assessment against exact `ffdb522` is complete. Both
native versions reproduce large-grid Float32 magnetic-energy errors, with
exact rational integer-multiplicity controls. This is tracked separately as
[ENERGY-010](MAGNETIC_ENERGY_ACCURACY_PLAN.md); it does not change the ordinary
component compatibility contract of this composition repair.


### Next timing-protocol experiment

Prepare a new isolated protocol that changes only batch calibration: require
five warmed samples to each contain at least 50 ms of measured calling-thread
execution, using the same SDK clock shim already checked above. Retain the
48 alternating pairs, elapsed ratio-of-medians and paired-median 1.5 limits,
reciprocal identical-call bounds, all 60 fixtures and four deliberate-slowdown
controls. Record both CPU and elapsed calibration values; keep the original
paired sampler unchanged. Calibration CPU time must not replace any acceptance
measurement. Run preflight and the complete controls on both runtimes after
particle-growth verification, without overlapping our other native work.

This is a new, unaccepted measurement procedure. The earlier failed controls
and rejected composition candidates remain rejected. Only complete controls
can authorize a subsequent composition acceptance snapshot.
