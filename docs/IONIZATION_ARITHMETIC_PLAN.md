# Ionization energy and optical-depth arithmetic

Status: COLL-001 confirmed on Julia 1.12.7 and 1.13.2; an isolated repair is
under verification and has not been promoted. This is separate from the GROWTH-001 storage preflight
and GROWTH-002 migration repairs.

## Verified baseline

The frozen reference source at
`/tmp/hybrid-rank-allocation-v2-20261020/source` matches pushed `b3e2eb6`
across all 92 source/extension files. `src/particles/collisions.jl:1029–1042`
specifies strict eligibility `m|v|²/2 > E_iz`, probability
`1-exp(-nσ_iz|v|dt)`, and primary cooling by `E_iz`.

**VERIFIED:** `/tmp/hybrid-ionization-energy-range-review-20261020/probe.jl`
passes 132 assertions and fails 30 on each runtime, with no fixture errors.
The matrix covers Float32/64 and 1D–3D; every velocity component, mass and
energy input is finite. Exact rational energy ratios and a serial 512-bit
square-root oracle establish eligibility and finite cooled components.
Ordinary controls pass. When both represented speed and threshold are
infinite, six eligible cases return zero births. In six further cases, the
birth occurs but the primary retains every original velocity component.
Complete records, failure-line counts and all log/source hashes are retained
in `energy-summary.json`.

**VERIFIED cause:** the owner first collapses the norm to `T` at line 1101
and compares it with the collapsed threshold from line 1097. The recorded
`Inf > Inf` comparison skips eligible particles. At lines 1168–1176, a
finite threshold divided by an infinite norm produces a zero loss ratio and
unit cooling factor. The recorded unchanged primary values reproduce this
path.

**VERIFIED:** the independent `probability.jl` probe passes 144 assertions
and fails 18 per runtime. With three `floatmax(T)` components and rate
`floatmin(T)/128`, exact arithmetic bounds the finite optical depth below
`1/16`. A uniform sample of `1/2` must reject the event, but all six extreme
cases create a pair and advance both counters. Ordinary rare-event and
saturated controls pass. `probability-summary.json` retains all records and
hashes. The source passes the already infinite norm to
`_collision_probability` at line 1103; its finite-factor reconstruction at
lines 67–90 cannot recover the discarded speed.

The first additional threshold probe used an ambiguous parametric RNG
overload. Its failures are retained and are not accepted evidence. The revised
probe at `/tmp/hybrid-ionization-boundary-review-v2-20261021` uses explicit
Float32/64 methods with native RNG controls, preserving all behavior assertions.
**VERIFIED:** both runtimes retain 77 passes and 19 behavior failures, with
every error result `nothing`; `summary.json` checks all hashes, failure lines
and twelve records per runtime. Float32 creates a birth at exact equality
`v=(1,0,0), m=1, E_iz=0.5`; Float64 skips a birth for `nextfloat(1)`. Tiny
orthogonal components expose additional lost eligibility/cooling, including
cases where a completed primary component is representable while a standalone
cooling scale would underflow. Below-threshold controls pass.

## Repair scope and acceptance

Repair the existing ionization owner, preserving the public signature,
validation order, newborn sampling order, ID semantics and storage guards.
Retain complete finite products and compare the strict energy threshold
before narrowing range. Cooling must preserve representable completed
components even when an intermediate speed or scale is unrepresentable.
Keep ordinary cases allocation-free at the scalar arithmetic layer and
measure full-call cost against the accepted parent.

The neutral-MCC caller checks that its relative speed is finite before using
the shared probability helper (`collisions.jl:976–980`). Do not silently
change its documented range-error behavior as part of this repair.

Before production edits, complete a source-backed design brief. Check existing exact-arithmetic helpers before
adding a primitive. Avoid BigFloat transcendental calls in concurrent
production paths while NUM-001 remains open.

Required gates:

- Preserve all original failed probes and verify the repaired public call.
- Cover equality, immediately adjacent energies, tiny orthogonal components,
  subnormal/overflow intermediates, sign/permutation symmetry and empty/no-event
  cases. Verify complete populations and counters.
- Check error paths, selected-batch validation and storage rejection together
  with the revised GROWTH-001 guard.
- Verify any new arithmetic primitive independently with exact midpoint and
  range boundaries, then threaded controls on both Julia versions.
- Run collision, particle-growth and integration regressions; concrete and
  complete package compiler comparisons; formatting; full-output and resource
  comparisons. Record any behavior or allocation change explicitly.
- Promote this repair separately after its parent and require all six GitHub
  jobs green for the exact commit.

## Design evidence and pre-edit brief

The isolated design at `/tmp/hybrid-ionization-arithmetic-design-v1-20261021`
passes 35,094 exact square-root rounding-cell checks and 2,400 threaded checks
per runtime. The primitive compares integer squares with exact squared
midpoints, including subnormal and overflow endpoints; every finite positive
Float16 square is covered. It does not call MPFR transcendental functions.

The scalar design also passes 10,296 boundary/reference checks and 144 ordinary
compatibility/resource checks per runtime. All 36 ordinary cases retain the
old probability and cooled-component bits with zero scalar allocations.
These are prototype checks, not public-call or full-step acceptance.

- Mode: repair the established ionization contract.
- Evidence: the three verified public probes above; the current speed/threshold
  comparison and unit cooling factor lose information before event selection
  or mutation. The shared exact owner is `src/utils/exact_arithmetic.jl`,
  included before `particles/collisions.jl` by the module.
- Contract: strict exact eligibility for stored finite components, probabilities
  from the complete speed/rate/time expression, and finite completed cooling
  components. Preserve validation, RNG draw order for eligible events, storage
  guards, ID reservation, newborn assembly and public signatures. Borderline
  events that violated the stated strict threshold intentionally change.
- Strategy: enclose the non-negative energy expression with outward adjacent
  floats. Decide eligibility only when those bounds prove it; otherwise compare
  exact stored values. Retain ordinary cooling only when an outward component
  enclosure certifies at most six machine eps relative error. Compute uncertain
  completed component squares exactly and round their square roots directly.
  Reconstruct complete optical depth only when the represented speed overflows;
  the existing finite-speed probability owner remains shared with neutral MCC.
- Ownership: private ionization helpers in `collisions.jl`; one exact
  square-root primitive beside existing exact conversions. No new public API,
  alternate collision implementation, catch, diagnostic suppression or flag.
- Compatibility: require ordinary bit/allocation comparisons and native
  Float16/32/64 plus BigFloat controls before accepting the new primitive.
  Full-call allocation and timing remain to be measured. No scalar benchmark
  establishes end-to-end performance.
- Verification: original public probes, exact midpoint/range cells, neighboring
  and sign/permutation cases, threaded arithmetic, validation/storage failures,
  existing collision/growth suites, compiler and format gates, then combined
  integration and exact-commit GitHub acceptance.

An infinity-only patch misses the finite threshold failures; narrowing a
cooling scale before multiplying loses the wide-exponent completed component;
saturating every overflowed norm fails the rare-event probability witness.
Using exact arithmetic for every ordinary input would fail the zero-allocation
gate. These cases are retained as independent acceptance tests.

## Isolated implementation

V1 at `/tmp/hybrid-ionization-arithmetic-v1-20261021` changes only the collision
owner, exact-arithmetic utility, two new permanent test files and their runner
includes relative to the revised growth/rank parent. It preserves the backing
storage guard. Formatting passes 351 assertions; both four-thread native runs
pass 38,043 assertions, including every original public witness, exact-root
cells, Float16 calls, BigFloat caller-context checks and stored BigFloat inputs
whose precision exceeds the working precision. All command/log and frozen
source hashes were checked.

V2 at `/tmp/hybrid-ionization-arithmetic-v2-20261021` retains all 92 production
source/extension files byte-for-byte. It makes the additional type controls
permanent and expands public coverage with complete population/counter/RNG
checks, coordinate permutations, exponent samples and concurrent near-threshold
births. Both runtimes pass 53,519 focused assertions, 5,879 neighboring
integration assertions and 34 concrete compiler checks. All nine native
command/log hashes and 360 frozen files were checked. Four complete
unloaded/fixture-loaded compiler comparisons retain the parent's 22/21
reports without changed groups; all eight capture/comparison logs and 710
frozen files were checked. Four resource commands pass 720 assertions each:
all 120 complete outputs per runtime match, with zero allocation change in
all 24 ordinary profiles and zero measured compilation/recompilation.
All resource logs and 195 frozen files were checked. These remain reference
results because the endpoint probe below rejects V2. Unpaired elapsed-time
observations do not establish speed or timing acceptance.

An additional reduced-exponent BigFloat probe checks the new primitive's
rounding at exponent endpoints under several caller rounding modes.
**VERIFIED:** `/tmp/hybrid-ionization-bigfloat-range-review-20261021`
retains 78 passes and 12 failed nearest-root assertions on each runtime.
All exact-boundary and caller-context controls pass; there are no fixture
errors, and both processes restore their exponent ranges. `summary.json`
checks every command/log/source hash and all 30 records per runtime.
V1 and V2 now have explicit rejection markers and must not be promoted.

**VERIFIED cause:** the initial primitive rounds a significand before final scaling. That loses
which side of the half-minimum boundary the exact root occupies. In addition,
Julia's installed `base/mpfr.jl:841–844` passes the caller's rounding mode to
the final `ldexp`, producing finite maxima where this primitive specifies
nearest overflow to infinity. The probe reproduces both paths.

Follow-up brief: repair `_sqrt_from_exact` itself. Compare roots below the
minimum positive value with the exact squared half-minimum boundary; decide
overflow explicitly after integer rounding. Keep ordinary finite scaling
exact, with no caller-setting changes or new MPFR transcendental calls.
Preserve the failed probes and add the endpoint test permanently in an
isolated child process, so its reduced exponent range cannot affect other
tests. Require the original rounding cells, type/context controls, public
matrix and full remaining gates on the revised source.

V3 makes that endpoint decision explicit and passes 53,521 focused assertions
per runtime, including the original 90 endpoint assertions in an isolated
child. A further compact-range probe at
`/tmp/hybrid-ionization-bigfloat-scale-review-20261021` then retains 42 passes
and six failures per runtime. At 96/128-bit precision with maximum exponent
63, `_sqrt_from_exact(BigFloat, 4//1)` returns an incorrect value despite the
result `2` being representable. All caller-context and range controls pass,
and both child processes restore their exponent ranges. The failed assertion
locations, twelve records per runtime and hashes are checked in `summary.json`.
V3 is also rejected.

**VERIFIED cause:** V3's `T(whole)` conversion precedes scaling in
`src/utils/exact_arithmetic.jl`; the rounded integer has more exponent bits
than this compact BigFloat range permits. Follow-up: retain the integer
rounding policy and construct the final BigFloat with `mpfr_set_z_2exp`,
which accepts the integer and binary exponent together. The installed MPFR
header declares this operation; Julia's `base/mpfr.jl` supplies the matching
`BigFloat`/`BigInt` reference layout and rounding enum. Use explicit nearest
rounding and the caller's precision without changing caller settings.
Extend the permanent child fixture to cover the wider precisions and normal
values, then rerun all primitive/public and remaining acceptance gates.

V4 at `/tmp/hybrid-ionization-arithmetic-v4-20261021` changes only the exact
square-root helper and its child fixture relative to V3. It constructs the
rounded BigFloat integer and exponent together, and compares the squared
half-minimum boundary using integers. The expanded fixture extracts stored
dyadic values with Julia's `Base.decompose`, avoiding floating scaling in its
independent endpoint oracle. It covers 24/53/96/128-bit precision, three caller
rounding modes, normal values and both exponent endpoints.

**VERIFIED:** formatting passes 352 assertions; both four-thread focused
runs pass 53,521 assertions, with 252 additional assertions in each isolated
endpoint child. All three command/log hashes and 360 frozen files were checked.
The unchanged original endpoint and integer-intermediate probes also pass,
with 90 and 48 assertions per runtime, changing only their asserted package
path. Their complete logs, frozen sources and restored exponent ranges are
checked in `/tmp/hybrid-ionization-bigfloat-repair-v4-20261021/verified-results.json`.

The separate exact-fallback resource matrix passes 1,576 assertions per
runtime. It covers Float32/64, near-threshold, wide-component, overflowing-norm
and rare-event cases at 1/16/128 incident particles, with five samples per
profile. Complete populations, IDs, counters, RNG draws and independent
8192-bit cooling references are checked. All 120 samples per runtime have
zero measured compilation/recompilation and stable full-output hashes within
each profile. Observed full-call allocations range from 8,336 to 3,897,584
bytes; exact fallbacks allocate. Both logs and 362 frozen files were checked
in `/tmp/hybrid-ionization-cold-resources-v4-20261021/verified-results.json`.
These cost observations do not establish elapsed-time acceptance.

Both neighboring integration suites now pass 5,879 assertions. The four
ordinary-resource commands pass 720 assertions each, with all 120 full
outputs per runtime matching the parent and zero allocation increase across
all 24 profiles. All measured calls have zero compilation/recompilation;
the complete resource logs and 195 frozen files were checked.

Concrete compiler checks and full compiler comparisons remain in progress.
No V4 acceptance, performance improvement or production promotion is claimed
yet.
