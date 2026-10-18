# Polytropic energy range and closure precision

Status: pointwise pressure is pushed as `0d8c8e8`; its latest-Julia CI job
failed resource-fixture byte-equality assertions. The corrected fixture,
EMPIC initialization fixture and algebraic dispatch/compiler repair are pushed
as `15662e1`, with [exact-revision CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39137529999)
running. Both native solver and shock lanes pass; checkpoint lanes continue.
The complete polytropic energy repair is integrated from verified V5. Its
core gates pass on both runtimes; native broad and exact-SHA CI acceptance
remain pending. This plan covers ENERGY-004, ENERGY-007 and CLOSURE-001/002.

## Confirmed behavior and implementation contract

**VERIFIED — ENERGY-007:** a Float64 gamma immediately above or below one
becomes exactly one when converted to Float32. The diagnostic then takes its
pressure-weighted logarithmic branch. At unit density/reference density and
`pe0 = 0.5`, it returns zero instead of the nonzero stored-closure integral.
The independent oracle is rational because `(n/n0)^gamma = 1` exactly.

The six 1D–3D witnesses have finite nonzero reference values from approximately
`4.50e15` to `3.60e16` in magnitude. Both runtimes pass 40 witness assertions,
including Float64-grid and exact-gamma controls. The witness diagnostic body
is byte-identical to `63f2420`: narrowing and branch selection there are at
`src/integrators/extrapolated_leapfrog.jl:1402` and `:1405`; the pressure law
retains stored gamma at `src/electrons/closures.jl:85`. Commands and hashes:
`/tmp/hybrid-polytropic-complete-design-v1/gamma-evidence.json`.

The earlier [ENERGY-004 witnesses](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#algebraic-closure-range--remaining-energy-004-scope)
also establish power overflow/underflow and density-ratio overflow for gamma
two and one-half. The repair must retain stored gamma, density/reference ratio,
coefficient, denominator and each spacing through the complete expression.
Exact gamma one uses the coefficient `pe0/n0` in the logarithmic free energy;
rounding that ratio first can also lose a finite complete integral.

At unit density the raw polytropic reservoir is proportional to `1/(gamma-1)`
while isothermal free energy is zero. The isothermal documentation overlay
removes the false claim that the latter is the exact limit of the unshifted
former. It preserves the polytropic implementation while its separate repair
is developed.

## Pointwise pressure — CLOSURE-001

**VERIFIED:** 44 assertions pass per native runtime against the pressure owner
byte-identical to `63f2420`. Six complete finite pressures are replaced by Inf
or zero after intermediate power/ratio range loss. Two normal-output cases use
`n = 2`, `n0 = prevfloat(2)` and `gamma = 2^precision(T)`: the routine returns
approximately 7.389 while the complete-input power enclosure is near 2.718.
The range witnesses have exact rational or square-root controls. Evidence:
`/tmp/hybrid-polytropic-pressure-witness-20261018/evidence.json`.

A separate eight-assertion classification per runtime encloses powers of the
already-rounded density ratio. Every observed native pressure correctly rounds
that altered power, isolating ratio rounding as the source of the large normal-
output discrepancy. Evidence: `ratio-classification-evidence.json` in the witness
directory. The six shock expressions listed below have separate public witnesses.
Their permanent regressions pass in public pressure v2; broad solver/restart
and exact-revision GitHub gates remain pending.

The canonical API is `electron_pressure!` in `src/electrons/closures.jl`,
with the checked Float32/64 owner in `src/electrons/polytropic_pressure.jl`. Its result feeds
`gradient!` through `src/electrons/ohms_law.jl:85`–`:87`, as well as the energy
diagnostic. Repairing only the diagnostic would leave this pressure defect.
The implementation must share complete power arithmetic while retaining the
unrounded factors needed for the integrated energy. Preserve broadcast/storage
contracts, density flooring, nonfinite/generic behavior and validation order;
verify ordinary accuracy as well as range, then test forces and solver/restart
behavior and measure whole-grid cost.

The pointwise cold evaluator in `/tmp/hybrid-polytropic-pointwise-v1` passes
508 rounding/range assertions, two concrete compiler checks and 614 parent-checked
threaded assertions per native runtime. Its `rounding-evidence.json` and
`quality-evidence.json` retain source/log evidence. It treats a single rational-power tie by
comparing numerator and denominator powers with the exact boundary. Bit-length
checks bound each constructed integer by the existing target size, including
very large stored gamma. This is distinct from the remaining multi-term
all-rational energy tie problem. No public pressure code has changed.

## Scalar pressure and shock coverage

Scalar v2 in `/tmp/hybrid-polytropic-ordinary-v2` passes 13,850 assertions per
Julia 1.12.7/1.13.2, including independent rational/root bounds, invalid and
nonfinite values, output precision, signed zeros, resources and parent-checked
threads. Ordinary calls allocate zero bytes on the measured grids through
4,194,304 cells. Its four-output-ulp ordinary controls do not establish universal
correct rounding; exceptional range and midpoint cases use complete intervals.
The 12 concrete scalar/compiler checks are clean on each runtime.

The original ownership harness retains 1,552 passes and eight failed old-bit
comparisons for positive ordinary values. A 19-assertion classification checks
all eight against the declared independent bounds. The corrected harness uses
old-bit equality for zero/base-zero behavior and independent bounds for positive
values. Candidate arithmetic is unchanged. Evidence: `native-evidence.json`.

Six shock expressions also require the canonical pressure arithmetic:
`reflecting_wall_shock.jl:251`, `:305`, `:480`, `:481`,
`reflecting_wall_shock2d.jl:206`, and `reflecting_wall_shock3d.jl:308`.
Both corrected native witnesses pass 108 assertions for force preparation,
initialization, reservoir updates and signed densities. Unlike the closure,
these callers retain raw signed powers. The first minimum-version witness has
14 passes and one error: the 2D spectral operator rejects overflowed pressure
instead of returning a nonfinite field. The corrected witness requires that
exception and preserves its exact finite-pressure oracle. Original source/log
and operator hashes remain in `/tmp/hybrid-polytropic-shock-witness-20261018`.

## Public pressure candidates

V1 in `/tmp/hybrid-polytropic-pressure-public-v1` passes 13,603 permanent
assertions on each native runtime. Earlier array/shock preflights pass 448 each;
formatter checks preserve 216 source/test syntax trees, 22 arithmetic helpers
match their components, and the generic owner retains its original expression.
Its 24 measured ordinary calls allocate zero bytes. Exceptional Float64 controls
allocate 301,312–318,976 bytes at 32 cells and 2,410,496–2,551,808 bytes at 256
cells; exceptional Float32 controls allocate zero bytes.

The first minimum-version concrete compiler gate has 100 passes and eight
failures; latest passes all 108. A 56-assertion classification reproduces all
eight reports in the unchanged original pressure owner, in Base's matrix-view
`unaliascopy` path, and checks runtime copying/shared-parent calls. The final
comparison retains every report and compares the complete unchanged copy
subtree: 124 assertions pass on minimum Julia, 108 on latest. The eight existing
minimum-version reports remain visible. Evidence is in
`hybrid-polytropic-pressure-public-v1-verification/classified-focused-results.json`
under `/tmp`; the original failed gate is retained.

V2 in `/tmp/hybrid-polytropic-pressure-public-v2` combines the checked scalar
and array optimizations below. It is the active public acceptance candidate;
V1 remains a checked reference. Both are based on the isothermal
public v3 snapshot, now pushed as `5323c4d`. V2 changes six Julia files relative to V1, with no change
to integrated polytropic energy. Formatter checks preserve 219 source/test syntax
trees; all 24 component helpers match. Its 308-file frozen snapshot passes
27,704 numerical/resource assertions on each Julia 1.12.7/1.13.2. Concrete compiler gates pass 172 assertions on minimum Julia and 156 on latest,
including the 48 new cache calls; the eight classified minimum-version Base
view-copy reports are retained. Commands, source/log hashes and exit status are in
`hybrid-polytropic-pressure-public-v2-verification/focused{112,1132}-evidence.json`
under `/tmp`. Broad acceptance is partitioned into the
same complete 22 files: two checkpoint, eleven solver/diagnostic and nine shock
files, followed by quality, threaded and full compiler-tree gates. No assertion
or file is omitted. The original combined gate artifacts remain retained. The
19-file source/test pressure correction, including the independent binomial
regressions below, is pushed as `0d8c8e8` to run exact-SHA GitHub CI concurrently.
No pending native or GitHub gate is treated as accepted.

A supplementary independent binomial-series oracle passes 1,314 assertions on
each runtime, including exact polynomial/square controls, independently searched
rounding boundaries, mixed precision and actual stored non-integer gamma.
It uses [NIST DLMF 4.6.7](https://dlmf.nist.gov/4.6#E7) and a derived rational
geometric-tail bound, without candidate logarithm/exponential calls or MPFR.
Evidence and the derivation are in `/tmp/hybrid-polytropic-binomial-oracle-v1`.
The permanent regression overlay passes the same 1,314 checks per runtime;
parsed syntax verifies that wrapping and formatting preserve all assertions.
It is separate from the frozen candidate until promotion, and its final suite
will also require exact-SHA GitHub CI.

## Measured pressure optimizations

Scalar v3 (`/tmp/hybrid-polytropic-ordinary-v3`) passes 26,535 assertions per
native runtime, including 12,685 new rational Float32 residual/remainder checks.
All-Float32 calls may use native Float32 powers with a precision-specific guard;
unresolved calls use Float64 before complete evaluation. Mixed/Float64 outputs
retain wider arithmetic. A latest-version paired run passes 120 assertions.
A nine-pair process-CPU/wall run passes 168 more: Float32 CPU ratios versus v2
are 0.885–0.982, and Float64 ratios are 0.969–1.007. Ordinary measured calls remain
allocation-free. These measurements do not establish a universal speedup.

Array reuse v1 (`/tmp/hybrid-polytropic-array-reuse-v1`) passes 1,658 assertions
per runtime: 1,152 scalar/API controls, 264 resource checks and 242 compiler/
threaded checks. A call-local cache reuses consecutive equal finite densities,
distinguishes signed zero and passes each nonfinite input to the scalar owner.
Measured exceptional Float64 uniform/four-block arrays use 9,416–39,872 bytes
at each actual size: 32, 256 and 4,194,304 cells. Float32 controls use zero bytes.
This is a repeated-density result, not an allocation bound for arbitrary data.
Varied ordinary calls retain zero allocations; the 168-assertion latest CPU
comparison measures large-grid ratios of 1.001–1.008 versus direct scalar calls.

The first reuse harness has 39 passes and one error because its raw negative-
infinity reference raises `DomainError`. The corrected harness compares exact
exception type/message as well as returned values; candidate code is unchanged.
Both component directories retain sources, commands, hashes, logs and measurements
in `native-evidence.json`. Public integration must preserve these contracts and
verify all affected callers before promotion.

## Completed component checks

Counts below are per native runtime unless stated otherwise. Component success
does not establish public final rounding, performance or solver acceptance.

| Component | Checked result | Evidence |
| --- | --- | --- |
| Reduced rational exponential | 868 assertions against an independent alternating expansion and exact controls; complete tail bounds and negative-argument reciprocals. | `complete-design-v1/exponential-evidence.json` |
| Scaled exponential | 175 assertions, including binary exponents of magnitude `10^300`; powers remain scaled. | `complete-design-v1/scaled-exponential-evidence.json` |
| Rational-power intervals | 909 assertions, including exact integer/half/quarter/eighth powers and stored Float32/64 gamma boundaries. | `complete-design-v1/power-evidence.json` |
| Scaled positive sums and endpoint conversion | 730 assertions for exact sums, signs, rounding and enormous exponents. | `complete-design-v1/scaled-sum-evidence.json` |
| Complete polytropic interval | 111 assertions for stored gamma, ratios, coefficients, spacings, floors and mixed precision. | `complete-design-v1/polytropic-bounds-evidence.json` |
| Rational-root classification | 2,157 assertions. Repeated exact square roots decide binary rational roots without forming the numerator power. | `complete-design-v1/rational-root-evidence.json` |
| Scaled integer powers | 820 assertions, including exact binary-power controls with integer exponents of magnitude `10^300`. | `complete-design-v1/scaled-integer-power-evidence.json` |
| Combined interval v2 | All 909 power and 111 integral assertions remain green; 24 concrete compiler checks and 2,042 parent-checked threaded assertions pass. All four worker threads participate. | `complete-design-v2/integrated-power-evidence.json`, `quality-evidence.json` |
| Exact gamma-one coefficient | 466 assertions for complete rational `pe0/n0`, mixed precision, density floors, dimensions and range restoration. The existing exact logarithm arithmetic is unchanged; its temperature type additionally accepts `Rational{BigInt}`. | `/tmp/hybrid-polytropic-log-energy-v1/corrected-evidence.json` |

The abbreviated component directories above have prefix `/tmp/hybrid-polytropic-`.
Their physical-input scope is Float32/64 arrays and closure coefficients. V2 uses
scaled integer powers when rational-root classification succeeds, and retains
the bounded logarithm/exponential path for irrational powers.

A derivation in `complete-design-v1/positive_power_rationality.md` shows that a
positive sum containing any irrational rational-power term cannot equal a
rational midpoint. It uses the fixed-field equivalence in Theorem 3.10(a,b) of
[J. S. Milne’s Fields and Galois Theory](https://www.jmilne.org/math/CourseNotes/FT.pdf),
checked on October 18, plus equality in the complex triangle inequality.
This is a derivation from the theorem, not a statement quoted from that source.
It does not settle the all-rational sum case.

## Retained oracle failure

The first gamma-one test run has 420 passes and 12 failures on minimum Julia;
latest Julia is not run after that failure. A separate 60-assertion classification
checks every failing fixture: the complete integral is exactly zero by a binary-
log identity, but independently bounded per-cell logs leave a nonzero-width
interval larger than the Float64 zero rounding cell. The corrected oracle
combines the exact `log(2)` coefficient for binary powers. Both corrected lanes
pass all 466 assertions with unchanged candidate arithmetic. Original source,
failed log and classification remain in `/tmp/hybrid-polytropic-log-energy-v1`.

## Exact rational-sum equality prototype

`/tmp/hybrid-polytropic-sum-equality-v1` passes 988 numerical/resource assertions
on each native runtime. It clears rational denominators, evaluates the integer
difference modulo distinct deterministically checked primes, and proves equality
once the divisor product exceeds a strict numerator bound. A nonzero residue
proves inequality. Uniform inputs reuse the checked pointwise equality owner.
The derivation is in `brief.md`; no production energy integration is included.

The tests include exact Float32/64 midpoint sums and 997-bit exponents, with
independent direct-power or monotonic controls. The 4,096-cell repeated tie uses
8,755,912 allocated bytes; measured elapsed times are 0.553 s and 0.171 s.
The first minimum-version quality run passes six compiler checks and 134
threaded assertions but fails 12 assertions that assumed worker IDs were 1:4.
The observed fourth default worker has ID 5. The corrected harness checks
actual default-pool membership and all four workers; it passes six compiler and
147 parent-threaded checks on each native runtime. All four default workers
participate; complete evidence is retained in `native-evidence.json`.
The candidate arithmetic is unchanged and the failed run is retained in
`quality-failure.md`. This component has a finite exact decision
procedure, but large-grid and difficult-tie costs remain unaccepted. It does not
yet close integrated-energy rounding or resource requirements.

## Current energy-owner recheck after pressure correction

The energy diagnostic body used by the new witness is byte-identical to
`0d8c8e8`. Both native runtimes pass 34 witness assertions: four exact dyadic
integrals remain finite while the public diagnostic returns zero or infinity,
and six mixed-precision stored-gamma cases still take the wrong branch. This
confirms that pointwise pressure acceptance does not close ENERGY-004/007.
The expected Float32 magnitudes are `2^-100` and `2^100`; Float64 uses `2^-600`
and `2^600`, with inverse box extents. Per-cell pressure rounds before the volume
can restore the finite complete integral.

A separate 216-assertion observation per runtime records zero/sign behavior over
108 combinations of precision, cell count, coefficient sign, gamma and density.
Sources, exact commands, source/log hashes and results are in
`/tmp/hybrid-polytropic-energy-contract-20261018/evidence.json`.
`design.md` there records the next canonical integration strategy and its
unverified fast-path and resource requirements. No complete integral code is
included in the pressure push.

A further 20-assertion witness on each runtime constructs three complete
integrals strictly above a Float64 midpoint. The actual diagnostic returns the
lower float in all three cases. Stored gamma is `2^52`, `2^54` or `2^56`; a
positive density contribution rounds to zero pointwise. Exact coefficient
identities and a rational upper bound on the positive increment establish the
expected adjacent float without constructing the enormous power. Evidence is
`tiny-term-evidence.json` in the current energy-contract directory.

The isolated strict-scaled-sum component retains the current numeric enclosure
and additionally records strict inequalities from quantization remainders.
It passes 1,509 numerical/resource assertions per runtime, including exact
midpoint decisions at zero, subnormal, normal and overflow boundaries. Tested
exponents up to 997 bits allocate 5,152–5,296 bytes. Its first compiler harness
errors because `typeof` widened a singleton type to `DataType`; the corrected
harness preserves exact argument types and retains every analysis. It passes
10 compiler/signature assertions and 99 parent-checked concurrency assertions
on each runtime, with all four default-pool workers participating. Sources,
failed/corrected harnesses and logs are in
`/tmp/hybrid-polytropic-strict-scaled-sum-v1`. This remains a component and is
not included in the pressure commit or an accepted complete-energy repair.


The isolated complete cold-integral component passes 273 assertions on each
runtime: independent rational powers, range restoration, stored gamma, exact
positive/negative midpoints, the tiny-term witnesses and measured zero signs.
The three huge-gamma midpoint samples allocate 87,464–96,336 bytes; these are
small-grid measurements, not a large-grid resource acceptance. The component
shares the pressure power/rounding owners and the checked strict-sum/equality
components. It also passes 12 concrete compiler checks and 195 parent-checked threaded
assertions per runtime. Gamma one, ordinary-path performance and public
integration remain unimplemented in this cold component.
Evidence: `/tmp/hybrid-polytropic-complete-rounding-v1/numerical-evidence.json`.


The separate ordinary-integral probe passes 1,429 numerical assertions per
runtime against exact rational and independent binomial integrals. Its range
certificate bounds the complete positive integral using exact stored-input
binary exponents; it does not assume that finite intermediates certify the
complete result. The selected ordinary accuracy target is four ulps, without a
universal correct-rounding claim. Each runtime also passes 24 concrete compiler,
36 resource and 291 parent-checked concurrency assertions. All 54 measured calls
allocate zero bytes across uniform/varied/weak profiles and 32–262,144 cells.
Timings and combined/public integration remain open. Source, proof and logs:
`/tmp/hybrid-polytropic-energy-ordinary-v1`.

A fixed-width scaled-rational converter passes 7,836 independent rounding
assertions per runtime, including signed Int128 limits, UInt64 denominators,
normal/subnormal/overflow ties and extreme integer exponents. It is intended to
apply the gamma-one reference density before rounding the shared near-unit
isothermal enclosure. Both runtimes also pass two concrete compiler calls,
64 resource assertions (192 measured calls at zero bytes) and 194 threaded
assertions. No production converter has been changed. Source and logs:
`/tmp/hybrid-scaled-rational-conversion-v1`.

The shared gamma-one integration v1 fails 12 strict-rounding assertions on the
minimum runtime (420 pass in that testset). That run is retained. V2 preserves
every test and requires both general-probe endpoints to round identically
before accepting that probe; otherwise it uses the shared exact logarithmic
fallback. Classification records errors of two, four and fourteen ulps in V1;
six cases also fail the additional four-ulp hypothesis (67 assertions pass,
six fail). The original range intervals contain the independent reference in
all twelve cases. That failed classification remains retained.

V2 passes the original 466 strict output assertions, but six near-unit probes
return unresolved. A separate 31-assertion classification verifies that their
stored pressure/reference coefficients are equal and that exact cancellation
before multiplication settles every case. V3 applies that general algebraic
identity. It passes all 826 original/additional numerical assertions, 12 concrete
compiler calls, 72 resource assertions and 67 parent-checked concurrency checks
on each runtime. All 108 measured near-unit calls allocate zero bytes, through
262,144 cells. General Float64 resource measurements are recorded below. Sources,
failed trials and logs: `/tmp/hybrid-polytropic-gamma-one-shared-v{1,2,3}`.


The ordinary V1 paired benchmark passes 252 assertions per runtime but measures
uniform-state latency at 8.6–36.8 times the current public diagnostic. V2 reuses
consecutive positive equal-density pressure terms while keeping every per-cell
check and compensated addition. Both runtimes pass 3,433 numerical/equivalence
and 351 quality/resource/threaded assertions. Its direct 11-pair comparison
passes 184 assertions per runtime; uniform latency is 0.032–0.065 of V1, with
identical results and 88 zero-byte calls per implementation/runtime. Varied
profiles remain mixed, and this is not a universal speedup claim. Full samples:
`/tmp/hybrid-polytropic-energy-ordinary-v2/paired-reuse-evidence.json`.

The gamma-one V3 general-state resource run passes 28 independent-oracle and
immutability checks per runtime. On 1.12.7, 4,096 uniform/alternating cells
allocate 6,658,784–7,222,024 bytes; 32 distinct Float64 densities allocate
1,330,976 bytes. These costs remain open. The near-unit zero-byte measurements
do not establish acceptable general-state performance.

An isolated canonical integration now exists at
`/tmp/hybrid-polytropic-energy-public-v1`. Six production Julia paths change;
shared conversion, isothermal and scaled-sum owners are extended, and a new
polytropic energy owner connects to the public diagnostic after validation.
Nonfinite and generic public branches retain the previous path. Permanent
component and public tests are added, with one shared binomial test oracle.
The initial migration check catches two error strings altered by textual
renaming; restoring them passes all 38 helper comparisons. The initial smoke
harness then errors on `g.ncells`; the reported real field is `g.n`, and the
harness is corrected. Both failed runs are retained. The final migration check
passes 41 assertions and both public smoke runs pass 25 assertions. V1 passes
1,679 concrete/quality and 1,608 threaded assertions per runtime. Its focused
minimum-version suite records 69,720 passes and one failed legacy expectation:
the stored-gamma witness still required the old narrowed-gamma answer. V2 retains
that fixture and checks independent rational concavity bounds for its actual
stored exponent; it also checks that the old answer is rejected. No tolerance
is widened. V2 passes 69,903 focused assertions on each native runtime.
Remaining integrated acceptance is pending. No integrated energy code is in main.

Gamma-one V4 counts repeated positive densities before exact rational log work.
It passes 1,006 numerical, 151 quality/resource/threaded and 28 general-resource
assertions per runtime. At 4,096 uniform/alternating cells it allocates
53,408–84,552 bytes across the two runtimes. At 32 distinct Float64 densities it
allocates 1,341,048 bytes on 1.12.7 and 1,337,920 on 1.13.2; distinct-density
costs remain open. The near-unit path retains its zero-byte controls. Source,
complete measurements and failed predecessors are in
`/tmp/hybrid-polytropic-gamma-one-shared-v4`. Public V2 includes this shared-owner
change and 180 independent repeated/distinct-density assertions.

## Algebraic storage dispatch — CLOSURE-002

**VERIFIED:** the pushed pointwise-pressure path fails for union-element input
or output storage. The corrected observation harness in
`/tmp/hybrid-algebraic-package-review-20261018/union_storage_corrected.jl` passes
69 assertions on each native Julia version: the legacy broadcast succeeds in
all 27 combinations of Float32, Float64 and union input/output/coefficient
types, while the new owner raises MethodError in the 15 combinations with union
input or output storage. `src/electrons/polytropic_pressure.jl` selects a cache
that calls `zero(T)`/`zero(S)` for those nonconcrete types. The original harness's
unrelated union-output initialization error is retained separately.

The isolated correction at `/tmp/hybrid-algebraic-dispatch-v1` restores legacy
broadcast for nonconcrete input/output element types. It exposes the actual
Array check at the cache call site and retains the scalar type explicitly in
raw shock broadcasts. It also generates the two selected Int64/Int128
isothermal polynomial methods from one arithmetic body. Full-package analysis
of the integrated V1 snapshot records six new diagnostic groups on 1.12.7
(22 to 28 reports): two mixed-width polynomial findings, one hidden-Array-check
finding and three raw-broadcast findings. These are not waived. The report
comparison harness initially fails on a new source file; its corrected
missing-file map preserves every report and marks new source locations for
review. Both the failure and complete corrected trees are retained in
`/tmp/hybrid-algebraic-package-review-20261018`.

V1 passes 274 union Array/view/empty-storage assertions per runtime. Its first
raw-pressure test extension incorrectly uses views for an Array-only helper;
that failed harness is retained. Corrected fixtures pass another 42 assertions
per runtime for shared Array storage, self-alias, broadcast expansion, signed
and nonfinite inputs, empty arrays and explicit rejection of unsupported views.
Full compiler comparison still fails with two new polynomial diagnostic groups
on each runtime (24/23 versus 22/21 baseline); all other groups match. V2 exposes
the existing concrete coefficient-tuple contract at its two indexing sites.
Thirteen syntax/method assertions verify the two integer widths, those exact
tuple contracts and the unchanged arithmetic body. Formatting preserves all
221 source/test executable syntax trees. V2 numerical, full compiler,
concurrency and paired allocation/latency gates remain running. Verification
logs and immutable hashes are in
`/tmp/hybrid-algebraic-dispatch-v{1,2}-verification`.

Energy public V3 merges that dispatch correction into V2, retaining all prior
energy fixtures and adding the 316 storage/broadcast checks. All 41 helper
migration assertions pass, including explicit validation of the generated
polynomial methods. Formatting preserves all 231 source/test executable trees.
Focused, concrete, threaded and complete compiler checks are running outside
main at `/tmp/hybrid-polytropic-energy-public-v3-verification`.

## Remaining acceptance work

1. Finish public pressure v2's native solver/restart, compiler-tree, quality and
   threaded gates and its exact-revision GitHub CI. Both remain required after
   the separate pressure push. Preserve the held radiation files.
2. Complete integrated polytropic final rounding and termination, including
   all-rational multi-term midpoint decisions and resource bounds for large
   integer powers. Pointwise pressure success does not settle this work.
3. Integrate exact gamma-one handling and preserved stored gamma with an accurate
   ordinary energy path. Check nonfinite propagation, generic dispatch,
   validation order, complete budgets and actual large-grid memory/latency.

The final dispatch correction is implemented with the two CI test repairs. Both native complete compiler
trees match baseline (22/21); concrete compiler gates pass 316/300 assertions,
four-thread checks pass 3,007 per runtime, and Aqua passes nine each. Paired
polynomial checks pass 528 per runtime, preserving outputs and 72 zero-byte
calls per implementation/runtime. Timing ratios are mixed (0.776–1.617 on
1.12.7; 0.507–4.383 on 1.13.2); no speedup is claimed. Each native focused suite
passes 57,165 assertions. The final
six-file source/test overlay exactly matches the checked snapshots, and all
other 217 source/test paths match `0d8c8e8`. Whole-repository formatting passes.
Remaining native broad gates and the new exact-SHA CI remain required after
push. Evidence: `/tmp/hybrid-algebraic-ci-promotion-v1-verification`.

The reviewed dispatch/CI overlay is pushed as `15662e1`; its
[exact-SHA GitHub run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39137529999)
is running. Integrated polytropic energy remains isolated and unpromoted.

Both dispatch focused suites now pass 57,165 assertions. The recovered
isothermal V2 native lanes each pass 44,809 assertions across checkpoints,
all eleven solver files, quality and concurrency. Their complete compiler
comparisons retain the same two polynomial findings already corrected by
`15662e1`; those original failed comparisons remain preserved.
The remaining 22 solver/shock/checkpoint files are running against the final
pushed source in `/tmp/hybrid-algebraic-final-native-verification`, with both
oracle modules loaded verbatim. No source/test is changed during these runs.

Energy V3 passes 70,219 focused assertions on each runtime; both complete
compiler trees match baseline, and all concrete and concurrency tasks pass. V4 rebases exactly 18 energy source/test paths onto `15662e1`, retaining
both CI fixture corrections; those energy files remain byte-identical to V3
after a full formatter scan. V4 is not promoted. Distinct-density logarithmic
resource costs and final native broad acceptance remain open.

The first public-resource harness records 12 zero-byte failures caused by its
heterogeneous measurement context; the concrete typed/warmed replay passes
all 108 assertions on each runtime, including every zero-byte control. Both
runs are retained. The typed 1.12 replay confirms 2,272,925,040 cumulative
allocated bytes for 4,096 distinct Float64 densities, not peak resident memory.
Energy promotion remains open on this measured resource cost.

An isolated fixed-width logarithm prototype in
`/tmp/hybrid-isothermal-fixed-log-v1` now passes 15,086 independent assertions
per runtime: 13,041 exact wide-product/fraction/rounding checks and 2,045 directed
logarithm enclosure checks across the binary input range. Its 34-term atanh
remainder bound is verified by an exact integer inequality. It reuses the
canonical dyadic sum, coefficient and reference-rounding owners, accepting a
probe only when both endpoints round identically; unresolved cases retain the
exact owner. Existing gamma-one fixtures and public resource gates are running.
No prototype logarithm code is promoted.

The fixed-log prototype now also passes the 862 retained gamma-one assertions,
120 public numerical/resource assertions and 223 concrete/resource/concurrency
assertions on each runtime. All 24 measured public profiles per runtime allocate
zero bytes, including 4,096 distinct densities in Float32 and Float64. These
measurements establish those cases, not a universal constant-cost exact fallback.
The unmodified exact evaluator still handles unresolved certificates and exact
cancellation. Evidence: `integration-results.json`, `quality-results.json` and
associated logs in `/tmp/hybrid-isothermal-fixed-log-v1`.

Public energy V5 integrates the six fixed-log helpers into the shared isothermal
arithmetic file and inserts their probe into the canonical gamma-one owner.
There is no prototype wrapper in production. Its brief and gates are in
`/tmp/hybrid-polytropic-energy-public-v5-verification`. The corrected formatter
verifier passes all 233 source/test syntax comparisons and independently checks
the three original unformatted inputs affected by line annotations. Original
harness failures are retained. Numerical, compiler and concurrency acceptance
must complete before promotion.

The final V5 migration checks pass all 45 assertions; the complete Julia 1.12
compiler tree matches baseline at 22 reports with no added/removed diagnostic
groups. Both public smoke runs pass 25 assertions. The full focused, latest
compiler and quality/concurrency gates continue. A separate public audit uses
512-bit directed MPFR bounds across mixed coefficient precisions, dimensions,
extreme densities and cell spacings; its results remain pending. The 22-file
native integration lanes are queued after all core gates succeed.

Both V5 complete compiler trees now match baseline (22/21 reports), and the
independent public logarithmic audit passes 577 assertions per runtime. All
144 sampled combinations per runtime retain the declared rounding and sign
behavior against independent directed MPFR bounds, without input mutation or
ambient precision changes. These are sampled checks, not an exhaustive proof.

V5's Julia 1.12.7 focused suite passes 85,433 assertions. Both runtimes pass
1,679 quality/concrete/Aqua assertions, 1,608 threaded assertions and 223
additional fixed-log compiler/resource/concurrency checks. The latter includes
explicit unresolved-cancellation and invalid-input controls. Latest focused
and subsequent native solver/shock/checkpoint acceptance remain pending.

V5's final focused suites both pass 85,433 assertions (9m33.8s on 1.12.7;
6m44.6s on 1.13.2). All core numerical, migration, complete compiler, quality,
concurrency, independent public range and read-only formatting gates pass.
The reviewed 20-file energy overlay is copied to main with matching SHA256s;
all three held radiation files retain their hashes. Native 22-file broad lanes
have started; their results and exact-revision GitHub acceptance remain required.
Promotion manifest: `/tmp/hybrid-polytropic-energy-public-v5-verification/promotion-plan.json`.

The verified overlay is pushed as `7104c6190b09daaa01349d7dfd4d299ac9571da3`. Its [exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39141208484) is running, with formatting and all MPI ranks green. The committed 24-path whitelist and remote SHA are verified in `promotion-15662e1-whitelist.json`.

The final V5 native 1.12.7 solver, shock and checkpoint lanes all finish with
exit code zero. The 1.13.2 lane follows under the same frozen-source launcher;
its final result remains required. Evidence is retained in
`/tmp/hybrid-polytropic-energy-public-v5-broad-verification/results112.json`.
