# Polytropic energy range and closure precision

Status: the pointwise pressure correction passes native focused and concrete
compiler checks on Julia 1.12.7 and 1.13.2 and is entering exact-revision CI.
Broader native gates remain pending. The complete integrated polytropic energy
repair is not implemented. This plan covers the polytropic part of ENERGY-004,
ENERGY-007 and CLOSURE-001.

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
regressions below, is being pushed to run exact-SHA GitHub CI concurrently.
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
