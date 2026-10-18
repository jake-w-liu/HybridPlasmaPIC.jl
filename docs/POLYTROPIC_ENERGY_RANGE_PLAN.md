# Polytropic energy range and closure precision

Status: confirmed defects; isolated arithmetic components pass on Julia 1.12.7
and 1.13.2. The complete public repair is not implemented. This covers the
polytropic part of ENERGY-004, ENERGY-007 and CLOSURE-001.

## Confirmed behavior and implementation contract

**VERIFIED — ENERGY-007:** a Float64 gamma immediately above or below one
becomes exactly one when converted to Float32. The diagnostic then takes its
pressure-weighted logarithmic branch. At unit density/reference density and
`pe0 = 0.5`, it returns zero instead of the nonzero stored-closure integral.
The independent oracle is rational because `(n/n0)^gamma = 1` exactly.

The six 1D–3D witnesses have finite nonzero reference values from approximately
`4.50e15` to `3.60e16` in magnitude. Both runtimes pass 40 witness assertions,
including Float64-grid and exact-gamma controls. The diagnostic body is
byte-identical to `63f2420`: narrowing and branch selection are at
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
directory. The same helper arithmetic must also be reviewed at the separate
`PerpShock` expression in `src/initial_conditions/reflecting_wall_shock.jl:250`–`:251`;
that caller has not yet received a public regression or repair.

The canonical owner is `src/electrons/closures.jl:85`. Its result feeds
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

## Remaining acceptance work

1. Complete final rounding and termination, including all-rational midpoint
   decisions and bounds on resource use for large integer powers.
2. Integrate exact gamma-one handling with an accurate ordinary path; check
   nonfinite propagation, generic dispatch, validation order and complete budgets.
3. Measure large-grid memory and latency. The interval and logarithmic cold
   components are not accepted as the ordinary public traversal.
4. Complete both native solver/restart, compiler-tree, quality and threaded
   gates. Commit/push the verified public delta and require green GitHub CI for
   that exact SHA. Preserve the separate held radiation changes.
