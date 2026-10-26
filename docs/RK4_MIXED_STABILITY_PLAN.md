# Joint Hall-wave and resistive RK4 timestep

## Contract and repair brief

`recommended_dt` owns the periodic hybrid magnetic timestep estimate. Its
RK4 estimate must bound the combined mode `z = (-ν + iω)Δt_B`, with
`ν <= η k_max² + ηH k_max⁴`, `|ω| <= ω_W(k_max)` and `Δt_B = Δt_p/NB`.
This is a frozen linear magnetic stability estimate. It does not certify
nonlinear evolution, pressure advection/compression or ion-acoustic accuracy.
The canonical implementation is `src/integrators/abstract_integrator.jl`;
the two affected stage owners are `_rk4_B!` and `_rk4_Bp!` in
`src/integrators/extrapolated_leapfrog.jl`.

**Mode:** repair. **Baseline:** on a 33-cell, length-2π grid, a circular k=16
Hall mode with η≈0.998283 grows by 1.4668 at the old default recommendation.
The independent RK4 polynomial and production `_rk4_Bp!` agree on Julia
1.12.7 and 1.13.3, for Float32 and Float64. Evidence is retained in
`takeover/rk4-mixed-stability/before-v1-*` under the dated external audit
archive, with pre-edit source pins and complete witnesses.

**Verified cause:** the old helper takes the minimum of separate imaginary-
and negative-real-axis caps. Their intersection includes points outside the
RK4 stability region. For Float64 the accepted mode has
`z = -2.0950156106519744 + 2.0986190386582613im` and polynomial magnitude
`1.4668001004038733`; the production update measures `1.4668001004038727`.

**Compatibility:** retain the public arguments, validation, gyro cap, subcycle
meaning, log-space range handling and all non-RK4 paths. Zero damping retains
the original arithmetic exactly. Damped RK4 recommendations change; the old
test equality to the separate-axis minimum is replaced by a geometric oracle
for the corrected contract, retaining its parameter coverage and tolerance.
There is no state or archive migration and no new dependency.

## Certified inner polygon

For `a = νΔt_B >= 0`, `b = |ω|Δt_B >= 0`, use the origin and vertices

`(0, 2.8), (1, 2.5), (1.85, 1.85), (2.5, 1), (2.785, 0)`.

Its four edges are

| Constraint | Rate bound before multiplying by safety and NB |
| --- | --- |
| `0.3a + b <= 2.8` | `2.8 / (0.3ν + ω)` |
| `0.65a + 0.85b <= 2.775` | `2.775 / (0.65ν + 0.85ω)` |
| `0.85a + 0.65b <= 2.775` | `2.775 / (0.85ν + 0.65ω)` |
| `a + 0.285b <= 2.785` | `2.785 / (ν + 0.285ω)` |

All normals are positive, making the polygon downward closed: bounding the
two rates independently bounds their whole rectangle. Scaling by `safety`
stays within the polygon. Pure-axis intercepts retain the existing 2.8 and
2.785 constants. At equal rates the new boundary is 1.85 per coordinate,
compared with about 1.3925 for the simpler triangle joining only the axes.
This is an inner bound, not the optimal RK4 stability boundary.

`test/test_rk4_mixed_stability.jl` constructs the exact classical polynomial
`R(z) = 1 + z + z²/2 + z³/6 + z⁴/24`. On each triangle from the origin to
adjacent vertices, it expresses `|R(-a+ib)|² - 1` in degree-16 triangular
Bernstein form. Every coefficient is nonpositive using `Rational{BigInt}`
arithmetic. The Bernstein basis is nonnegative and sums to one on a triangle,
so this certifies every point, rather than sampling its boundary. The test
certifies vertices enlarged by `1 + 2^-20` as well, leaving numerical margin
for coefficient conversion and the log-space evaluation. The final certificate
has 624 assertions on both supported runtimes, including
eight checks of the positive-normal/downward-closure geometry. The initial
616-assertion polynomial certificate passed before the production repair.

The helper normalizes both rates by their maximum logarithm. Each normalized
rate lies in [0,1], at least one equals one, and each denominator remains
bounded and positive. It takes the minimum of the four ratios and returns
its logarithm minus the scale. Thus neither an overflowing physical rate nor
an underflowing reciprocal is formed. Existing final representability and
gyro-limit semantics remain in the public helper.

## Acceptance tracker

- [x] Preserve baseline scalar-polynomial and production-stage failures on
  actual Julia 1.12.7 and 1.13.3.
- [x] Exact polynomial certificate, including the enlarged polygon.
- [x] Frozen V2 source: JuliaFormatter 1.0.62 preserves syntax and is idempotent.
- [x] New rectangle sweep, both RK4 stage owners, original helper/scale tests,
  and existing transport/coupled stability regressions on both runtimes.
- [x] Concrete analysis, Aqua, full compiler-tree comparison and matched
  allocation profiles; no diagnostic suppressions or timing claim.
- [x] Root registry-only integration on Julia 1.12.7 and 1.13.3: all six
  focused/resource/compiler commands pass; executable files match frozen V2.
- [ ] Separate main push and all six exact-SHA CI jobs.

Adversarial checks distinguish separate-axis minima, a damping-only cap,
omitted ηH or NB, a blanket smaller timestep and a repair confined to one
stage owner. Exact certificates and the independent BigFloat geometry oracle
complement the mode rectangle and production tests. Failed proof experiments
and test harnesses remain in the evidence archive; they are not accepted
production candidates.

V2 passes 23,432 focused assertions per runtime. All 48 matched profiles
preserve three warmed allocation samples at 0 B; undamped and non-RK4
recommendations are bitwise identical to the frozen parent. V1 changed the
association of one non-RK4 logarithmic sum, producing a last-bit compatibility
regression. Its failed receipts remain retained; V2 restores the original
operation order. Full 29/27-report compiler aggregates, all 252/203 virtual
frame records and messages match the parent after exact source-root replacement
only. Twenty-four concrete analyses and Aqua pass per runtime.

An independent 144-case geometric oracle accepts V2 and rejects each of five
mutants: separate axes, damping only, omitted hyperresistivity, omitted NB and
a blanket half-step. Original witnesses, source freezes, failed V1 receipts,
V2 terminal commands, compiler review and mutant outputs are retained under
`takeover/rk4-mixed-stability`. No latency improvement is claimed.

A separate pre-existing caller-context defect remains open: at BigFloat
precision 11, η=10^200 on a 129-cell grid makes the original recommendation
grow an actual RK4 mode by 1.296709 on both runtimes. The current low-precision
BigFloat branch is not certified by the floating-point margin above. Its
precision, rounding, exponent-range and return-type repair has a separate
owner and acceptance campaign; this change does not close that defect.
