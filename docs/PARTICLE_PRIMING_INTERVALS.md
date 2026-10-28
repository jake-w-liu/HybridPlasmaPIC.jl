# Residual-certified relativistic priming

Relativistic priming stores a velocity that the following Boris step converts
back to momentum. Near light speed, a one-ulp velocity error can appreciably
change that next conversion. The native priming guard therefore requires
a resolved input Lorentz factor and a rounded native-output normalized squared
speed `<=1/2`. This is a comparison of the computed native output, including
its rounded normalization; it is not an exact physical speed classifier.
Other IEEE cases enter the
existing repair owner, which first attempts the bounded interval evaluation
below, then retains the existing BigFloat replay if the interval is unresolved.

The physical expression is unchanged:

```
gamma0 = c / sqrt(c² - dot(v,v))
m = gamma0*v - h*qm*(E + cross(v,B))
gamma = sqrt(1 + dot(m,m)/c²)
primed = m/gamma
```

The center arithmetic is Float64 double-word arithmetic from the shared
`src/utils/error_free_arithmetic.jl` owner. Each center is accompanied by a
Float64 absolute error radius in `src/particles/relativistic_priming_dd.jl`.
The certificate uses exact residual identities rather than a fitted global
DD error coefficient. It does not require an assumed normalization or error
theorem for the DD center routines: an inaccurate center widens its interval.

## Arithmetic assumptions and bounds

This is a CPU IEEE round-to-nearest calculation with gradual underflow,
explicit single-rounded FMA and no unsafe reassociation. For nonnegative
bounds, one `nextfloat` after addition/multiplication/division encloses the
exact operation above; one `prevfloat` encloses a subtraction below. Infinite
or negative radii are unresolved.

TwoSum and grow-expansion preserve the exact sum of finite Float64 words
without intermediate overflow. Absolute-summing every expansion term with
outward rounding bounds the true residual. No expansion ordering/nonoverlap
claim is needed for that triangle inequality. A nonfinite word rejects the
bound. TwoProduct's rounded FMA residual can lose at most one
minimum-subnormal quantum; four such losses cover all four products of two
word-pair sums, including the low-by-low product omitted by the DD center.

Let `A=a+da`, `B=b+db`, with `|da|<=ra`, `|db|<=rb`, and let `z` be the exact
sum of a computed center pair. All error-bound arithmetic rounds outward.
`U(a)` bounds the magnitude of pair a above; `L(a)>0` bounds it below.

| Operation | Residual and propagated error |
| --- | --- |
| `A+B` | Expand `a+b-z`; bound by `|residual|+ra+rb`. |
| `A*B` | Expand four products minus z; bound by `|residual|+U(a)rb+U(b)ra+ra*rb`. |
| `A/B` | Require positive lower bound l on `|B|`. Expand `a-b*z`; use `(A/B)-z=(A-B*z)/B`, giving radius `(|residual|+ra+rb*U(z))/l`. |
| `sqrt(A)` | Require A>0 and positive center z. Expand `a-z²`; use `sqrt(A)-z=(A-z²)/(sqrt(A)+z)`. Since the denominator is at least `z>=L(z)>0`, radius is `(|residual|+ra)/L(z)`. |

The division identity includes the uncertain final gamma, so momentum
perturbations are propagated through normalization. Structural-zero outputs
are admitted only when each source factor proves the exact component zero.

## Admission and publication

Every nonzero stored operand must be between `2^-100` and `2^100`. The
positive sum of momentum-term scales must be finite and at most `2^300`.
These limits bound working range, without limiting physical speed or clipping
any value. The canonical scaled Lorentz-gap owner exposes the same square
words already used by Boris; priming encloses their residual before the
initial sqrt/division. Nonpositive or uncertain denominators, nonfinite
centers/radii and unresolved gaps reject the fast attempt.

For each final component, neighboring T values define its proposed rounded
value's half spacings. The component interval must lie strictly inside that
one rounding cell. The subtraction from its center has its own exact-residual
bound. Strict comparisons remain conservative when the cell endpoints round:
a representable Float64 number strictly below RN(boundary) is also below the
true boundary. Ties, ambiguous cells and underflowed Float64 half spacings
retain BigFloat replay. Every admitted fast output is thus certified to round
to the proposed T component under the arithmetic assumptions above.

## Compatibility and verification

The original 100,000-input fixtures for every IEEE type and both priming
owners are retained. An independent 4096-bit reference classifies physical
output squared speed/c²: the well-conditioned slice `<=1/2` retains historical
bits; the repaired slice `>1/2` equals the independently rounded stage. This
changes the former blanket priming bit-equality requirement. The old failing
fixtures and optimized candidates remain in the takeover evidence. Whole-step
input coverage and tolerances were not weakened.

Permanent controls test exact rational interval endpoints, underflow,
cancellation, uncertain denominators, ambiguous rounding cells, independently
rounded stages and actual EMPIC full steps. The broader audit also checks
minimum/latest Julia, concurrent caller BigFloat scopes, actual Metal,
compiler reports and warm allocations. See
[KINEMATICS_AUDIT_20261027.md](KINEMATICS_AUDIT_20261027.md) for current gate
status and limitations. A finite corpus and source review are not a formal
machine-checked proof. Accurate component rounding alone can put a limiting
velocity on the public strict-cone boundary. The separate
[publication repair](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md) keeps that domain
and rejects such outputs; it does not move them inward or preserve an otherwise
unrepresentable physical Lorentz factor. Its acceptance status is recorded in
that plan.

A separate eight-step scalar magnetic-field stress audit (`Bz=1/sqrt(eps(T))`,
`dt=0.125`, loaded speeds 0.9c/0.99c followed by priming) distinguishes those
scopes:
V5 passed all 32 local actual-stored-input four-ulp comparisons, but 28 points
differed from a sequential independently rounded trajectory by more than an
exploratory accumulated 64 ulp/step envelope. The cold replay matched that
trajectory. This is a material multistep conditioning limitation; the local
contract and first-step witnesses do not establish a trajectory-wide ulp
guarantee. Defining a stronger trajectory/state contract remains open.

At the threshold, exact physical and rounded native classifications can differ
by roundoff. For example, Float64 priming from rest with E=(-0.6,-0.8,0),
h=qm=c=1 and B=0 has exact output squared speed/c² slightly above 1/2 but
native rounded value 0.49999999999999994. V5 keeps the historical x component
0.42426406871192845, one ulp from its independently rounded value. This
well-conditioned boundary slice retains the local normwise contract; it is
not covered by a blanket exact-physical-speed correct-rounding guarantee.
The retained 100k fixture classification is still independent and unchanged.
