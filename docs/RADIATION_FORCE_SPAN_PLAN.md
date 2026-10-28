# Radiation force-span completion

## Pre-edit brief, 2026-10-27

- Mode: repair the canonical frozen-rate radiation evaluation.
- VERIFIED requirement failure: both Julia 1.12.7 and 1.13.3 fail six of seven retained exponent-span witnesses in takeover/radiation-span-before-*.log. With v=(1/2,0,0), E=(2^-e,2^(e-1),0), B=(0,0,2^e), K=dt=2^e, c=1, the exact large transverse force cancels and K*dt*chi²=3/4. The result must be approximately `0.23598626362439432` for every exponent. Baseline output is `0.5` in the failing cases.
- Cause: _radiation_scaled_factor rescales all force components together and squares the remaining small components. Those squares can underflow to zero; the early chi2V4<=0 return loses a physically finite optical depth after K*dt amplification. A small fraction of the original force scale is not an error bound for the final velocity.
- Canonical owners: radiation_reaction.jl owns force/invariant/factor evaluation; exact_arithmetic.jl already owns correctly rounded rational and square-root conversions. Reuse those cold conversion owners. No new public API, dependency or second radiation model.
- Contract: the same frozen-rate substep, finite converted fields and subluminal input; hardware velocity factor in[0,1], normwise few-ulp accuracy; normal paths preserve operation order/allocation; BigFloat/generic policy unchanged. Cold exact work is per call, with no MPFR or global precision/cache changes.
- Strategy: keep the native error certificate. Before uncertified scaled double-double arithmetic, bound its field-scale error propagated through K*dt*gamma. A conservative2^-80 normalized-invariant envelope (106-bit double-double, fewer than2^24 bounded elementary-error coefficients) is below2^-56 when its exponent is at most24. Complete larger amplifications or an unresolved Lorentz gap from exact rational stored operands. Round the positive optical-depth square root with the existing exact owner, retain its residual as a second word, and evaluate the same bounded velocity factor using normalized rational r and a two-word log2 when needed. No silent force clipping or unsupported-input narrowing.
- Range details: exact r is positive and at most1 after existing validation. For Float64 stored operands r>2^-4196, so alpha>=4096 makes even floatmax(v)*F smaller than the minimum subnormal; alpha<2^-60 changes a velocity by less than the normwise budget. Positive optical-depth values within those endpoints fit Float64; exact residual/(2alpha) supplies the second word. Normalizing r before its logarithm avoids losing an unrepresentable raw ratio. Derive the log2 split against a serial high-precision test before using it.
- Verification: unchanged before witnesses; all axis/sign/field-orientation and neighboring exponent cases with a serial8192-bit momentum oracle; zero-force and weak-drag monotonicity; genuine damping, underflow and log-domain boundaries; ordinary/previous exceptional zero-allocation gates unchanged; measure new cold allocations explicitly; task isolation; concrete/full compiler comparison; formatter and both-version neighbors. Review source and all failure paths before promotion.
- Wrong patches rejected: use the dominant field and discard a small component; return one for zero scaled invariant; clamp K*dt; add a BigFloat fallback with shared precision; widen only the already rounded invariant. Exact original-operand references and task/allocation gates distinguish these.

## Acceptance evidence

The frozen numerical V1 passes 54,863 affected assertions, 108 new four-thread
exact-tier checks, nine clean concrete compiler analyses and Aqua on each Julia
1.12.7 and 1.13.3. Complete local package report trees match the parent after
explicit checkout-root replacement (29/27 existing reports). All earlier 27
ordinary/balanced-exceptional profiles remain allocation-free. The new amplified
exact cases allocate 19,384–23,040 bytes per particle, exactly linear across the
measured 1/16/64-particle calls; these are cold correctness costs, not a speedup.

Independent exact-invariant/8192-bit checks pass 9,728 assertions over 2,432 cases
per runtime; public near-cone checks pass another 1,719. Maximum observed error
is 1.799 ulps. The final receipt pins every raw case, command and source hash in
`radiation-span-independent-audit/final-review-v1.json`. Root verification is in
`takeover/radiation-span/quality-review-v1.json`. All evidence is retained under
`HybridPlasmaPIC-audit-evidence/2026-10-27`.

The integrated source retains V1's executable expression tree; comments, this
proof record and a permanent amplified-path task regression are added separately.
All eight post-integration commands pass: formatter AST/idempotence, executable
AST comparison (documentation strings and line nodes excluded), and both-runtime
focused, quality and cold-path campaigns. The integrated affected suite passes
54,877 assertions per runtime, including the permanent task regression. Nine
concrete calls are clean; complete 29/27 package trees are byte-identical to the
parent without normalization. All 21 cold-resource profiles per runtime match V1
exactly. `takeover/radiation-span/integration-review-v4.json` records these checks;
failed verification-harness attempts remain in `harness-corrections-v1.json`.
`07792ab` passes all six [exact-SHA CI jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40374929835); complete compiler acceptance is recorded in the [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress). The frozen-rate physics approximation and
generic BigFloat policy are unchanged.

The independent hand derivation below is tied to frozen numerical source
`64fb33b5f1eb0d518793304fa03f12e80411b9d6fe0c345103be58475235c442`.
It is not a machine-checked theorem. The unchanged executable expression tree
connects that derivation to the integrated source.

## Shared arithmetic integration

The V5 particle repair moves the pure error-free and double-word routines to
`src/utils/error_free_arithmetic.jl`, under `_eft_*`, `_dd_*` and
`_DoubleFloat64` private names. Radiation retains its complete physical and
admission policy, including the exact-invariant completion and amplification
guard. No compatibility wrappers or duplicated arithmetic bodies are added.

The integration source proof reconstructs all 35 top-level expressions from
radiation commit `07792ab5b7ac5854948fe85fcb6ab196d43a3366` using the 16 moved
definitions and 19 retained policy expressions. After the sealed bijective
symbol mapping and line-node removal, the complete tree is equal. Guard,
product and dropped-expression mutants are rejected. This proof changes no
executable expression or documentation string in radiation. Fresh integrated
checks with registered SpectralOperators 0.1.6 pass the same 54,877 affected,
108 cold resource/task and 11,447 independent assertions on each runtime.
All 2,432 raw independent cases match bitwise and all 21 cold allocation
profiles match the original repair. Remaining integration gates are listed in
[the tracker](KINEMATICS_V5_INTEGRATION_20261027.md); the historical receipts
above retain their original source and dependency scope.

## Conservative arithmetic inventory

This is a hand error derivation tied to the literal frozen operation order,
not a machine-checked floating-point theorem. It assumes round-to-nearest
binary64 with gradual underflow and FMA, the arithmetic used by the candidate.
Let delta=2^-1075 be half the minimum-subnormal spacing.

The double-word addition is the normalized `AccurateDWPlusDW` algorithm. Its
published bound and normalization requirement are documented in the primary
[Boldo/Jeannerod/Melquiond/Muller survey](https://doi.org/10.1017/S0962492922000101),
Algorithms 11/13 and pp. 268–272. I compared the literal addition order with
Algorithm 13. The candidate's nested-FMA product is a different ordering from
the survey's product; its bound must be derived for that ordering.

Conservative building blocks used here are:

- DD+DD: absolute arithmetic error at most 4u²(A+B), plus small additive
  underflow losses, for normalized high-word magnitude bounds A,B.
- DD+scalar: use 4u²(A+B), although the normal-range algorithm has a tighter
  bound. The literal TwoSum/low addition/FastTwoSum order is retained.
- DD*DD: use 8u² AB plus 8delta. The two nested FMAs contribute at most about
  2u² AB and 3u² AB; the omitted low*low contributes at most u² AB. The product
  residual loses at most delta. The retained bound leaves slack for higher
  terms and final normalization. DD*scalar uses 4u² AB plus 4delta.
- Reciprocal of the exact normalized c_hat² in [1,4): use 24u² relative. The
  first quotient has an O(u) residual; the first multiply/subtract leaves at
  most about 13u² absolute residual error. The second residual correction
  cancels the second quotient's own rounding error; the first residual's error
  remains. The final normalization contributes at most 2u² relative. A 24u²
  bound leaves slack for the smaller terms. This bound is for the literal
  three-quotient routine and this bounded reciprocal, not arbitrary division.

These building-block bounds require normalized pairs. The source creates them
from TwoProduct and normalized DD operations; arbitrary two-word inputs are not
covered. Underflowed products/scalings contribute absolute delta terms instead
of an unjustified relative error. There is no overflow inside the bounded
invariant computation.

For the scaled invariant, the exact box is:

| Quantity | Magnitude bound |
| --- | --- |
| each vtilde | <2 |
| each Etilde | <1 |
| each Btilde | <1/2 |
| each w | <3 |
| D | <6 |
| V2 | [1,12) |
| each n=V2*w-D*vtilde | <48 |
| chi2V4 | <7344 |

Using the block bounds and outward slack for higher-order errors gives this
absolute-error coefficient inventory in units of u²:

| Quantity | Coefficient used |
| --- | ---: |
| r | 256 |
| each w | 64 |
| D | 64 |
| V2 | 128 |
| each n | 2048 |
| sum of three n² | <691204 |
| r*D²*V2 | <134799 |
| chi2V4 after their addition | <855379 |

For example, V2*w contributes at most
`12*64 + 3*128 + 8*36 + 1 = 1441`; D*vtilde contributes at most
`2*64 + 4*12 + 1 = 177`; their subtraction adds at most
`4*(36+12)`, keeping n below 2048. Squaring and summing three such n terms
uses `3*(96*2048+8*2304+1) + 4*(4608+6912) < 691204`.

The exact normalized chi² is at most |w|²<27. The V2² denominator has error
coefficient below 4226 and is at least one. Including this denominator, its
normalization, the final division and remaining small terms stays below 2^21
u²=2^-85. A generous 2^24 delta allowance is negligible beside that bound.
This is over five binary orders tighter than the proposed 2^-80 envelope.

For `gain<=24`, allow an extra factor two for the rounded r-high binade when
bounding the actual gamma. The tighter inventory still gives an optical-depth
error below `2^-85 * 2^25 = 2^-60`, leaving margin to the stated 2^-56 budget.
An underflowed numerator/A cannot be made significant: the same gain bounds
its lost absolute contribution, while V2>=1. No clipping of a finite large
K*dt is involved.

## The gap denominator needs its own bound

The invariant budget does not alone justify replacing the Lorentz denominator
by its DD value. For the coarse sequential-gap branch, the error is absolute
O(u²), not uniformly relative O(u²). Three bounded subtractions, the bounded
reciprocal and a multiplication fit the 256u² absolute bound above. The branch
condition ensures the true r is safely above 2^-41.

For `F=(1-r+r*exp(2alpha))^-1/2`, with alpha=kappa/sqrt(r) and kappa held fixed,
write L=log(1/r)/2 and t=r*exp(2alpha). Direct differentiation gives

`dF/dlog(r) = -F^3 * (t*(1-alpha)-r)/2`.

Since the denominator is at least both one and t,
`|dF/dlog(r)| <= 2 + log(1/r)/4`. One proves this by splitting alpha at L+2:
below the split use `F^3*t<=1`; above it use `F^3*t<=exp(L-alpha)`.
Thus even at r=2^-41 the coarse gap error contributes under about 6.4e-17
absolute factor error, below one normwise ulp.

Below the gap switch, the exact nonoverlapping expansion is accumulated from
small to large, avoiding subtraction of an already truncated gap. With at most
eight terms, a 128u² relative budget covers accumulation, reciprocal and
multiplication. Lost square tails are at most a few delta; with r>=2^-800 their
relative contribution is below about 2^-273. The factor sensitivity is below
150 in this range, so this denominator contribution remains far below a
normwise ulp. Smaller/unresolved gaps enter the exact branch.

This derivation does not depend on a stronger optional lower bound for the
smallest physically constructible three-component gap.

## Exact factor and exponential/logarithmic branches

- The exact branch reconstructs c² and K*dt from normalized exact products and
  powers of two, and force/chi² from the original stored operands. For a
  subluminal velocity chi² is nonnegative by its parallel/perpendicular
  decomposition. It does not widen an already rounded force invariant.
- The early endpoints are conservative: |dF/dalpha|<=1 handles alpha<2^-60;
  the loose r>2^-4196 bound makes alpha>=4096 safely underflow even after the
  largest finite velocity is multiplied. These checks use exact alpha².
- The exact square-root residual divided by 2alpha provides the next word;
  the remaining Newton error is O(u² alpha), with alpha in [2^-60,4096).
- In the `2alpha<=708` branch, a subnormal rounding loss in the two-word r
  may be amplified by exp(708), but its absolute factor error remains below
  roughly 2^-54.5, within the normwise budget. A relative-only claim for that
  subnormal case would be wrong.
- The normalized logarithm preserves ratios too small for Float64. Its log2
  split is retained through the exponent multiplication; the only libm input
  is a mantissa in [1,2). The scalar logarithm's existing libm accuracy is an
  inherited arithmetic-library assumption, not an exact-integer theorem.
- Below argument=-60log(2), returning one costs under 2^-61 in F. Above
  argument=50log(2), dropping the remaining log1p term costs under 2^-76;
  rounding the final summed exponent is attenuated by F<2^-25. In the middle
  branch `2alpha>708` and its small argument force r to be tiny, so
  `exp(argument)-rh` cannot incur a dangerous subtraction near one.

The exact-invariant/conversion path for Float64 calls integer/rational code;
its shared converters enter MPFR only in their separate `T===BigFloat`
branches. Constant tuples are immutable, and all new integers/rationals are
local to a call. Root's independent thread/resource gates cover actual task
execution and allocation; this review does not substitute for them.

## Review conclusion and limits

No numerical counterexample or changed-contract failure was found. The new
certificate is conservative under the stated normalized-DD/IEEE arithmetic
assumptions; the separate gap-denominator sensitivity term should be retained
in the repository proof record. The operation-count comment alone is not that
proof. The derivation is not a formal Coq/Gappa proof, and the numerical sweep
is not exhaustive over binary64 bit patterns. The existing native/libm bounds
and root's exact-SHA CI remain their own verification responsibilities.
