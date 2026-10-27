# Coulomb population-weight normalization repair

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: V2 passes final local verification on Julia 1.12.7 and 1.13.2;
promotion and exact-SHA CI review remain. V1 is retained unaccepted after the
coarse-BigFloat follow-up.
This is a separate repair after the sampled-angle/update slice.

## Evidence and owner

**VERIFIED:** `coulomb-angle/weight-mean-followup-v1.jl` records four passes
and eight failures on each of Julia 1.12.7 and 1.13.2, testing both frozen
`a18f09d` and the sampled-angle candidate. Two Float16 populations incorrectly
select isotropic scattering when an exact rational weight/count reference
requires the Gaussian branch. A zero Gaussian draw should preserve velocities;
both implementations instead alter them.

- With 4,096 uniform weights, the sequential Float16 total stops increasing
  before the full count. The computed mean is too small.
- With 65,536 particles and only two positive unit weights, converting the
  particle count to Float16 produces infinity, so the computed mean is zero.

`collide_coulomb!` in `src/particles/collisions.jl` owns normalized-weight
validation, `Wtot` accumulation and `wbar = Wtot / T(N)`. The sampled-angle
owner consumes that mean; correcting its downstream arithmetic cannot recover
an already incorrect population average. Local multispecies cell-density
normalization is a separate caller/contract that must be inspected before
sharing an implementation.

## Repair contract

- Compute the population's normalized mean without loss from a narrow
  accumulator or converting a large count to the particle type too early.
  Preserve original count/rate factors through any unrepresentable mean.
- Retain finite/nonnegative weight validation, positive normalized-weight
  rejection, all-zero no-op behavior, RNG order and transactional publication.
- Keep an allocation-free ordinary CPU path. Measure large-population and
  exceptional costs; do not allocate an exact rational per particle by default.
- Treat changed collision rates for affected populations as the intended fix.
  Retain all ordinary profiles and record changes instead of preserving an
  incorrect rate for historical floating-point equality.

## Pre-edit brief

- **Mode:** repair of the whole-population normalization owner.
- **Baseline:** `coulomb-weight-mean/before-v1-commands.json` records two passes
  and four failures on each runtime at `4f26e1d`, reproducing both public
  witnesses above. Their exact rational variance requires Gaussian sampling.
- **Root cause:** the owner accumulates normalized weights in `T`, then divides
  by `T(N)`. Float16 loses additions at its significand boundary and cannot
  represent the larger count. The downstream sampler cannot recover this loss.
- **Canonical owner:** a private normalized-mean helper in `collisions.jl`,
  called once by the existing public whole-set collision path. The pair and
  shared sampler must accept its accumulator type. The separate cell-local
  candidate uses physical partner density, count and reuse factors, not this
  population mean; it must not substitute this normalization.
- **Invariants:** keep normalization in the particle type (including the
  existing rejection of a positive weight that normalizes to zero). Include
  zero weights in the population count. Validate before scratch publication or
  RNG consumption. Keep pair order, rejection draws and all-zero no-op rules.
- **Strategy:** sum the native normalized weights with Neumaier compensation
  in `promote_type(T, Float64)`, divide in that accumulator type, and retain the
  mean in that type through variance construction. For hardware floats this
  supplies normal finite range for every `Int` population count and its mean:
  at least one normalized weight equals one, so the mean is at least `1/N`.
  Cast the final weight quotient/scaled binary variance to the particle type;
  retain the wide mean in the existing exact regime-boundary comparison.
  BigFloat retains its precision and exponent policy; this slice makes no
  new arbitrary-MPFR-range or exact-summation claim.
- **Compatibility:** no public API or persistence change. Private sampler
  signatures accept a distinct mean type. Changes to rates caused by more
  accurate means are intentional; seed determinism within the new version,
  draw ordering for the same regime, and transactional failure stay required.
- **Acceptance:** independent rational means of native normalized inputs;
  uniform, sparse, mixed, subnormal and scaled populations; neighboring and
  odd counts; public branch/kick checks; unchanged rejection/no-op behavior;
  zero warmed allocations for ordinary Float32/64 calls with supplied scratch.
  Preserve every differential profile and raw compiler report. BigFloat
  ordinary precision is a regression surface, not an inferred range proof.
- **Wrong patches to reject:** widening only the divisor leaves stalled sums;
  counting only positive weights changes the physical rate; rounding the mean
  back to Float16 loses sparse populations; dropping small addends loses a
  collectively representable contribution. Tests must distinguish each.
- **Verification commands:** both installed Julia versions use the retained
  `env-112`/`env-1132` manifests. Run focused tests, the existing collision/angle/
  sampled-update suites in a `Pkg.test` sandbox, ordinary and exceptional
  allocation profiles, concrete JET signatures and full unsuppressed package
  reports, then the repository JuliaFormatter gate and exact-SHA GitHub CI.

## Tracker and verification

- [x] Reproduce both public failures on the parent and current candidate.
- [x] Identify the population-average owner and downstream sampler dependency.
- [x] Inspect equivalent weight-sum/mean owners and local-collision callers.
- [x] Freeze a bounded accumulation/range strategy and resource targets before
  implementation; do not patch the two visible population sizes specially.
- [x] Verify neighboring counts around significand and exponent boundaries,
  mixed/zero/subnormal/extreme weights, order permutations and all supported
  particle precisions against independent rational means and scattering rates.
- [x] Reject wrong patches that clamp the mean, use the positive-particle count
  instead of all particles, omit tiny positive weights, or change error/RNG
  behavior. Include no-op, rejected-call and odd-triplet regressions.
- [x] Reverify the collision sandbox, compiler, resource and formatting gates.
- [ ] Commit and push this repair separately, then review all six exact-SHA CI
  jobs and complete compiler logs. Fixed latency remains open optimization work.

The same observations do not establish a defect in every other weighted
reduction. Those owners need their own input/output witnesses before changes.

## Verification in progress

The first candidate passes all 3,152 new permanent regression assertions on
both runtimes, plus six assertions in the unchanged original public witness.
The permanent checks cover independent exact means, permutations, large and
sparse populations, native subnormal weights, collective small contributions,
public sampled kicks and validation-before-mutation. A lazy 33,554,433-entry
weight vector verifies a positive mean that would round to zero in Float16
without allocating an artificial particle population of that size.

`focused-v1-source.json` and `focused-v1-commands.json` preserve the tested
source, commands and log hashes. The completed V1 matrix below did not settle
the coarse-precision follow-up and does not establish acceptance.

## Expanded precision repair brief

The first candidate is **not accepted**. Its complete standard matrix passes
on both runtimes (58,107 sandbox assertions, 12,096 ordinary checks, 44 concrete
signatures, 64 mean-resource/wrong-patch checks and 720 exceptional checks), but
`coarse-bigfloat-followup-v1` independently records 19 passes and nine failures
per runtime. At BigFloat precisions 4, 8 and 16, a uniform population of
`4 * 2^precision` has computed mean 0.5 instead of one, chooses isotropic
sampling and changes velocities for an otherwise zero Gaussian draw.
The two compensated accumulators both stop accepting unit additions.

The same canonical owner must handle this before promotion. Extend its brief:

- Preserve native normalized-weight validation and use a BigFloat accumulator
  with explicit local precision `p + bitlength(N) + 2`. For nonnegative inputs
  this reserves count-dependent guard bits for accumulation; the final mean
  is rounded at the caller's active precision using an exact count operand.
- Use fresh local MPFR output objects for in-place addition/division, following
  Julia's own `sum(::AbstractArray{BigFloat})` calling convention. Do not call
  `setprecision` or change process-wide rounding/exponent settings inside the
  production operation. Reuse the normalized-weight validation owner.
- Keep hardware-float behavior and its zero-allocation gate. This is not an
  arbitrary compact-MPFR exponent-range repair: existing finite-total rejection
  remains until that separate range contract is implemented.
- Prototype the accumulator against the unchanged public failing witness,
  rational mixed-weight means, precisions around hardware boundaries and all
  supported rounding modes. Verify settings restoration and allocations.
- Add permanent low-precision and wide-mean/sampler regressions, repeat relevant
  regression/compiler/resource checks after the implementation change, and
  retain the entire unaccepted V1 freeze and logs.

`final-v1-source.tar.gz`, `final-v1-unaccepted-review.json` and the coarse
follow-up source/command receipts preserve the earlier candidate.

The V2 prototype passes the unchanged public precision witness (28 assertions),
2,161 independent mixed-mean/rounding/settings checks and 11 compact-exponent
integer-divisor checks on both runtimes. The implementation now uses one shared
normalized-weight validator, a hardware compensated sum and a BigFloat local
accumulator with count-dependent guard bits. Division by `BigInt(N)` uses
Julia's native `mpfr_div_z` path so the count remains exact even when a floating
count would overflow. No per-particle rational is constructed. All new checks,
including the 390 wide-mean/sampler checks, are registered permanently.

### Guard-precision bound

For the exact values of the already normalized native inputs, let `p` be the
active precision, `N ≥ 1`, `k = bitlength(N)`, and `q = p + k + 2`.
The [MPFR rounding contract](https://www.mpfr.org/mpfr-current/mpfr.html#Rounding)
rounds each addition to its destination precision. In the absence of exponent
clipping, each positive operation has relative error at most `εq = 2^(1-q)`
in any of the five tested rounding modes. Set `ε = 2^(1-p)`; then
`N*εq < ε/4`. For a nonnegative sum, the accumulated relative error is bounded
by `γ = N*εq/(1-N*εq)`. Division by the exact integer count and rounding to
`p` bits therefore gives relative error at most
`γ + ε*(1+γ) < 5ε/(4-ε) ≤ 5ε/3 < 4ε` for `p ≥ 1`.
This establishes the permanent test's `4ε` bound for this domain; it does not
claim exact summation or include error from the preceding native normalization.

`guard-precision-bound-v1.json` records 258,048 exact-rational checks of these
inequalities for precisions 1–1,024 and every 64-bit positive `Int` count
bitlength, evaluated at each bitlength's largest count. Compact MPFR exponent
limits remain a separate range contract: a total may overflow, or a mean may
be unrepresentable, even with sufficient significand precision.

### Final local verification

`final-v2-source.json` and `final-v2-commands.json` retain the exact source,
commands and log hashes. All twelve commands exited successfully, six per
runtime. Each runtime records:

- 60,697 assertions in the collision sandbox, including 5,742 population-mean
  regressions and Aqua checks;
- 12,096 ordinary checks over 576 profiles, with baseline and candidate both
  allocating zero warmed bytes when supplied scratch; three varied-weight
  Float32 profiles have intentional rate-dependent trajectory changes, and
  every full state is retained;
- 48 concrete compiler assertions, 64 mean-resource/wrong-patch assertions,
  and 720 exceptional-resource assertions across 72 profiles;
- complete package compiler reports: 28/26 trees unchanged, plus one private
  pair-entry signature accepting the wider mean. Its entire `hypot` callee
  chain and existing abstract `Number`/bitcast diagnostic are unchanged.
  `compiler-v2/explicit-signature-review-v1.json` records the full comparison;
- the unchanged hardware exceptional allocation maxima; BigFloat's measured
  maximum is 2,393,680 bytes, below V1's 2,467,552. This is allocation evidence,
  not a latency claim.

JuliaFormatter 1.0.62's AST-preservation/idempotence gate passes, and all frozen
production/test files still match their recorded hashes after verification.
Final source/caller review finds one normalization owner, no public API change,
no test-specific production branch and unchanged validation-before-publication.
The coarse-precision public path and independent rational-mean path both pass;
compact exponent-range completion and local physical-rate integration remain
explicitly outside this repaired population-mean scope.
