# Particle kinematics takeover audit (2026-10-27)

Base: `cd558315956df0debb14b4a2fb7790eaef2821e6`. Work is isolated from the
primary checkout. The existing normwise contract in `DEEP_DEBUG_TRACKER.md`
remains the acceptance standard.

## Current progress

- [x] Reproduce the accepted near-light Lorentz-gap error on Julia 1.12.7/1.13.3.
- [x] Preserve correctness-only and optimized V1/V2 candidates, including their
  allocation/composition failures.
- [x] Preserve V3/V4 and the four old blanket ordinary-bit failures; explicitly
  migrate the retained fixtures to independent physical output conditioning.
- [x] Replace the provisional global DD coefficient with exact residual interval
  propagation; independently review every arithmetic enclosure identity.
- [x] Pass latest independent stage/whole-step witnesses and 3,600 full-expression
  references; pass four-thread caller-context checks on both supported runtimes.
- [x] Finish frozen V5 compiler, allocation, caller-context and actual Metal gates
  on minimum/latest Julia.
- [x] Finish the full retained-input broad regression lane on Julia 1.12.7
  (all 16 affected files, including 100k fixtures).
- [x] Finish that unchanged broad lane on Julia 1.13.3 (all 16 files).
- [ ] Integrate/retest with current main and root's independent radiation repair.
- [ ] Resolve representable output-domain closure after auditing all gamma owners.
- [ ] Define and verify the desired near-light multistep trajectory contract.
  A separate scalar-B probe passes all 32 local actual-stored-input four-ulp
  checks in V5, while 28 sequential-reference points exceed its exploratory
  accumulated 64 ulp/step envelope. Cold exact replay matches that trajectory;
  the local normwise contract does not promise such trajectory-wide accuracy.
- [ ] Complete the separate gathered-midpoint alias repair (root-owned track).

The production certificate is described in
[PARTICLE_PRIMING_INTERVALS.md](PARTICLE_PRIMING_INTERVALS.md). The global
coefficient proposed in earlier pre-edit notes below was superseded by V5; it
is not an accepted error theorem. Neither V1 nor V2 passed the full-step gate,
and the old blanket priming bit-equality gate is not relabeled passed. No
latency acceptance is claimed under shared-machine contention. The source
publication trigger compares rounded native normalized squared speed with 1/2.
A separately tested exact-physical threshold boundary may retain native bits
one ulp from correctly rounded output; its local normwise bound is preserved.
The independent 100k physical fixture classification is not changed to follow
production decisions.

## Pre-edit brief: Lorentz-gap conditioning

- **Mode:** repair of the canonical relativistic Boris and priming guards.
- **Requirement:** accepted finite subluminal inputs must agree with the same
  exact-gap expression within the documented few ulps of the operation scale;
  scalar and certified batch paths must agree. Preserve the native operation
  order on well-conditioned ordinary inputs and the accepted velocity domain.
- **Reproduction:** `nearlight-oracle-v1.jl` in the external takeover evidence
  directory fails all 10 scalar/batch/priming checks on both Julia 1.12.7 and
  1.13.3. Float32 `(prevfloat(0.6f0,2),0.8f0,0)` and Float64
  `(prevfloat(0.6,3),0.8,0)` at `c=1`, `q/m=1`, `dt=2`,
  `Bz=1/sqrt(eps(T))`, zero electric field differ from an independent 4096-bit
  Cayley reference by about 0.208 and 0.168 respectively. The current guards
  admit both. Original registered range suites and actual latest-Julia Metal
  uniform tests pass, so their coverage does not settle this condition.
- **Root cause:** `full_pic.jl` computes `1-β²` from rounded division, products
  and sums. Its guards check finite range without bounding the conditioning of
  that subtraction. Batch certificates also bypass those scalar guards. The
  same rounded Lorentz factor initializes both priming implementations.
- **Contract/invariants:** preserve existing input validation; no clipping,
  relaxed tolerance, hidden device transfer, or new physical formula. Restrict
  native admission to a well-conditioned Lorentz gap; existing exact-gap
  replay owns the remaining accepted evaluations.
- **Canonical owners:** `_boris_rel_in_range`, `_rel_prime_in_range`, their
  batch certificates in `kinematics_range.jl`, and their existing callers in
  `full_pic.jl` and `extrapolated_leapfrog.jl`.
- **Compatibility:** the helpers are private. Scalar ordinary states with
  rounded `β² <= 1/2` retain native bits; larger accepted speeds use existing
  replay. Batch certificates admit only component maxima `<= c/4`, proving
  `β² <= 3/16` without a second particle reduction. Batches outside this
  conservative condition still preserve native results for ordinary particles
  through per-particle admission. Other number types keep native semantics.
- **Strategy:** add the Lorentz-gap condition to both scalar guards; carry the
  already computed `β²` to the Boris guard; prevent both batch certificates
  from bypassing it. The one-half bound keeps the subtraction gap at least
  one-half and its relative sensitivity bounded, rather than selecting a
  fixture-specific band near one. A separate optimization may later reduce
  replay cost with the existing compensated/scaled arithmetic owners.
- **Verification:** retain the original failing logs; rerun the unchanged
  oracle on minimum/latest Julia; add near-boundary, moderate fast-particle,
  sign/permutation/scale, mixed-batch and actual EMPIC step controls; rerun the
  complete registered kinematics and affected EM/priming suites, Aqua and
  concrete/full JET checks. Record allocations and high-speed replay cost as
  observations under contention, with no timing acceptance claim.

## Audit matrix

| Surface | Correctness and robustness | Resources | Status |
| --- | --- | --- | --- |
| CPU nonrelativistic Boris/drift | Registered ordinary/range/public-step tests | No new formula or arithmetic owner | Minimum/latest broad pass |
| CPU relativistic Boris/batches | Independent exact-gap scalar/batch/full-step controls | Both-runtime targeted warm samples zero; timing observation only | Targeted and minimum/latest broad pass |
| Relativistic priming | Residual intervals, independent stage rounding, retained physical classification | Same; high-speed arithmetic cost remains | Targeted gates and both-runtime 100k fixtures pass |
| Metal uniform | Actual functional device, independent oracles and atomic parameter errors | No latency acceptance | 824 assertions pass on each runtime |
| Gathered aliases/atomic validation | Verified midpoint alias defect is separate | Root owns repair | Open |
| Shock coordinates/reflection | Registered finite-drift/reflection and unresolved-step tests | No latency acceptance | Minimum/latest broad pass |

The package-wide audit and simulator roadmap remain open. This record does not
claim that all finite-range physics, device states, or performance gates are
complete.

The frozen V5 campaign completed with all 16 command receipts returning zero.
Its full 16-file regression logs have SHA256
`e5fc8def5afc3ce0660983176b92a227b2f579977b3580a5c240f34001632260`
(Julia 1.12.7) and
`874d435106c83ee278b5e43673d4280aa53bf43023f1ba91002a5321c23cbafa`
(Julia 1.13.3). All 352 frozen Julia source/test files and seven campaign
harnesses were unchanged at completion. These gates used registered
SpectralOperators 0.1.3; integration with root's newer dependency and independent
repairs remains a separate required gate.

## Pre-edit brief: compensated Lorentz-factor optimization

This is a separate candidate; the correctness-only tree remains frozen.

- **Mode:** optimize the verified high-speed replay path without changing the
  physical expression or ordinary native operation order.
- **Baseline:** the frozen repair passes the original exact-gap oracle and
  full-step controls. Its 0.9c Float64 scalar path allocates 5,115,904 bytes per
  256 calls; high-speed 17-particle batches allocate 86,867,688 bytes per 256
  calls. These are allocation observations, not latency acceptance.
- **Existing owner:** radiation supplies pure error-free transformations,
  expansion summation and Float64 double-word arithmetic. Its radiation force,
  error admission, exponential formula and tiny-gap clipping policy belong to
  radiation and will not move. An external prototype using only its arithmetic
  passed 900 independent exact-gap Cayley cases, with 821 fast admissions and
  79 range replays, at the unchanged four-ulp bound.
- **Canonical extraction:** lift the pure arithmetic to
  `src/utils/error_free_arithmetic.jl`, with `_eft_*`, `_dd_*` names and the
  private `_DoubleFloat64` pair type. Migrate radiation call sites mechanically;
  compare every moved definition's normalized AST to its original after the
  exact name mapping. There will be one formula owner and no compatibility
  wrapper or duplicated numerical kernel.
- **Lorentz owner:** `_conditioned_lorentz_factor` returns the historical
  factor for rounded `β² <= 1/2` and other number types. Faster IEEE states
  use a power-of-two scaled `c²-sum(v²)` expansion and double-word sqrt/division.
  Gaps below `2^-900`, nonfinite factors or failed existing range checks retain
  the existing exact-gap replay. No velocity, gap, or result is clipped.
- **Range reasoning:** scaled `c` lies in `[1,2)` and accepted velocity
  components are bounded at that scale. Products cannot overflow. Tiny product
  residuals or scaled components may underflow; their absolute error is a few
  minimum-subnormal Float64 quanta. Requiring the positive gap at least
  `2^-900` makes that error far smaller than the double-word relative scale.
  The sqrt and correction terms down to roughly `2^-1006` then stay normal;
  the factor is below roughly `2^450`. Unresolved cases replay instead of
  substituting a bound. Downstream momentum, rotation and normalization keep
  their native range checks; the ill-conditioned subtraction has been removed
  from the admitted factor. Independent native-formula/oracle and invariant
  tests remain the acceptance evidence for propagated error.
- **Compatibility:** well-conditioned ordinary bits and the input domain stay
  unchanged. Only the already repaired high-speed path gains a nonallocating
  admission. Conservative batch admission remains; no batch can bypass factor
  resolution or scalar range checks.
- **Verification:** exact-factor comparisons, original and widened Boris/
  priming/full-step tests on minimum/latest Julia, all radiation arithmetic and
  relativistic-validation regressions, pure-move AST checks, Aqua/concrete and
  complete compiler review, ordinary/high-speed allocations, shared-machine
  cost observations. Reject wrong patches that repair only scalar Boris,
  round the gap from divided components, or clip a small gap. Performance
  acceptance remains open while the machine is shared.

## Pre-edit brief: certified native-loop specialization

- **Mode:** optimize a redundant branch; preserve optimized V1 unchanged.
- **Evidence:** actual Julia 1.13.3 LLVM for `_empic_kicks!` with
  `Val(true), Val(false)` contains one `_scaled_lorentz_factor` call in V1 and
  none in clean cd55831, for both Float32 and Float64. Raw LLVM is archived.
  Both versions have no vector body; no vectorization-loss cause is claimed.
- **Requirement:** a batch already certified componentwise `<= c/4` must use
  the historical well-conditioned factor directly. Guarded scalar/batch and
  priming calls must retain the compensated factor and unresolved replay.
- **Canonical strategy:** give `_boris_rel_parts` the existing guard mode and
  select the native factor for its certified native caller; use the same
  compile-time guard in both priming loops. `_native_lorentz_factor` owns the
  historical scalar expression, shared with the conditioned owner. No second
  physical formula or new public API is introduced.
- **Compatibility:** only internal specialization changes. Exact accepted
  domain, ordinary bits, compensated high-speed arithmetic and replay policy
  remain unchanged.
- **Verification:** compare native-loop LLVM to parent, all frozen V1 source
  hashes unchanged, unchanged independent witnesses and fast-factor tests,
  minimum/latest affected EM/radiation regression, complete compiler and
  resource reruns. The target remains zero warm allocations; elapsed-time
  samples remain observations under contention.

## Pre-edit brief: accurately rounded priming publication

- **Confirmed failure:** optimized V1/V2 pass individual normwise priming but
  fail all eight full-step comparisons on both runtimes. A strong magnetic
  priming kick leaves a velocity near c; one-ulp publication differences then
  change the next Lorentz factor and final gyro angle appreciably. Original
  failures and both frozen candidates remain retained.
- **Requirement:** publish the correctly rounded primed stage in this
  ill-conditioned slice, so subsequent exact-gap evaluation sees the same
  stored operands as the independent sequential oracle. Keep ordinary native
  bits, no weakened whole-step tolerance and no changed cone/validator policy.
- **Canonical owner:** the existing `_rel_prime_repair` gains a bounded DD
  attempt before its BigFloat replay. Native priming admission requires a
  well-conditioned output, in addition to its existing input/range checks.
  The shared scaled factor owner exposes its DD words without prematurely
  rounding gamma to T. All physical terms follow the same priming expression.
- **Bounded fast slice:** nonzero stored scalar operands are between 2^-100
  and 2^100; a positive bound on all momentum terms is <=2^300. These bounds
  keep force products and normalization inside Float64, with unresolved cases
  using the unchanged BigFloat replay. This is arithmetic admission, not a
  physical speed cap. Structural-zero components are proved from zero source
  factors and do not gain an artificial uncertainty interval.
- **Publication certificate:** propagate a conservative DD absolute error
  radius proportional to eps(Float64)^2 times the sum-of-terms momentum scale
  divided by the final gamma. Require its pre-normalization radius <=c/64.
  For each component, the DD result plus that radius must lie strictly inside
  the rounding cell of one T value; otherwise retain cold replay. The proposed
  2^20 coefficient is deliberately above the primitive/error-count derivation,
  which must be reviewed independently before promotion. Tests alone do not
  establish a uniform error theorem.
- **Verification:** unchanged failed full-step oracles, first-stage exact
  rounding against 4096-bit references, cancellation/range and cell-boundary
  controls, zero-allocation resource slices, all affected runtime/compiler and
  radiation regression gates. The output-domain closure issue remains open;
  this change neither clips velocity nor broadens public validators.

## Pre-edit brief: priming fixture compatibility migration

The V3 implementation passes the unchanged sequential near-light oracle but
fails the old blanket priming equality fixtures for Float32/64 general and
planar priming. The four original failures and the complete V3 source are
frozen outside this tree. This is an explicit compatibility change: priming
whose physical output is faster than c/sqrt(2) must publish the accurately
rounded stage because a normwise one-ulp perturbation there can corrupt the
next Lorentz gap. The native publication threshold remains one half in beta
squared, not a fitted difference threshold.

The existing 100,000-particle fixtures for all three IEEE types and both
priming owners retain every input. An independent 4096-bit exact-gap stage
reference classifies each physical output by squared speed/c squared. The
well-conditioned slice (<=1/2) must retain exact historical bits. The repaired
slice (>1/2) must equal the independently correctly rounded stage. A large
4096-bit comparison margin is checked around the classification boundary; an
ambiguous classification fails the fixture. Classification never calls a
candidate helper and never depends on observed changed bits. Nonrelativistic
fixtures and every end-to-end tolerance remain unchanged. The old blanket
equality gate remains failed; the replacement gate has a new name. The
predeclared zero-allocation targets remain separate and are not relaxed.

## Pre-edit brief: residual-certified priming intervals

V3/V4 remain frozen. Their global `2^20 eps² S/gamma` radius passes the
finite oracle corpus, but the primitive coefficient/normalization and
underflow proof is not sealed. V5 replaces that heuristic with a posteriori
residual certificates. Every double-word operation remains the same center
formula. Exact grow-expansions bound center residuals, with outward-rounded
Float64 radii propagated through addition and products. Division uses
`A-B*z` divided by a certified positive lower bound on `|B|`; sqrt uses
`A-z²` divided by the positive center `z`, a lower bound on `sqrt(A)+z`.
Thus input uncertainty and final-gamma sensitivity are included directly.

The scaled Lorentz-gap owner exposes the same existing square words to
priming; Boris still uses the unchanged sum/sqrt/divide formula. Product
residual certificates allow one minimum-subnormal quantum per FMA product
residual. Exact TwoSum grow steps require finite additions, which the bounded
operand and momentum admission ensures; any nonfinite center/radius or
nonpositive denominator rejects the attempt. The rounding-cell subtraction
also gains its own residual certificate. No DD relative-error theorem or
normalization assumption is required for the certificate: a poor center
merely widens the interval and falls back to existing BigFloat replay.

Verify each interval primitive with independent exact rational endpoints,
including cancellation/zero/underflow; retain actual scalar/batch/full-step
oracles, the explicitly migrated fixture, unchanged resource targets, and
min/latest compiler/radiation/Metal regressions. Arithmetic cost remains an
observation under contention.
