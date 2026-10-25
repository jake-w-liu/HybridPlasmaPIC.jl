# Local Coulomb physical-rate integration

Status: the integrated candidate repairs the confirmed V10 branch discrepancy
on both Julia targets. Numerical boundary and range checks pass; minimum-runtime statistics and thermal checks pass. The broader factor matrix
confirms finite-draw rejection, and ordinary allocation/private-caller gates
fail. Compiler, fixed-cost and promotion gates remain open. This follows the
prescribed-rate sampled-angle and population-mean repairs.

## Verified public witness

`coulomb-weight-mean/local-rate-boundary-followup-v1-source.json` freezes the
candidate and witness. Both Julia 1.12.7 and 1.13.2 report eight passes and two
failures, with every input, exact native-π rational variance, 512-bit physical-π
variance, RNG branch and output-identity result retained in the runtime TOMLs.

Two equal unit particles occupy a unit-volume cell. Their relative speed,
charge, mass, permittivity, weight conversion and Coulomb logarithm are one.
The physical variance is `dt/(2π)`:

- Float64 `dt = nextfloat(2Float64(π))` has variance strictly above one against
  both physical π and the exact rational expression using the native π value.
  The candidate draws a Gaussian angle instead of an isotropic angle. A zero
  Gaussian draw incorrectly preserves the original velocities.
- Float32 `dt = 2Float32(π)` is above the physical-π boundary but equals the
  boundary when π is first rounded to Float32. The candidate draws Gaussian.
  This distinct constant-rounding case must remain explicit in the numerical
  contract; it is not interchangeable with the Float64 arithmetic-loss case.

## Canonical integration requirements

The candidate's `src/particles/coulomb_local.jl` constructs `log_common`,
`log_mu`, `log_base`, then `base = exp(log_base)`. `_coulomb_local_pair!`
passes the rounded base and its logarithm to the older sampler. The published
`src/particles/coulomb_scattering.jl` now owns the shared sampled expression,
original rate factors, sharp branch comparison and complete velocity update.
Simply passing `exp(log_base)` as an original factor would preserve the
observed information loss.

**VERIFIED:** `local-rate-trace-v1-{112,1132}.toml` and six passing trace
assertions per runtime isolate the Float64 loss: `log_base = 0`, `base = 1`,
whereas the original native-factor rational is
`7074237752028441/7074237752028440`, which rounds to `1.0000000000000002`.
Passing that independently rounded factor to the same sampler selects
isotropic scattering. The loss therefore precedes angle sampling.

- [x] Reproduce the public discrepancy on both supported runtime targets.
- [x] Trace the prefactor loss with original-factor counterfactuals; retain the
  physical-π and native-π distinctions and quantify the branch margin.
- [x] Complete a pre-edit contract for numerator/denominator factors, reduced
  mass, cell volume, self-count/triplet factors and unequal-count reuse. Freeze
  the treatment of the mathematical constant π before implementation.
- [x] Migrate local scheduling to the single shared sampler and complete-update
  owner. Remove superseded duplicated arithmetic; preserve pairing, rejection,
  validation, scratch ownership and transactional publication.
- [ ] Reverify physical normalization, sharp boundary neighbors, rate/tangent
  range, all prior cancellation and identity witnesses, unequal-mass/weight
  invariants and statistical relaxation against independent references.
- [ ] Regenerate ordinary full-state differentials and resource/compiler gates;
  preserve intended trajectory changes and verify fixed performance contracts.
- [ ] Recompose the collision-aware stepper/restart candidate with the accepted
  canonical operator, then push and review all exact-SHA CI jobs.

Neither this witness nor the existing V10 passing tests establishes acceptance
of the full local operator. The population-mean helper belongs to the separate
whole-set prescribed-rate API; it must not replace local physical density or
reuse normalization.

## Integration pre-edit brief

- **Mode:** repair and integration of the unmerged CAP-001 physical operator.
- **Baseline:** V10's 20-file permanent source record still matches. Its public
  branch witness fails two of ten checks per runtime. The six-assertion trace
  identifies information lost in `log_base`, before angle sampling. Integration
  starts from accepted population-mean revision `861e562` and preserves V10.
- **Requirement:** retain the original finite input factors through range
  scaling and the sharp variance-one decision. Both physical-π public witnesses
  must pass. Keep the current empirical partner-density model and canonical
  angle/update law; this is not a change to pairing physics.
- **Contract:** for positive-weight pairs, use
  `dt*number_scale*logLambda*(qa*qb)^2*(ma+mb)^2*wmax*count /
  (8π*epsilon0^2*ma^2*mb^2*prod(dx)*reuse*triplet*speed^3)`, where
  `speed=max(hypot(g),u_floor)` is the stored relative-speed value. Self channels
  use `count=N-1,reuse=1`; triplet edges use `triplet=2`, other edges one.
  Cross channels use `count=max(Na,Nb)` and the exact integer reuse count.
  Zero-weight particles remain excluded from local density, unlike the
  prescribed-rate whole-population mean.
- **π policy:** the physical prefactor uses mathematical π. Ordinary floating
  evaluation approximates it; near the branch boundary an exact rational of
  the other stored factors is compared against certified π bounds. The
  existing native azimuth calculation and draw order remain unchanged. This
  intentionally distinguishes the Float32 rounded-π equality from the
  physical boundary; no assertion is weakened to preserve that equality.
- **Canonical owner:** `coulomb_local.jl` owns physical prefactors, count/reuse
  scheduling, workspace and transactional publication. The accepted
  `coulomb_scattering.jl` owns the sampled angle, increment and complete
  velocity expression. Add narrow physical-rate dispatch to its existing
  variance/boundary interface; remove V10's duplicated mass-update and complete
  expression functions. Migrate private tests instead of retaining the old
  sampler as a compatibility shim.
- **Strategy:** construct a binary-scaled channel rate from normalized positive
  factors; retain charges, masses, units, widths and counts separately. A stable
  `(1+min(m)/max(m))^2/min(m)^2` approximation avoids materializing the reduced
  mass. Renormalize the positive mantissa after each product/division. Preserve
  original masses for the exact boundary comparison. Count and reuse remain
  integers until their approximate fast evaluation, with exact values retained
  for the boundary. Define a forward-error guard from the bounded operation
  count, including spatial dimension, before selecting the cold comparison.
- **Boundary certificate:** use exact alternating rational bounds for
  `π=16atan(1/5)-4atan(1/239)`. Each arctangent's remainder lies between zero and
  its next signed term. Refinement stops only when the rational target is
  outside the enclosing interval; π is irrational, so a finite rational target
  cannot equal the limit. This avoids changing or depending on global MPFR
  precision/exponent settings. No floating equality or silent iteration cap
  substitutes for a certificate.
- **Compatibility:** the new local API is unmerged; keep its documented
  validation, species order, scratch ownership and failure semantics. Historical
  trajectories can change with corrected rate assembly and the already accepted
  increment law. Preserve every old differential output and explain changes.
  Released prescribed-rate behavior must remain bitwise unchanged from `861e562`.
- **Failure and resources:** validate before sampling; stage every population
  until all channels finish. Rejection retains public velocities, while RNG and
  scratch may advance as documented. Ordinary warmed local calls with workspace
  must retain zero allocation. Inventory exceptional allocations, complete-step
  costs and retained storage; keep V10's fixed 1.5 latency ceilings unchanged.
- **Verification:** run the unchanged public boundary witness, neighboring
  rate/count/mass/volume cases, independent physical rates and all V10 permanent
  numerical/statistical/range/normalization/resource suites on Julia 1.12.7 and
  1.13.2. Rebuild full ordinary outputs, concrete JET, complete unsuppressed
  package trees, Aqua, formatter and the public example. Recompose CAP-010 and
  require exact-SHA CI after separate promotion.
- **Wrong implementations to reject:** passing `exp(log_base)` as an original
  factor; using rounded reduced mass or rounded count/reuse in the cold branch;
  replacing mathematical π with native-π equality; counting zero-weight local
  particles; retaining duplicated V10 angle/update code; publishing one species
  before a later channel fails.

The independent π prototype V2 passes 42 assertions per runtime, including
directed-MPFR endpoints at 24–1,024 bits, an exact tangent-addition identity
check and a compact-exponent check where floating π overflows. Its comparator
uses 4–219 series iterations in those profiles; allocation is recorded, not
accepted as a performance claim. V1's two failures per runtime are retained:
the test converted 512/1,024-bit MPFR endpoints using an active-precision
conversion, losing their direction. A stored-bit decomposition trace identifies
that oracle defect; V2 leaves the comparator unchanged.
Evidence: `2026-10-26/local-coulomb-rate/pi-boundary-prototype-v2-commands.json`
and `pi-prototype-oracle-review-v1.json` in the external audit directory.

## Integration verification

The candidate starts at `861e562`, preserving the unmerged V10 snapshot. Its
physical factors and integer count/reuse values reach the canonical sampler;
the duplicate V10 angle and velocity arithmetic is removed.

| Gate | Observed result on Julia 1.12.7 and 1.13.2 |
| --- | --- |
| Unchanged public branch witness | 10 passes each |
| Initial public physical/conservation/rollback suite | 1,260 passes each |
| New physical-π, integer-count and original-mass boundary matrix | 3,728 passes each; three deliberately rounded-factor alternatives rejected |
| Prior complete-angle, update, cancellation and normalization range suites | 147,009 passes each, with original tolerances retained |
| Whole-tree formatter | JuliaFormatter 1.0.62; 24 AST/idempotence checks pass |
| Statistical, thermal, ordinary allocation, neighboring and prescribed-rate differential suites | Running against frozen integration source |
| Concrete and full package compiler checks | Running; complete report comparison required |
| Fixed latency, final acceptance and exact-SHA CI | Pending |

The boundary matrix also tests finite maximum-magnitude Gaussian draws. These
cases pass; the suspected tangent overflow did not reproduce in this matrix.
That observation is bounded evidence, not a proof over every admissible input.

Evidence: `2026-10-26/local-coulomb-rate/{focused-v1-commands.json,
numerics-v1-commands.json,integration-v1-source.json,integration-v1-commands.json}`.
The frozen archive retains all production files and migrated private tests.
Changing a private sampler call does not remove the associated range assertions.

## Retained integration failures

The first complete minimum-runtime integration suite has 1,782 passes, six
allocation failures and one private-call error. All 192 empirical Landau and
14 thermal assertions pass. The error identifies an unmigrated nine-argument
private `_coulomb_mass_update` test; the canonical owner now takes the complete
finite increment. That test must be migrated without weakening its oracle.

The six allocation failures are the Float32, four-species, 4,096-particle
profiles for every dimension/collection form: the first measured call uses
33,448–34,712 bytes, followed by two zero-byte calls. Same-state traces on both
runtimes identify original-factor rational conversion and the π certificate.
The copied RNG used for stack attribution also grows its integer buffer;
profiled inclusive byte totals are not interchangeable with the original
`@allocated` measurements. A Float64 certificate using original factors is
being investigated before the exact rational comparison; no gate is relaxed.

The broader independent factor matrix passes 8,980 checks and fails 23 on each
runtime. The forward-error and regime comparisons pass, but finite maximum
Gaussian draws throw when the approximate variance lies two ulps above one
and exact variance is below one. The same defect is reproduced in six public
prescribed-rate cases per runtime and is being repaired separately in
`HybridPlasmaPIC-coulomb-finite-draw-20261026`. V1's missing diagnostic write
and V2's tuple/TOML serialization error are retained; V3 preserves every
original input bit and failure without changing the candidate or assertions.

Minimum-runtime neighboring regressions pass 49,726 assertions. Other frozen
integration commands remain active. Evidence: `factor-guard-v3-*`,
`allocation-trace-v1-*` and `integration-v1-commands.json`.
