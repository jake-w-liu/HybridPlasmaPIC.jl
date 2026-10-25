# Local Coulomb physical-rate integration

Status: confirmed branch discrepancy in the unmerged V10 candidate;
implementation and acceptance remain open. This follows the prescribed-rate
sampled-angle and population-mean repairs.

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
- [ ] Complete a pre-edit contract for numerator/denominator factors, reduced
  mass, cell volume, self-count/triplet factors and unequal-count reuse. Freeze
  the treatment of the mathematical constant π before implementation.
- [ ] Migrate local scheduling to the single shared sampler and complete-update
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
