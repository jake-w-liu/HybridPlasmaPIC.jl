# Coulomb population-weight normalization repair

Status: confirmed defect and pre-edit plan; implementation has not started.
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

## Tracker and verification

- [x] Reproduce both public failures on the parent and current candidate.
- [x] Identify the population-average owner and downstream sampler dependency.
- [ ] Inspect equivalent weight-sum/mean owners and local-collision callers.
- [ ] Freeze a bounded accumulation/range strategy and resource targets before
  implementation; do not patch the two visible population sizes specially.
- [ ] Verify neighboring counts around significand and exponent boundaries,
  mixed/zero/subnormal/extreme weights, order permutations and all supported
  particle precisions against independent rational means and scattering rates.
- [ ] Reject wrong patches that clamp the mean, use the positive-particle count
  instead of all particles, omit tiny positive weights, or change error/RNG
  behavior. Include no-op, rejected-call and odd-triplet regressions.
- [ ] Reverify compiler/resource/package gates; commit and push this repair
  separately, then review all six exact-SHA CI jobs and complete compiler logs.

The same observations do not establish a defect in every other weighted
reduction. Those owners need their own input/output witnesses before changes.
