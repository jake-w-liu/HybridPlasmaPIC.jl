# Finite Gaussian Coulomb draws at the regime boundary

## Pre-edit brief

- **Mode:** repair the released canonical sampled-angle path.
- **Requirement:** a finite Gaussian draw and finite representable physical
  update must not be rejected because the approximate tangent overflows.
  Nonfinite draws must still reject transactionally. Ordinary full state, RNG
  order and warmed allocations must remain unchanged.
- **Baseline:** the independent local physical-factor matrix has 8,980 passes
  and 23 finite-draw failures per runtime. The released prescribed-rate API
  independently reproduces six cases (three each Float32/Float64) per runtime;
  its 19 trace assertions confirm rejection and unchanged public velocities.
  The equal-weight search did not reproduce the issue; its failed search gate
  is retained. The unequal-weight search freezes all original input bits.
- **VERIFIED cause:** `coulomb_scattering.jl::_coulomb_sample_increment` selects
  Gaussian scattering using exact original factors, while its approximate
  variance is two ulps above one. `sqrt(var)*floatmax(T)` overflows; the current
  unconditional finite-tangent check rejects it. Exact variance is below one,
  and the physical rotation remains bounded. The trace retains exact rational
  variance, approximate variance, draw, exception and rollback.
- **Canonical owner:** the shared sampler already stores a signed binary
  fraction/exponent tangent and routes difficult cases to complete physical
  velocity evaluation. Both public collision APIs use this owner.
- **Contract:** validate the random draw itself. A finite draw with an infinite
  approximate tangent uses the existing scaled tangent and complete-update
  representation. The reciprocal-tangent direction formula remains bounded;
  complete evaluation retains tiny physical transverse terms. Do not clamp
  variance, synthesize a draw, bypass rejection, or add a second sampler.
- **Compatibility:** previously rejected valid finite draws may now complete.
  Invalid random samples retain ArgumentError and public rollback. No public
  API or checkpoint schema changes. The local feature remains independently
  unaccepted; this repair starts from main861e562 in a separate worktree.
- **Verification:** retain frozen main witness, add original and neighboring
  signs/magnitudes plus nonfinite cases, compare complete physical updates to
  a high-precision independent oracle, rerun prescribed-rate regressions and
  ordinary full-state/allocation differentials, concrete/full compiler trees,
  Aqua, formatter and required exact-SHA CI after separate push. Recompose
  the frozen local candidate later and rerun its 23 failing rows.
- **Wrong patches to reject:** accepting nonfinite draws; clipping variance or
  tangent to an arbitrary finite value; returning the unchanged state or only
  the rounded infinite-tangent direction; consuming additional RNG draws.

Evidence directory: `2026-10-26/local-coulomb-rate`,
`factor-guard-v3-{112,1132}.toml` and
`prescribed-extreme-witness-v2-{112,1132}.toml`. Source and command records
remain alongside them. The implementation is present; acceptance is pending.

## Focused verification

The permanent regression reproduces rejection on unchanged main861e562 with
19 passes and one error on each runtime. After the canonical repair, 1,284
finite-update assertions and 54 invalid-draw/rollback assertions pass on both
Julia 1.12.7 and 1.13.2. The oracle uses original factors and 4,096-bit physical
expressions; it retains small transverse updates that a rounded backscatter
direction would discard. Both draw signs, neighboring coefficients, zero draws
and near-maximum magnitudes retain the expected RNG sequence.

The first latest-runtime fixture run had two metadata comparison failures:
positions were left uninitialized. `ParticleSet` explicitly allocates position
storage with `undef`; the corrected fixture initializes finite coordinates
before capturing state. Physical assertions and tolerances are unchanged. The
original failed source/logs and the exact limits of that attribution remain in
`fixture-review-v1.json`.

JuliaFormatter 1.0.62 passes six AST/idempotence checks. The frozen final
regression, resource, ordinary-output and complete-compiler campaign passes
on both runtimes.
Evidence: `2026-10-26/coulomb-finite-draw/final-v1-source.json` and
`final-v1-commands.json`. No latency measurement is claimed. Separate promotion and exact-SHA CI remain.

## Final local acceptance

All twelve final commands exit zero, with unchanged frozen production/test
hashes. Each runtime passes 62,035 sandbox assertions across 50 testsets,
16,704 ordinary differential assertions, 48 concrete compiler assertions,
576 new resource/replay assertions and 720 prior exceptional resource checks.
The new resource harness also repeats all 1,338 focused regression assertions.

All 576 ordinary profiles retain complete velocity/scratch/RNG outputs and
zero warmed allocation relative to main861e562. The 36 new exceptional profiles
retain deterministic outputs and stable measured allocation across identical
inputs; the largest is 2,858,496 bytes for 64 pairs. These cold paths require
exact-boundary/complete-expression work; they are not ordinary zero-allocation
or latency claims. The 72 prior exceptional profiles retain their gates.

Complete package compiler comparison preserves every report frame, signature
and diagnostic: all 29 minimum-runtime and 27 latest-runtime trees are
identical after checked unchanged-source line maps. Both complete normalized
trees remain in `compiler-v1/complete-review.json`; no count-only acceptance
or diagnostic suppression is used.

The final source/caller audit finds one sampler and no clipping, new fallback
API, extra RNG draw or publication-path change. The independent physical
formula tests cover tiny transverse effects that a rounded infinite-tangent
rotation would erase; nonfinite draws still reject transactionally. The local
physical-rate feature and its remaining allocation/validation gates are
separate. Evidence: `final-v1-verified-results.json`, complete raw logs and
`final-v1-source.json` in the directory above.
