# Bounded BigInt whistler counters

`run_whistler` converts a positive `BigInt` step count to `Int` when it fits
exactly and the timestep is a built-in Float16, Float32 or Float64. It then uses
the existing fixed-counter loop. Every amplitude still undergoes every ordered
multiplication; no multiplier powers or additional arrays are used.

Generic timesteps retain their original evaluation order and call counts.
Zero, larger BigInts and custom counters retain the original path. The helper
never clamps a count. Boundary tests inspect counts around `typemax(Int)` and
larger values without attempting astronomically long simulations. Validation,
frequency construction, seeded draws, band filtering and energy reductions keep
their order. This changes no coupled solver or persistence format.

## Measured owner and acceptance

Frozen parent `f93e17d` allocates an extra 19,056 bytes for a 100-step BigInt
packet on both local Julia versions. Full traces identify 397 additional
BigInts in range construction and iteration. The canonical driver can avoid
that work by using an exactly representable fixed-width counter.

The separate [CI investigation](WHISTLER_ALLOCATION_CI_DIAGNOSIS.md) retains
both failed package runs and all four focused Ubuntu artifacts. Their original
tests pass, but retained repeated allocation samples vary substantially on both
source revisions. Native usable-size accounting is an untested explanation for
that variation. This optimization addresses the measured iterator cost; it
does not claim to establish the historical failure's precise cause.

The original 4,096-byte gate is unchanged. A new permanent gate requires every
declared positive bounded-BigInt call to stay within 512 bytes of its matching
Int call. All three samples are tested. Independent histories preserve full
output bits, including signed zero and NaNs, rather than only energy or norms.

Evidence is under `takeover/whistler-bigint-count` in the dated audit archive.
The pre-edit brief fixes the compatibility and resource targets. The original
formatter comparison failure is retained; its read-only successor checks all
three Julia files after normalizing line-location metadata and passes nine
syntax/format/idempotence assertions. The original packet test is byte-identical
to the parent.

## Progress

- [x] Verify the allocation owner and preserve the parent and failed CI logs.
- [x] Implement exact bounded-count normalization at the canonical driver.
- [x] Pass 12,956 focused assertions on Julia 1.12.7 and 1.13.3.
- [x] Prove the complete source outside the intended count change unchanged,
  including the generic loop, and reject a changed-energy control.
- [x] Demonstrate the new cost gate fails on the frozen parent on both runtimes.
- [x] Compare all original history fixtures plus BigInt prefixes 0–100,
  exceptional frozen endpoints and generic BigFloat caller contexts.
- [x] Pass all 148 complete-call resource profiles per runtime: the 76 original
  fixtures plus 72 counter fixtures, with three alternating paired samples,
  exact output hashes, no compilation in measured windows and no byte increase.
- [x] Pass 78,801 covered assertions per runtime and all concrete/Aqua/caller
  gates. Complete 29/27-report compiler forests, renders and all 252/203 actual
  frames match the parent under repository-root substitution only.
- [x] Rerun the resource supplement with actual argument, counter and timestep
  types: all 296 typed profiles and 888 paired samples across both runtimes pass.
- [x] Load the original and new test extensions before the complete quality
  campaign: all four matched commands pass, with every compiler tree, render
  and actual frame field equal under repository-root substitution only.
- [x] Complete independent final source/evidence review and docs; preserve all
  generated coverage files and verify the primary worktree files unchanged.
- [x] Push this optimization separately to main as `cd8b94d`.
- [ ] Complete all six exact-SHA jobs; see the [current audit tracker](DEEP_DEBUG_TRACKER.md).

The old resource fixture's 95% reduction target compared against the pre-hoist
implementation. This follow-up retains its complete fixture inventory and
checks nonregression against the already optimized `f93e17d` parent. It adds the
bounded-counter target above. No latency, whole-simulator optimization or
completeness claim follows from these allocation measurements.
