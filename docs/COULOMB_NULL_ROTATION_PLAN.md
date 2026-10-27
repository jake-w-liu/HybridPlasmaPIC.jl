# Coulomb exact identity rotation repair

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: canonical repair is pushed as `a18f09d` after local numerical, RNG,
resource and compiler checks on Julia 1.12.7 and 1.13.2. Exact-SHA CI and full
compiler-log review remain.
The separate nonzero-small-angle defect below and broader CAP-001 remain open.

## Pre-edit brief

- Mode: repair the prescribed-rate whole-set Coulomb pair update.
- Requirement: an exactly sampled identity rotation (`cosθ == 1`, `sinθ == 0`)
  leaves all original velocity bits unchanged, including signed zeros and tiny
  components. Preserve angle/azimuth/rejection draw order and error behavior.
- Baseline: main `839c109`; the retained `legacy-null-v1-1132.log` records eight
  failed identities in parent/local candidates. Reproduce against this main
  before production edits.
- Root cause: `_coulomb_pair_scatter!` reconstructs accepted velocities as
  midpoint plus/minus half the rotated relative velocity. That algebraically
  equivalent reconstruction can round away an original tiny component even
  when the sampled rotation is identity.
- Canonical owner: `src/particles/collisions.jl` pair update, called by the
  public whole-set `collide_coulomb!` for pairs and odd-population triplets.
  Existing transactional velocity scratch and input validation remain owners.
- Compatibility: preserve every nonidentity path, draw and error. No API/schema
  change. Exact identity outputs intentionally replace the incorrect rounded
  reconstruction. General nonzero small-angle accuracy remains a separate item.
- Strategy: recognize identity after the existing validation and rejection
  draws, before writing velocities. This is the exact algebraic identity, not
  a tolerance cutoff, fallback, or failure suppression. Preserve zero-weight,
  equal-velocity and public no-op validation behavior.
- Verification: retain before/after tests on Julia 1.12.7 and 1.13.2; vary
  dimension, Float16/32/64/BigFloat, pair/triplet/even counts, weights, signs,
  axes and tiny/signed-zero components. Compare RNG events, unchanged metadata,
  scratch transactionality, invalid random samples and nonidentity rotations.
  Differentially compare ordinary full outputs/RNG/allocation against frozen
  main, run collision neighbors, Aqua, concrete/full compiler checks and
  repository formatting. Require unchanged ordinary warmed allocations.
- Wrong patches: returning before azimuth/rejection draws changes the RNG;
  testing only `sinθ == 0` hides a 180-degree swap; an epsilon cutoff erases
  valid small rotations. Explicit controls cover each distinction.
- Integration: independent main fix after applicable gates; do not edit the
  frozen CAP-001 candidate. Its later integration must compose this repair and
  reverify the changed legacy compatibility baseline explicitly.

## Tracker

- [x] Reproduce current-main identity failure on both Julia runtimes.
- [x] Implement the canonical identity repair with permanent regressions.
- [x] Verify targeted, neighboring, error, RNG, resource and compiler contracts.
- [x] Independent final source/evidence review for this bounded repair.
- [x] Individual main push as `a18f09d`.
- [ ] All six exact-SHA CI jobs with full compiler review.
- [ ] Separate nonzero-small-angle repair and remaining CAP-001 acceptance.

The first current-main test run exposed a fixture dispatch ambiguity between
its generic `randn` extension and Random's hardware-float methods. Both failed
runs and their exact fixture are retained as `before-v1-*`; those method errors
are not counted as production failures. Concrete test RNG methods remove the
ambiguity before the baseline replay. No production source changed.

The corrected `before-v2` fixture has 6,912 passes and 3,456 failed velocity/scratch
identity checks on each runtime, with no errors. All metadata and draw-order
checks pass. The exact null repair is placed after the original finite-result
checks and acceptance draws, preserving those boundaries.

## Final local verification

Evidence is retained in `HybridPlasmaPIC-audit-evidence/2026-10-23/coulomb-null/`.
Both runtimes pass the following unchanged gates:

| Check | Result per runtime | Retained evidence |
| --- | --- | --- |
| Actual package test sandbox | 10,404 new + 512 existing collision + nine Aqua assertions | `final-v1-sandbox-*.log` |
| Ordinary differential replay | 8,064 assertions across 576 profiles; exact outputs, scratch, RNG and warmed bytes match frozen main | `final-v1-ordinary-*.log` |
| Concrete compiler checks | Eight signatures pass | `final-v1-concrete-*.log` |
| Executed wrong-patch controls | Early return, sine-only identity and epsilon cutoff are all rejected | `final-v2-mutants-*.toml` |
| Complete package compiler trees | 29/27 existing reports; every normalized tree matches frozen main with checked unchanged-line mappings | `compiler-v1/*-compiler-review.json` |

Whole-tree JuliaFormatter 1.0.62 verification passes six checks, including AST
identity and two clean formatting passes. Final source freezes retain 104/105
hashes. The first mutation-control launcher failed on Julia world age; the
corrected launcher uses `invokelatest`, and both the failure and repaired run are
retained. This was a harness error, not a production failure or a passing gate.
The final review rechecks source and command/log hashes; no test or compiler
finding is suppressed.

## Independent neighboring finding

VERIFIED: `legacy-variance-underflow-v1.jl` fails four of eight assertions on
each runtime, against both frozen main and this null candidate. For Float32
`g=1e12, gcoeff=1e-12`, and Float64 `g=1e100, gcoeff=1e-100`, with dt=1 and
a unit Gaussian draw, the variance division rounds to zero before its square
root. The independent BigFloat formula gives finite transverse kicks near
-1e-12 and -1e-100, respectively; both implementations return zero. This is
separate from a genuinely zero Gaussian draw. Logs, commands and source hashes
are retained in `coulomb-null/variance-underflow-review-v1.json`.

The next small-angle repair must carry the complete angle/kick scale before
rounding and update by a stable velocity increment. It must cover variance
underflow, zero and tiny Gaussian draws, oblique cancellation, overflow restored
by final scaling and ordinary seeded compatibility. The existing local
CAP-001 scatter/increment owner is an integration dependency to inspect, not
a reason to copy its arithmetic into a second sampler.

Fixed latency acceptance remains open under the shared 1.5 screening ceiling;
the current concurrent process inventory does not provide an exclusive window.
This repair makes no latency or whole-Coulomb-correctness claim.
