# Local physical-rate allocation repair

The held V3 candidate passes the numerical, allocation, state/RNG and compiler
checks listed below on Julia 1.12.7 and 1.13.2, including neighboring tests and
executable examples. The final audit checks 9,017 source/evidence invariants.
Landing verification on `main` after `3a9c85d` (2026-10-26): all seven new
test files, ten neighboring collision/sorting/normalization files, Aqua, JET
and the example pass on Julia 1.12.7 and 1.13.2 with registered SpectralOperators
0.1.3. JET reports (24 / 21) match the parent apart from one line shift in
`sorting.jl`. The new files add about five minutes of test time per runtime.

Measured cost on one Apple M5 core (16³ cells, 16 particles per cell and
species, equal weights, warmed, shared machine): the local operator costs
109 / 206 ns per pair for one Float32 / Float64 species versus 73 / 121 ns per
pair for the homogeneous operator on the same particles, and stays
allocation-free (the homogeneous operator allocates per call). Most of the
time is spent in the shared compensated velocity update; per-pair coefficient
evaluation and binning are the scheduling overhead. Optimization is tracked
as follow-up work; exact-commit CI is pending.

## Scope and frozen sources

The repair removes ordinary allocations from the sharp variance-one decision
without changing the physical rate, sampler, RNG order, native coefficient
rounding, public API or error-publication contract.

- Frozen parent LW: `HybridPlasmaPIC-local-coulomb-rate-20261026`, on main
  `861e562d14eb2b80513ea40aa4ab3efb2fc21d94`.
- Source-identical dependency snapshot:
  `b58e3d3b9cce6a3d0a4e9675529ee6e6f5451c40`.
- Frozen local-only V2: `43a24e1f97874eaeab39487d27c9d5862c35557b`, in
  `HybridPlasmaPIC-local-coulomb-rate-v2-20261026`.
- Held V3: `HybridPlasmaPIC-local-coulomb-rate-v3-20261026`, based on accepted
  `85e75f2eb40d4ad055ee46c5c9db800d3b17ffa3`. Its shared sampler and finite-draw
  regression are byte-identical to that commit. All other production sources
  match V2. Root shared plans and trackers remain inherited from accepted main.

Evidence lives under
`HybridPlasmaPIC-audit-evidence/2026-10-26/local-coulomb-rate/`:
`optimization-v2/` retains the pre-edit proof, baseline, failures, source archives
and optimization-only patch; `composition-v3/` contains the composed campaign.
`composition-v3/assembly-v3.json` pins all 307 Julia source/test files.

## Pre-edit contract and verified cause

The exact pre-edit plan is retained as `optimization-v2/pre-edit-plan-v2.md`.
Its recovery from the recorded checklist edits matches the original pre-edit
SHA; `pre-edit-plan-retention.json` records that check and method.

The six failing Float32, four-species, N=4096 profiles cover dimensions 1–3 and
both tuple/vector collections. On both runtimes, their first measured calls
allocate 33,448 / 34,080 / 34,712 bytes by dimension, followed by zero-byte calls.
The original gate remains all 108 profiles with three zero-byte samples each;
no profile exclusion or additional warmup is permitted.

Source inspection and allocation stacks identify BigInt work in the local
original-factor/mathematical-pi comparator. The shared sampler invokes that
comparator inside the Float32 roundoff guard even when a binary64 interval can
certify the side of one. Copied-RNG profiling also observes separate integer
buffer growth; its totals are not the original `@allocated` measurements.
`baseline-trace-{112,1132}.toml`, logs and raw stacks retain the distinction.
The baseline records 18 trace-equivalence assertions per runtime.

The repair extracts the existing scaled channel, count/reuse and weight/speed
operations into common helpers. Existing wrappers preserve operation order and
native narrowing points. Only the Float32 local near-one decision reevaluates
all original factors in binary64. Unresolved cases retain exact rational factors
and certified mathematical-pi comparison. No cache, fixture branch, MPFR setting
change or shared sampler edit is introduced.

The old nine-argument private mass-update test call is migrated to the canonical
seven-argument finite-increment owner. Its exact independent `gp-g` oracle and
tolerances are retained, with finite and exact-representability assertions added.

## Certificate derived before implementation

The original physical variance is

```text
dt * number_scale * logLambda * (qa*qb)^2 * (ma+mb)^2 * wmax * count
----------------------------------------------------------------------------
8*pi * epsilon^2 * ma^2 * mb^2 * prod(dx) * wbar * reuse * triplet * speed^3
```

Here `triplet` is two for an odd self triplet and one otherwise; self count is
N-1, cross count/reuse are the original integers, and speed is the stored
`max(hypot(g), u_floor)`. Zero-weight particles do not enter the bins.

Float32 factors convert exactly to Float64, including subnormals. `frexp` carries
exponents separately; bounded mantissa products and quotients remain normal.
Integer conversions and Float64(pi) each contribute one rounding factor. The
mass expression `(1+min/max)^2/min^2` includes four factors for the squared
ratio-plus-one correction. The Float32 mass ratio remains normal in Float64.
The conservative n=32+D operation count covers the channel, count/reuse and
weight/speed calculation; removing Float32 narrowings reduces that count.
The fast certificate accepts only Float32 weight/speed and Float32/64 mean
arguments, whose conversions are exact. Other private types retain rational
fallback.

With u=2^-53, exact positive variance v and binary64 reevaluation vh,
`(1-u)^n <= vh/v <= (1-u)^(-n)`. Independent Python Fraction inequalities in
`certificate-inequalities-v2.json` prove both endpoint deviations are strictly
less than g=(64+4D)*eps(Float64), for D=1–3. Thus vh>1+g certifies v>1 and vh<1-g
certifies v<1. Subtraction near one is exact by Sterbenz. A final `ldexp`
overflow/underflow has an unambiguous side because its mantissa/exponent remain
bounded and accurate. Inside or on the guard, original rational factors and
Machin pi bounds resolve the decision. The source-level derivation assumes the
existing hardware contract: round-to-nearest and retained subnormals. It is not
a formal compiler or hardware proof.

Independent tests reject widening an already-rounded coefficient, native-pi
substitution, rounded integer counts, removal of rational fallback, changed
sampling/casts and relaxation of the 108-profile allocation gate.

## Implementation and verification tracker

- [x] Inspect the canonical local owner, shared sampler boundary, original guard
  derivation and raw allocation failures; freeze source-identical baseline.
- [x] Derive and check the binary64 certificate before editing.
- [x] Extract canonical scaled arithmetic and add the certified local decision.
- [x] Migrate the stale private test call while preserving its independent oracle.
- [x] Preserve V2 and assemble accepted finite-draw repair 85e75 into separate V3.
- [x] Revalidate physical boundaries, range, normalization, fallback and settings.
- [x] Verify complete particle-state/RNG histories and unchanged ordinary bytes.
- [x] Review concrete compiler entries, Aqua, full raw compiler trees and format.
- [x] Complete final neighboring tests and executable examples on both runtimes.
- [x] Audit the final source and evidence; freeze patch/docs in the final handoff.
- [x] Repeat verify-only cost preflight with the explicitly authorized reference.
- [ ] Run the unchanged fixed-cost campaign in a coordinated exclusive window.
- [ ] Parent integration, individual main push and exact-commit GitHub CI.

Results per runtime unless stated otherwise:

| Gate | Checked result | Evidence in `composition-v3/` |
| --- | --- | --- |
| Boundary/certificate/resources | 4,472 assertions; 108 profiles x 3 samples = zero bytes | `focused-*` |
| Original physical-factor matrix | 9,003 assertions, 3,600 rows; all original 23 finite-Gaussian failures now pass with accepted 85e75 | `factors-*` |
| Independent wide-factor matrix | 4,323 assertions over 1,440 valid-grid rows | `wide-*` |
| Local and accepted finite-draw suites | 1,860 local + 1,338 finite-draw assertions | `integration-*` |
| Complete local state/RNG differential | 2,916 assertions over 108 ordinary + 72 exceptional profiles; all calls succeed | `differential-*` |
| Prescribed-rate state/scratch/RNG | 16,704 assertions; 576 profiles x 4 states, equal outputs and zero bytes | `prescribed-*` |
| Prior range and normalization | 139,319 main + 7,690 isolated checks = 147,009 | `range-*` |
| Concrete compiler/Aqua | 127 concrete JET + 2 return-type + 9 Aqua checks | `quality-*` |
| Complete compiler trees | 29 minimum / 27 latest reports identical to reviewed reference after checked unchanged-line mapping | `complete-compiler-review-v3.json` |
| Neighboring regressions | 49,726 assertions | `neighbors-*` |
| Executable example | Completes successfully on both runtimes | `example-*` |
| Final evidence audit | 9,017 source/evidence invariants | `final-evidence-audit-v4.json` |
| Read-only formatter | 81 assertions over 20 owned Julia files, including AST/idempotence/source hashes | `format-readonly-v3.*` |
| New reference cost preflight | 104 mutations + 56 ordinary and 30 mutations + 8 exceptional profiles | `*-preflight-v4-*` |

The local differential reference is immutable unoptimized LW with the same
accepted 85e75 sampler. It has not been replaced by the optimized implementation.
Full compiler comparison means no new reports, not a zero-report package.

The independent exact Fraction replay uses a different identity,
`pi = 4*(atan(1/2)+atan(1/3))`, with 160 alternating terms and exact remainder
brackets. It certifies all 1,440 wide and 3,600 original factor rows per runtime,
checks original input/approximate-variance dictionaries against retained V1,
and rejects missing-row, one-ulp input and perturbed-output controls. The full
compiler auditor independently rejects changed-message, changed-frame and
missing-report controls per runtime. Both reviewer scripts are SHA-pinned.

## Retained failures and changes to verification

- V1's six allocation failures per runtime and stale private-call error remain
  retained. V2/V3 pass their unchanged allocation gate.
- V1's 23 finite-Gaussian failures per runtime were repaired by the separately
  accepted 85e75 sampler. The unchanged 3,600-row matrix passes after composition.
- The V2 extreme-factor generator stopped with 27 passes / 1 fixture error per
  runtime because length 2.2336249e-38 violates the canonical FourierGrid domain.
  V3 preserves proposed lengths and filters only the actual four-point-grid
  requirements: finite positive dx, scale and maximum wavenumber 2*scale.
- Initial exact replay V3 decoded Float32 decimal strings as Float64. Its files
  remain retained; replay V4 correctly decodes the declared binary32 type and
  is the accepted evidence checker.
- The original historical-reference preflight remains failed as described below.

## Fixed-cost contract and authorized reference decision

No acceptance timing has run. The frozen contract remains 304 timing rows per
runtime: 72 legacy + 216 local + 8 ordinary controls + 8 exceptional/control rows;
8 warmups, 20 ms parent/control calibration, 21 alternating samples, complete
state/scratch/RNG replay, zero measured compilation and strict byte gates.
Identical controls must satisfy both metrics in [0.8,1.25]; doubled controls must
exceed 1.5; legacy ratios must be <=1.5. Local workspace payload remains
`3*sizeof(T)*sum(N) + (2*sum(N)+S*(2*cells+1))*sizeof(Int)`, with a
1024*(S+1)-byte header allowance. No released equivalent physical-rate speed
baseline exists for the new local operator.

The unchanged historical 1097076 reference failed verify-only preflight on both
runtimes at row 9: D1 Float32 n=16, uniform weights, gcoeff=0.5, three steps.
Complete output/scratch/RNG differs despite zero bytes/compilation. An independent
immutable-85e75 replay passes 13 checks per runtime and equals V3 after every
step, while reproducing old-reference differences of 26/35/38 velocity
components on minimum and 26/36/37 on latest. No V3 production edit followed.

Root explicitly authorized a new versioned campaign using accepted
85e75f2 as the prescribed/legacy numerical, bitwise, allocation and latency
parent. Only the legacy harness include changes; all workload inputs,
mutant/control cases, row counts, byte/compilation gates and 1.5 limits remain.
Local and other component references remain unchanged.
`cost-v4/BASELINE_DECISION.md` and `baseline-decision-v4.json` record the
premeasurement decision and exact hashes. The new verify-only preflight passes
56 ordinary and 8 exceptional profiles on both runtimes, including unchanged
104 + 30 mutation checks. The original failed receipts stay failed.

Root must coordinate timing and continuously monitor foreign heavy workload
arrivals using process identities/start times. An arrival invalidates the
campaign; preserve it without signaling a foreign process or silently retrying.
Promotion and exact-commit CI follow acceptance; this handoff makes neither claim.
