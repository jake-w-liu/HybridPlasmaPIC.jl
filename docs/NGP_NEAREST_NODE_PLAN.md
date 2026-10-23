# NGP nearest-node boundary repair

Status: 2026-10-23; final local gates pass, pushed as `df75656`; CI remains.

## Pre-edit brief

- Mode: repair the shared particle-to-grid shape policy.
- Requirement: NGP selects the nearest represented node, with exact half-cell
  ties toward the larger index. CPU and Metal must use the same policy.
- Baseline: against exact `2993f5a`, Julia 1.12.7 and 1.13.2 each retain 118
  passes and 14 failures, including public gather/deposit in dimensions 1–3.
  The exact Metal helper AST retains 139 passes and five disagreements per
  runtime. Actual Metal 1.11.1 execution retains 77 passes and 13 failures.
- **VERIFIED cause:** `src/coupling/shapes.jl:40` forms `s + T(0.5)` before
  flooring; `prevfloat(T(0.5))` consequently selects node one instead of zero.
  `ext/HybridPlasmaPICMetalExt.jl:63` uses nearest-even for NGP, whereas
  `test/test_deposit.jl:99-110` requires one-cell translation of half-cell ties.
- Contract: retain unit shape weight, upward ties, periodic indexing, public
  APIs, validation and input/output ownership. CIC and TSC retain their policies.
- Canonical owner: a private NGP rounding helper in `src/coupling/shapes.jl`,
  used by CPU `_stencil1d` and Metal `_stencil_base`. CPU scalar/vector gathering,
  deposition, density, current, momentum and pressure moments already share
  `_stencil1d`; Metal scalar/vector kernels share `_stencil_base`.
- Compatibility: change only incorrectly assigned NGP boundary points. Keep
  index conversion and validation at the existing CPU/device boundaries.
  CUDA currently provides backend storage/telemetry, not a separate stencil.
- Strategy: delegate the rounding policy to Julia's `RoundNearestTiesUp`, whose
  controls pass the independent exact-rational baseline oracle on both runtimes.
  Share its floating result before the existing backend-specific integer conversion.
- Verification: exact rational node oracle across ties, adjacent floats, signs,
  exponent ranges and independent random samples; public scalar/vector gather,
  deposition and moments; translation, partition and adjoint identities;
  unchanged failure paths, complete coupling suites, concrete/package compiler
  checks, formatter and real Metal kernel execution where supported.

## Acceptance gates

- [x] Retain failing CPU and exact-helper baseline commands on both runtimes.
- [x] Observe the actual Metal baseline and classify any unavailable check.
- [x] Implement the common policy with permanent CPU and optional-device tests.
- [x] Pass numerical, integration, failure, compiler and formatting gates.
- [x] Measure complete public NGP calls against the parent: ordinary outputs and
  allocations must match; all fixed profile median latency ratios must be <=1.5.
  Boundary output differences must match the independent oracle.
- [x] Review exact source, preserve held changes, commit and push separately.
- [ ] Require all six exact-SHA GitHub jobs green.

Three tempting wrong patches are rejected explicitly: nearest-even breaks
half-cell translation; fixing only CPU leaves the Metal tie mismatch; changing
all shapes to upward ties changes TSC's established centered stencil. Boundary,
device and CIC/TSC controls distinguish them. Reintroducing `floor(s + 0.5)`
is rejected by the immediately lower half-cell neighbor.

Device coordinate wrapping and nonfinite validation remain separate audit
questions. This repair does not establish those contracts or GPU performance.

Persistent baseline evidence:
`../HybridPlasmaPIC-audit-evidence/2026-10-23/metal-shape-boundaries` and
`../HybridPlasmaPIC-audit-evidence/2026-10-23/new-baseline-witnesses.json`.

## Observed first candidate

The actual Metal baseline on Julia 1.13.2 / Metal 1.11.1 retains 77 passes and
13 failures, with no fixture errors. The common-policy candidate passes 7,472
focused CPU assertions on each Julia 1.12.7 and 1.13.2 and 516 actual latest-Julia
Metal-device assertions, including exact NGP outputs and tie translations.
Formatting passes ten syntax/idempotence assertions. An earlier formatter
harness incorrectly passed a string style keyword; its failed log is retained,
and the corrected driver delegates configuration discovery to JuliaFormatter.

The first combined candidate was composed onto accepted `7b50b51`. The four
production/test blobs were copied without modification, and both test includes
were added to the new parent's runner. Its combined gates and final composition are
recorded below. Original evidence is retained in
`../HybridPlasmaPIC-audit-evidence/2026-10-23/ngp-nearest-node`; combined evidence
is in the sibling `ngp-nearest-node-combined` directory.

## Final composition and untimed gates

The candidate is composed onto `bf2a8bc`, preserving the separately accepted
Metal uniform-push kernel. Source comparison against the prior combined NGP
candidate finds one additional CPU source difference: removal of the unused
recursive private Metal helper. That exact removal is already accepted in the
parent; the canonical CPU push body is unchanged. All other CPU source bytes
match the prior 49,590-assertion coupling/solver integration snapshot on each
runtime. Its full compiler comparison also passes.

Fresh final-composition runs on Julia 1.12.7 / 1.13.2 each pass 7,472 numerical
and public-invariant assertions, 46 concrete compiler checks and 1,410 actual
Metal 1.11.1 backend integration assertions. Both full package compiler trees
match `bf2a8bc` after checked source-line maps (30 / 27 existing reports).
Formatting passes ten parsed-syntax/idempotence assertions; 276 Julia files are
frozen. All 126 ordinary owner profiles per runtime preserve complete outputs,
scratch and inputs and allocate the same 0 bytes as the accepted parent in
three warmed samples (630 assertions per runtime). The separate fixed timing gate now passes on both runtimes.

The previous combined GPU command retains its parent uniform-push stack
overflow. Its post-overflow test results are not acceptance evidence. The fresh
final GPU commands run in separate processes and pass the complete backend,
uniform-push and NGP suites together.

Final evidence directory:
`../HybridPlasmaPIC-audit-evidence/2026-10-23/ngp-final`; see `composition.json`,
`native-*.json`, `preflight-summary.json` and `compiler-comparison.json`.

## Measured acceptance

Each version passes all eight identical/doubled-work timing controls and all
882 resource/state assertions across the fixed 126-profile matrix. Every
profile retains three 0-byte allocation samples and complete ordinary outputs.
Twenty-one alternating timing pairs follow a parent calibration of at least
20 ms. Maximum ratio-of-medians is 1.210119 on Julia 1.12.7 and 1.195323 on
1.13.2; maximum median paired ratio is 1.187697 / 1.199248. All are within the
predeclared 1.5 gate; this is a correctness repair with measured overhead, not
a speedup claim. Both identical-owner controls pass and doubled-work controls
exceed 1.5 on both statistics. Raw samples remain in `costs*.toml`.

An additional private-scalar compatibility probe passes 204 Float16/BigFloat
half-cell/neighbor assertions on each runtime. It does not claim complete
solver support for those scalar types. `final-verification.json` checks all
15 final commands, 276 frozen Julia files, 356 unchanged committed Julia files
and held hashes. Promotion verifies all 361 committed Julia files against the
validated candidate and unchanged held hashes; remote main matches
`df75656390ed416f7b0a37cb1173058cd9465e71`. All six jobs in the
[exact-SHA GitHub run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39912564815)
remain required.
