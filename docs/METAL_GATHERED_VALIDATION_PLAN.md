# Metal gathered-field preflight repair

Status: baseline confirmed on Julia 1.12.7 and 1.13.2 with actual Metal 1.11.1;
V1 rejected by the fixed allocation gate; V2 rejected by latency; V3 is untested.
Original parent: `ff4a973` (Julia sources match `df75656`); V3 starts from
`61356f1`, with identical committed Julia sources. Persistent evidence: `../HybridPlasmaPIC-audit-evidence/2026-10-23/metal-gathered-validation`.

## Pre-edit brief

- Mode: repair backend validation at the gathered-push boundary.
- Requirement: the public Metal overload rejects invalid charge/mass/timestep
  and nonfinite gathered fields before writing any particle or midpoint output,
  matching the existing CPU contract. Empty populations still validate host
  parameters and array shapes/storage.
- Baseline: each runtime retains 114 passes and 51 failed assertions, including
  45 passing finite CPU/device controls; no fixture errors. The existing CPU
  rejects all invalid cases and preserves complete state.
- **VERIFIED cause:** `ext/HybridPlasmaPICMetalExt.jl:643-675` uses raw `q/m`,
  converts `dt` without validation and proceeds to mutating broadcasts without
  a finite-field scan. CPU `src/particles/boris.jl:165-177` preflights these
  values; `test/test_particles.jl:246-290` explicitly exercises rejection and
  late-invalid-field rollback. Actual public baseline logs prove the mismatch.
- Canonical owners: existing `_validated_qm`, `_require_finite_real`, gathered
  length validation and scalar finiteness policy. Implement the device traversal
  at the Metal overload; reuse host policy rather than re-derive it.
- Compatibility: preserve signatures, ordinary broadcast arithmetic, storage,
  midpoint behavior and shape/backend error ordering. Invalid inputs now throw
  `ArgumentError` before mutation. The field error names the combined fields;
  exact message text is not a published contract. No new alias restriction or
  update-kernel rewrite is bundled into this repair.
- V2 strategy: after shape/storage checks, validate host values; empty valid sets
  return without device work. Scan all six fields in one private device kernel,
  atomically flag invalid samples, and transfer only the four-byte status before
  existing writes. UInt64 striding retains all host-representable indices while
  limiting launch coordinates to Metal's UInt32 intrinsic range.
- Verification: real-device baseline and expanded late/first/block-edge sample
  failures, complete input/output snapshots and recovery, CPU controls,
  parent differential finite outputs, optional midpoint, empty populations,
  D1–3, concrete compiler and full backend integration on both runtimes.

## Fixed acceptance gates

- [x] Preserve both actual-device baseline failures and CPU controls.
- [ ] Implement the preflight and permanent regression tests.
- [ ] Pass original 165 assertions and expanded rejection/rollback controls;
  all six fields, NaN and both infinities, invalid host values including conversion
  overflow, zero/negative timesteps as valid controls, length/storage failures,
  empty sets, repeated valid steps and successful reuse after rejection.
- [ ] Retain ordinary parent outputs and field inputs, including supported
  identical-array aliases and midpoint behavior; no user-visible reorder.
- [ ] Pass full actual Metal backend/NGP/uniform/gathered suites, CPU particle
  neighbors, formatter and applicable concrete/package compiler comparison.
- [ ] Fixed cost matrix: D1–3 × N=0,1,257,4096,262144 × midpoint absent/present
  (30 profiles per runtime), fixed nonzero magnetic rotation and zero E.
  Compare complete synchronized calls and complete outputs with the parent.
  Warm/compile first, parent-calibrate to at least 20 ms, then 21 alternating
  pairs in an exclusive window. Both ratio-of-medians and median paired ratio
  must be <=1.5; identical/doubled controls validate the protocol before acceptance.
- [ ] The newly required validation may allocate metadata: additional warmed
  host bytes are limited to 8 KiB per complete call, independent of particle
  count, compared with the parent. This is an explicit cost allowance for the
  validation, not a zero-allocation claim. Record all profiles and
  device reduction storage where the actual backend API exposes it; unavailable
  pool counters cannot establish a measured device-memory bound.
- [ ] Inspect exact final diff, unchanged CPU sources, all committed candidate
  blobs and held hashes; commit/push separately and require all six exact-SHA
  GitHub jobs green. Hardware checks remain separate from CPU-only CI.

Wrong patches rejected by the gates: host-only validation misses device NaNs;
checking fields after the push violates complete-state rollback; only checking
first samples misses late/block-edge invalid values; returning early before host
validation accepts invalid empty inputs; copying complete fields to the host
violates the traversal design and measured scaling requirement.

The existing broadcast push still duplicates scalar Boris algebra and allocates
intermediate device arrays. That is a separate optimization candidate; no
performance claim is inferred before profiling and compatibility review.

## Retained revisions

V1 passes the original 165 and complete 6,394 integration assertions per runtime.
Both fixed untimed cost matrices retain 120 alias passes and 234 passes/six
allocation failures: largest profiles add at least 9,344 host bytes (limit 8,192).
Full allocation profiling attributes all additional sampled bytes to the
reduction; unrelated allocation totals match the parent exactly.

V2 replaces that reduction with one per-call atomic status flag. Its independent
prototype passes 344 assertions per runtime, including forced repeated strides
and all-invalid fields; complete public integration, compiler and fixed costs
remain required. See retained `pre-edit-v2.md`, `allocation-attribution.json`,
`cost-preflight*.log`, and `atomic-probe*.log` in the evidence directory.

## V2 retained acceptance results and rejection

Both Julia 1.12.7 and 1.13.2 pass the original 165 assertions and 6,426 full
integration assertions. Each passes 120 parent-alias checks and 240 untimed
cost/state checks across 30 profiles. Maximum additional host bytes are 4,048
and 4,016, below the fixed 8,192 limit. Full package compiler trees match the
parent after source-root replacement (30/27 reports). Each complete device
call retains 24 matching parent diagnostic leaves plus three status-transfer/
kernel-launch reports in Metal/GPUCompiler. Raw signatures, complete caller
paths, dependency source and actual-device executions are retained; these
reports were not suppressed or asserted to be absent.

The minimum-version timing matrix passes all identical/doubled controls but
retains 254 passes and 46 latency failures across the fixed 30 profiles. V2 is
rejected; the latest timing lane is unrun. All raw samples are preserved in
`timing-candidate3-112.toml`, with the frozen source and command manifest.
Task agents held native work; unrelated user-owned workloads remained active
and were recorded without interruption. This is not an idle-machine claim.

A separate six-profile intervention adds a per-call fence to the parent and
preserves all 54 complete-state checks. Its fenced/queued median ratios range
from 1.26 to 2.49, while the candidate/parent ratios range from 3.68 to 4.05. This
is diagnostic evidence, not isolated acceptance or a complete causal
attribution. No threshold changed and no same-source retry was accepted.

## V3 proposed guarded device commit

V3 retains the existing Boris broadcast expressions and their evaluation order.
It initializes a per-call status on the device, scans all six fields, then queues
conditional output copies. Each copy kernel returns without a store when status
is invalid. A single final four-byte host transfer reports rejection. Valid
calls complete their device work before returning. Base.Broadcast preparation
retains scalar expansion and alias-copy semantics; no alternate Boris formula
or new public alias restriction is introduced.

The installed Metal source routes device fill, scan and copy launches through
the same default queue; the host read flushes and waits for that queue. Runtime
ordering, invalid-input rollback and compatibility still require verification.
This is an untested implementation in a separate worktree, not an accepted fix.

Additional gates: instrument actual stores so a rejection cannot pass merely
by writing old values; verify that the instrument detects those writes; expand
cross-component and shared-storage alias comparisons. Then rerun the original
behavior, full hardware/compiler and unchanged allocation/latency matrices on
both runtimes. Every checkbox in the final acceptance list remains open until
that final source passes. `pre-edit-v3.md` retains the full pre-edit brief and
wrong-implementation challenges in the evidence directory.
