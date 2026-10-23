# Metal gathered-field preflight repair

Status: baseline confirmed on Julia 1.12.7 and 1.13.2 with actual Metal 1.11.1;
V1 rejected by the fixed allocation gate; V2, V3 and V4 rejected by latency;
V5 fails its fixed latency gate. V6 passes native gates and minimum-Julia costs,
but latest-Julia timing calibration fails; no gathered-push candidate is accepted.
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
The original V3 proposal is retained here; its completed results follow below.

Additional gates: instrument actual stores so a rejection cannot pass merely
by writing old values; verify that the instrument detects those writes; expand
cross-component and shared-storage alias comparisons. Then rerun the original
behavior, full hardware/compiler and unchanged allocation/latency matrices on
both runtimes. Every checkbox in the final acceptance list remains open until
that final source passes. `pre-edit-v3.md` retains the full pre-edit brief and
wrong-implementation challenges in the evidence directory.


## V3 completed checks and rejection

Both runtimes pass the original 165 assertions, 6,459 full integration checks,
120 original alias checks and 240 resource/state checks. Additional tests pass
1,944 cross-component/overlapping-storage comparisons and 90 full-buffer
rejection checks per runtime. Actual store counters distinguish a rejected
no-store operation from a write of the previous value. Maximum additional host
bytes are 1,088 / 1,040 within the unchanged 8,192 limit.

Full package compiler trees match (30/27 reports). Every concrete device call
retains 21 matching typed parent leaves, three changed fill-kernel
specializations at the same upstream sites/signatures/messages, and four added
upstream kernel/status-transfer reports. Full captures and their explicit
classification are retained in `compiler-candidate4-structural-comparison.json`
and `compiler-candidate4-review.md`; no diagnostics are suppressed.

The minimum-version fixed 30-profile timing matrix passes its identical/doubled
controls but retains 256 passing assertions and 44 latency failures across 23
profiles. V3 is rejected; latest-version timing is unrun. Original samples,
source hashes and rejected revisions remain preserved.

## V4 disjoint-storage kernel

V4 starts from the accepted CN commit `ce52836`. For disjoint storage, a single
guarded per-particle kernel calls the existing canonical `boris_kick` and writes
velocity, midpoint and position. For overlapping storage, the existing ordered
guarded-broadcast body is retained, preserving cross-component/view semantics.
`Base.mightalias` selects the storage strategy before device updates. Both paths
share the same host validation, complete finite-field scan and final status
read. No new public alias restriction, exception fallback or copied scalar
Boris formula is introduced. The pre-existing broadcast algebra remains a
separate consolidation concern.

This revision remains unaccepted. Its new gates include exact parent-bit
comparisons with heterogeneous nonzero fields and actual store instrumentation
for the fused kernel. All existing hardware, compiler, allocation and unchanged
latency gates remain required. `pre-edit-v4.md` records the compatibility decision
and tests before editing; no speedup is claimed before measurement.


## V4 completed checks and rejection

Both runtimes pass 6,783 integration checks, 864 heterogeneous disjoint
comparisons, 865 range checks and 2,034 expanded alias/rejection checks, along
with the original 165 witnesses and 360 original alias/resource assertions.
Full package compiler trees match the CN parent (29/27 reports). All 12 complete
device captures are retained and classified, including the added fused-kernel
report. Ordinary host allocation decreases in every sampled profile.

The first minimum-runtime timing attempt stops in a parent-only doubled-work
control: one 0.830592-second calibration observation selects one repetition,
but subsequent parent batches are near 0.0013 seconds. No candidate profile is
measured. A separately registered five-observation median calibration passes
93 deterministic checks and all four real controls on each runtime. It retains
every original profile, 21 alternating pairs, allocation and median-ratio gates.
The original failed measurement and all earlier rejected revisions remain intact.

With that verified calibration, the complete minimum-version matrix retains
298 passes and two latency failures. D3/N1/midpoint has median ratios 1.798107 /
1.765316 against limits of 1.5; all other 29 profiles pass. V4 is rejected, and
latest-version timing is unrun. No V4 candidate is pushed.

## V5 single-group composition

V5 starts from `c810eb0`, retaining the accepted memory-status constructor repair.
When the existing launch geometry has one group and outputs are disjoint, a
small wrapper calls the existing field scan, a device-memory threadgroup barrier,
then the existing disjoint update. Every thread reaches the barrier after its
scan returns, including inactive threads and early invalid-field returns. Other
launches retain the separate scan and guarded update; aliases retain the ordered
broadcast path. The public owner allocates one private status flag per call.

This uses the existing launch group's count, with no new particle threshold,
fixture-specific condition, copied Boris formula or exception fallback. The
pre-edit contract is retained in `pre-edit-v5.md`. New instrumented tests exercise
actual stores across group sizes, inactive threads, forced strides and all six
invalid field components; valid zero-step controls still require real stores.
The 29 native commands covering both runtimes all pass. Each runtime records
165 original witness assertions,
16,935 integration assertions, 1,313 constructor assertions, 864 heterogeneous
trajectory comparisons, 865 range checks, 2,034 alias/rejection checks and 360
original alias/resource checks. The new single-group store instrumentation adds
10,152 assertions per runtime. Full package trees match (29/27 reports); all
12 device captures retain complete classifications of 30 versus 24 reports.

The unchanged calibrated 30-profile matrix passes all four parent controls on
Julia 1.12.7, then records 296 passes and four latency assertion failures:

| Profile | Ratio of medians | Median paired ratio | Limit |
| --- | ---: | ---: | ---: |
| D3, one particle, midpoint | 1.521545 | 1.532365 | 1.5 |
| D3, 257 particles, no midpoint | 1.550458 | 1.563715 | 1.5 |

All outputs and allocation checks pass. The runner stops before latest-version
timing as declared. V5 is rejected; it is neither retried nor promoted. Complete
commands, samples, frozen source and failure classifications are retained in
`metal-gathered-validation/candidate6-timing-rejection.json` and
`timing-candidate6-*` under the persistent audit evidence directory. A further
implementation requires a new concrete hypothesis and pre-edit contract; the
original profiles and acceptance limits remain binding.


## V6 shared status and retained timing failure

V6 uses a fresh initialized `MtlArray{UInt32,1,SharedStorage}` through Metal's
public constructor. The existing validators, scan, guarded disjoint update,
one-group barrier and ordered alias path remain the owners. Public copy routines
synchronize transfer and protect array lifetime; no unsafe host pointer, private
Metal API or shared global status cache is introduced.

The original new lifetime fixture left particle positions uninitialized. Its
minimum-Julia run recorded 462 passes / six failures, including NaN positions;
both original lifetime runs were invalidated. The corrected fixture explicitly
initializes every position, retains every original assertion and adds finite
initial/final state and actual-motion checks: 636 passes per runtime. This was
a fixture repair, not a production relaxation.

The final source composed onto `f850b84` passes **17,571 integration assertions
on each runtime**. The previously recorded exact, alias, range, store-instrumentation
and allocation gates also pass. Full package trees match the parent (29 / 27
reports). Each of the 12 complete device captures retains 30 reports versus 24
in the parent: three existing Metal launch leaves have the new kernel signature;
six additions are classified at GPUCompiler method-definition access and Metal's
optional-device branches for shared status copies. Every full typed path remains
in evidence; none is suppressed.

The frozen 30-profile matrix passes all four controls on Julia 1.12.7; all
candidate outputs, allocation bounds and latency gates pass, with maximum of
both ratio estimators **1.382485**. Julia 1.13.2 records four controls and **zero
candidate profiles**: doubled non-midpoint work has ratio of medians
**1.4907853205918302**, failing its required **>1.5** gate, although the paired
median is 1.9925238486913364. The runner exits 1 as declared. Latest-runtime
performance is unverified; V6 remains held. The failed calibration is preserved
in `candidate7-timing-review.json`, full samples and command logs. No rerun,
threshold adjustment or hardware-cause conclusion is accepted from these data.
