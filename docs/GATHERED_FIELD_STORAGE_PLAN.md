# Gathered field inputs and particle outputs

## Repair brief

The public `push_gathered!` treats gathered E/B as fixed inputs while writing
particle positions and velocities. Detectable overlap with those output arrays
must reject before mutation on CPU and Metal. E/B may share each other or
unchanged weight/metadata storage. Midpoint output admission retains its
separate [contract](GATHERED_MIDPOINT_STORAGE_PLAN.md).

**Mode:** repair. On the clean parent, v₁=[1,2,3], E₁=reverse(v₁), B=0 and
dt=0.5 gives v₁=[2.5,3,4.25] instead of the fixed-input [2.5,3,3.5]. Both
Julia 1.12.7 and 1.13.3 reproduce it for Float32/64 and D=1–3. The CPU loop
reads E₁[3] after writing v₁[1]. Evidence is retained at
`takeover/gathered-input-alias/before-v1-*` in the dated external audit archive.

**Canonical owner:** `_check_gathered_field_storage` in `src/particles/boris.jl`
already precedes publication in both the public CPU method and Metal extension.
Extend it with `_particle_arrays_mightalias` over x/v only. No kick, drift,
device broadcast, buffer, public signature or archive changes. Rejection is the
compatibility decision for fields sharing writable particle storage. CPU views
that the predicate proves disjoint and read-only sharing remain supported.
Device views retain the canonical predicate's
conservative behavior.

The check is bounded by dimension and field component count, independent of
particle count. That establishes its loop structure, not its measured latency.
No allocation or speed improvement is claimed. Detection follows the existing
predicate and does not promise detection of hidden unsafe pointer aliases.

## Acceptance tracker

- [x] Preserve actual clean-parent failures on both supported runtimes.
- [x] CPU E/B components against every x/v family, identity and distinct views,
  empty/nonempty populations, zero/nonzero dt, with/without midpoint output.
- [x] Explicit offending-argument error and unchanged arrays after rejection.
- [x] Shared read-only fields/weights and disjoint views remain valid; compare
  all public state against a frozen-field reference.
- [x] Actual Metal rejection/valid sharing on both supported runtimes.
- [x] Neighboring particle/range tests, matched allocation/state profiles,
  concrete analysis, Aqua and complete compiler-tree review.
- [x] CI-pinned formatter AST/idempotence and exact final diff.
- [x] Root integration onto current main, including both supported runtimes.
- [ ] All six exact-SHA GitHub CI jobs; current status belongs to the [audit tracker](DEEP_DEBUG_TRACKER.md).

Wrong repairs include identity-only checks, a velocity-only check, checking
only E, rejecting shared read-only fields, and checking after the kick. The
all-component/view matrix, state snapshots and safe-sharing controls distinguish
them. The original reversed-view numerical oracle remains an independent
reproduction; it is not replaced by a hard-coded production result.

## V2 allocation repair, fixed before source edits

V1 is preserved with its full candidate tree and failed gate. On Julia 1.12.7,
36 of 144 matched admission profiles regress: mixed E-component storage adds
496, 544 or 592 bytes for D=1, 2 or 3. The other 108 match the clean parent.
The baseline itself allocates on mixed tuples, so the gate is equality with the
matched parent, not a zero-allocation claim.

A controlled five-variant experiment retains the exact dimensions and midpoint
checks. Splitting the outer E/B loop still regresses; hoisting its indexed
components leaves 240 extra bytes. Passing each of the six literal tuple
components to a motion-component helper, scanning homogeneous x and v families
separately, matches all 144 parent profiles. Explicit inlining gives the same
result and is unnecessary. Allocation traces identify tuples and SubArray boxes
in the original mixed-component loop. The exact probe sources, complete raw
profiles and optimized typed IR are in the external audit archive under
`takeover/gathered-input-alias/campaign-v2/allocation-*` and `typed-*`.

Implement that component helper in the canonical admission path; do not change
particle alias predicates or mover arithmetic. Add distinct empty CPU views and
native nonempty contiguous Metal views to the rejection matrix, including zero
dt. Metal tests must execute the specialized native-array method, rather than a
CPU fallback. Distinct empty device wrappers require a separate measured
canonical-predicate review; no new hidden-alias guarantee is made here.

Frozen gates: all original 144 admission profiles must match the parent on both
runtimes; original and supplementary whole-public allocation/state profiles must
match their exact parent; focused CPU and actual Metal tests must pass; all
concrete compiler reports must be empty, Aqua must pass, and complete package
compiler trees/messages/frames may differ only by exact checkout-root replacement
and explicitly checked source-line movement. Compare the full executable mover
AST outside the two admission helpers, including mutation controls. Apply the
CI-pinned formatter and verify two idempotent passes. Record all raw results,
including failures, before assessing acceptance. Resource claims use allocation
bytes, not wall-clock timings under concurrent host load.

Supplementary gates fixed before their measurements: compare exact parent
admission AST with the canonical candidate on actual native Metal storage in
36 profiles per runtime (D=1–3, N=0/16, midpoint absent/present, independent
fields/read-only weights/shared fields). Require identical three-sample host
allocation bytes; this is not a device-memory or whole-device-push profile.
Independent admission mutants must distinguish identity-only, velocity-only and
E-only checks while retaining read-only sharing. These light follow-on checks
run serially after the full compiler campaign.

## Local validation completed

The isolated V2 candidate passes on Julia 1.12.7 and 1.13.3. Each runtime has
35,971 focused CPU assertions, 3,135 actual-Metal assertions, 144 exact matched
CPU admission profiles, 198 exact whole-push allocation/state profiles, and
36 exact native-Metal admission host-allocation profiles. All three allocation
samples match the parent in every profile. The supplementary whole-push matrix
includes mixed B components, tag-backed fields and different x/v view types.
Independent identity-only, velocity-only and E-only admission mutants are
separated by 18 witnesses; six read-only controls also pass per runtime.

All 90 concrete compiler reports per runtime are empty and Aqua passes. Complete
package reports contain the same 29/27 trees as the matched parent, including
all messages, unlimited rendered output and every recorded virtual-frame field.
Only exact checkout roots were replaced for comparison; no source-line,
message, frame or tree filtering was needed. Full executable mover AST equality
outside the two admission helpers passes, including negative mutation controls.
JuliaFormatter 1.0.62 passes two idempotent whole-tree formatting checks.

CPU read-only metadata tests use the generic vector API. Native Metal field
storage remains Float32; the actual device suite checks weight sharing, shared
fields, identity aliases and nonempty contiguous native-array view aliases,
with scalar indexing disabled. On both runtimes a distinct empty contiguous
Metal wrapper shares device data identifiers but the canonical alias predicate
returns false; its nonempty counterpart is detected. This is a measured
predicate limitation for zero-length storage, where no particle elements are
written. The repair does not broaden that predicate or claim detection of all
independently wrapped device pointers. Predicate-proven disjoint CPU views
remain supported; conservatively overlapping device identifiers need not be
admitted as disjoint.

The external receipt is
`takeover/gathered-input-alias/campaign-v2/final-v3/comparison-v1.json`, with all
38 checks true. The main campaign has 27 successful commands and the declared
light follow-on campaign has four. Complete failed V1 sources and evidence,
allocation-isolation experiments and two preliminary wrapper failures remain
preserved; wrapper corrections did not change the frozen production candidate.
At this isolated handoff, root integration, the separate main push and
exact-SHA GitHub CI remained pending. No full `test/runtests.jl` or latency
benchmark is claimed for this local slice.

## Current-main integration

Root integrated the exact sealed admission helper and both new permanent tests
onto `fd0a168`, retaining the landed near-light, midpoint, BigFloat, metadata
and RK4 owners. All 27 commands pass on Julia 1.12.7 and 1.13.3 with
registry-only SpectralOperators 0.1.6 and unchanged transitive dependency graphs.
Each runtime passes 57,632 focused and 3,135 actual-Metal assertions, 144 exact
admission profiles, 198 exact whole-push state/allocation profiles and 36 zero-byte
native-Metal admission host profiles. All original allocation samples match.

Aqua and all 90 concrete analyses pass. Complete 29/27-report package output,
unlimited rendering and 252/203 actual virtual frames match the current parent
with only exact checkout-root replacement. No line, label, path or diagnostic
filter is needed. Eight changed-message/render/report/frame/line/signature/path
reviewer controls reject per runtime. All 753 frozen source/harness pins remain
unchanged; executable arithmetic outside admission is AST-equivalent. The pinned
formatter preserves all Julia source bytes and is idempotent.

Receipts are retained in `takeover/gathered-field-root-integration/`, including
`integration-review-v1.json` and `complete-compiler-root-only-proof-v1.json`.
The actual device empty-wrapper limitation above remains unchanged. No whole
`test/runtests.jl`, latency improvement or complete resident-solver claim follows
from this focused integration. Exact-SHA GitHub CI remains an acceptance gate.
