# Distributed FFT in-place roundtrip diagnostic repair

## Pre-edit brief

- Mode: repair the existing allocating roundtrip diagnostic.
- VERIFIED baseline: `distributed-fft-contract/roundtrip-inplace-baseline-v1-*`
  records three failures on both Julia 1.12.7 and 1.13.2 at accepted 85e75.
  FFT!, RFFT! and R2R! all throw MethodError at parent(input). Source inspection
  confirms their AbstractManyPencilArray container has first/last views but no
  parent method. The existing helper also allocates an independent default
  output, whereas the dependency requires identical input/output objects for
  in-place plans (`PencilFFTs/src/operations.jl::_check_arrays`).
- Canonical owner: distributed_fft_roundtrip_error in the Pencil extension.
  Dependency owns transform application, container views, plan compatibility and
  inverse normalization. No new public helper or storage type is needed.
- Requirement: default and explicit compatible output calls return the actual
  global maximum absolute roundtrip error for supported in-place and out-of-place
  transforms. Preserve ordinary out-of-place numerical behavior and allocation.
  Caller provides valid compatible storage on every communicator rank; arbitrary
  rank-local storage errors and MPI transport failure are outside this utility's
  existing contract. Collective constructor repair is separate.
- Strategy: use dependency is_inplace(plan) to choose input itself as the default
  output when required, and compare the first container view before/after the
  roundtrip. Preserve the existing copy, subtraction, norm and Allreduce order.
  Keep explicit invalid separate in-place output rejection in the dependency;
  do not copy it to make an invalid call silently work.
- Verification: original baseline witness, default/explicit output for complex,
  real-to-complex and real-to-real in-place plans and corresponding ordinary
  transforms, uneven/asymmetric shapes, rank-dependent fields and global-error
  reduction, exact manual-path diagnostic parity, first-view restoration, invalid
  output no-mutation. Run 1/2/4/8 ranks on minimum/latest; resource differential
  for unaffected ordinary path, concrete/full compiler and formatter gates.
  Freeze source before every campaign; separate commit/push and exact-SHA CI.
- Wrong patches: merely choosing first(input) while retaining independent default
  output; treating every input as a container; returning a canned small error;
  accepting invalid output by copying; comparing the output view after inverse.
  Default/explicit cases, manual full-state path, nonuniform fields and invalid
  storage tests distinguish these changes.

## Progress

- [x] Reproduce and retain parent failures on both runtimes.
- [x] Inspect dependency container and alias contracts and complete pre-edit brief.
- [x] Implement canonical view/default-output repair.
- [x] Pass complete numerical, compiler and resource gates.
- [x] Update docs and tracker and integrate after constructor repair.
- [ ] Push main and verify all six exact-SHA jobs and complete compiler trees.

## Retained compiler finding and integration

V1 passes all 1/2/4/8-rank numerical matrices on both runtimes. Its initial
concrete-analysis harness incorrectly assumed an empty diagnostic report. All
12 entries report one possible invalid index from PencilArrays transpose_recv!
at line 520, where MPI.Waitany's possible `nothing` return is used as an index.
The failure is retained: 9 Aqua assertions pass and 12 zero-report assertions
fail on each runtime.

Six corresponding unchanged-parent forward-transform probes independently
produce that same report on both runtimes. Full report tails, signatures and
messages match all 12 candidate entries exactly; the new caller frames are
checked against their actual source. Generic package trees are unchanged at
24 / 21 reports. `complete-compiler-review-v1.json` retains every full report,
source hashes and rejected report mutations. No analyzer filter is installed.

The inspected dependency posts one receive for each nonself peer, processes the
local slot first, and calls Waitany exactly once per nonself receive. An active
unconsumed request therefore remains at each wait under the valid-transport
contract. The runtime 1/2/4/8-rank matrices exercise those paths. The compiler
acceptance gate is corrected to full baseline/source comparison; it does not
claim these reports disappear.

V2 integrates the identical helper body and regression into accepted constructor
revision `00f7922`. The CI MPI step now runs the complete extension test entry,
including both constructor and diagnostic tests. All candidate/source identities
are recorded in `assembly-v2.json`; the original V1 source and failed harness
remain unchanged. Combined runtime, allocation and compiler comparisons run
again before promotion.

## Integrated candidate verification

The V2 source freeze covers 307 source/test/workflow/harness files. The helper
body and permanent diagnostic tests are byte-identical to V1; constructor code
is unchanged from accepted `00f7922`. Source identities are in `assembly-v2.json`.

| Gate | Checked result | Evidence |
| --- | --- | --- |
| Complete extension | All 1/2/4/8-rank commands pass on both Julia 1.12.7 and 1.13.2; 8,955 assertions per runtime across ranks | `diagnostic-v2-commands.json` and rank TOML files |
| Compiler/Aqua | Nine canonical Aqua assertions per runtime; all twelve complete concrete report trees equal V1 after checked source-line and exact keyword-owner maps; generic 24/21 trees unchanged | `complete-compiler-review-v2.json` and complete raw logs |
| Actual CI step | Updated MPI step passes in a fresh environment on minimum Julia with two real ranks | `ci-distributed-fft-v2-command.json` |
| Formatter | Six assertions pass: three Julia ASTs preserved and whole-tree idempotence | `format-v2-command.json`, `format-v2.log` |
| Complete-operation resources | 540 rank/profile pairs per runtime; all 2,700 allocation, returned-error and restored-input comparisons equal the reference exactly | `RESOURCE_CONTRACT_V2.md`, `resource-review-v2.json` |
| Independent boundary pass | 288 assertions pass per runtime for singleton/thin axes, zero fields, all six transform families and both output-call forms | `edge-v2-*` |

The compiler comparison checks the entire report, including all argument types,
frames, signatures and messages. The only generated identifier map is the same
byte-identical keyword helper's local index after the separately accepted
constructor introduces closures. Changed/missing frame controls reject. No
report is removed, and the V1 zero-report assumption failure remains retained
as historical evidence with its verified dependency classification above.

The resource reference uses the accepted public helper for out-of-place plans.
For the broken parent in-place helper, it uses the complete valid manual
snapshot/forward/inverse/error-reduction sequence. This is explicitly a reference
operation, not a successful execution of the broken public call. No latency or
multi-node scaling claim is made. This diagnostic does not close CAP-003's
operator, halo, coupled evolution or restart stages.

All local gates pass. The largest allocating diagnostic sample is 5,914,544 bytes,
identical to its complete reference operation; this utility is not advertised as
allocation-free. Comparison controls reject changed allocation, returned-error
and restored-state samples. The final independent boundary/source pass found no
additional confirmed defect within this diagnostic repair. Exact-revision CI
and its full compiler-log review remain required.
