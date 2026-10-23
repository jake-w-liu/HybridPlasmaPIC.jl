# Backend memory-status constructor consistency

Status: pushed as `c810eb0`; both native correctness/compiler and fixed cost
lanes pass. [Exact-SHA CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40005434857)
remains required. Parent:
`ba3694c` (production GPU owner unchanged from `61356f1`). Evidence: `../HybridPlasmaPIC-audit-evidence/2026-10-23/backend-constructor`.

## Contract and scope

The exported `BackendMemoryStatus` represents byte counts. Its documentation
uses `nothing` for unavailable counters, and its existing three-positional
constructor requires each counter to be a nonnegative integer that fits `Int`.
All current backend telemetry callers use that checked constructor.

**VERIFIED baseline:** the parent's keyword and full-positional constructors
bypass those checks.
The retained public witness has 185 passing controls and ten failed negative
rejections on each supported Julia version. All five counters admit `-1`
through either bypass. Source: `src/parallel/gpu.jl:20–61`, checked against
`2993f5a`; its production bytes are unchanged at this plan's parent.

Compatibility decision: all public construction forms will enforce the existing
checked-constructor policy. This intentionally rejects negative, noninteger and
out-of-range counts previously admitted or implicitly converted by the raw
constructors. Valid integer conversions, `nothing`, defaults and metadata
conversions remain supported. Cross-counter ordering is not newly constrained:
independent telemetry samples can be inconsistent, and `memory_pressure`
already clamps the resulting ratio.

The canonical owner remains `_nonnegative_bytes`. The full-positional inner
constructor calls it for all five counters; outer wrappers delegate once.
This prevents generated keyword construction from bypassing validation and
avoids maintaining two copies of the policy. Keep the public fields, exports,
backend telemetry APIs and pressure calculation unchanged.

## Acceptance tracker

- [x] Review the exported type, all constructors, telemetry callers and consumer.
- [x] Retain the original failed public witness and valid controls.
- [x] Recheck the witness against the exact parent and capture raw-constructor
  metadata conversions before editing production code.
- [x] Implement canonical validation and document the intentional compatibility
  tightening; no sentinel values, silent coercion or caught validation errors.
- [x] Cover all three construction forms and five counters: missing, zero,
  positive and maximum values; signed/unsigned wide integers; fitting BigInt;
  negatives, overflow, fractional/integral floating values, NaN and infinity.
  Preserve defaults, note conversion and valid pressure behavior.
- [x] Run the complete backend helper suite on Julia 1.12.7 and 1.13.2, concrete
  constructor/compiler checks, Aqua and a complete package compiler comparison.
- [x] Compare valid complete status values and warmed allocation against the
  parent; require no additional allocation for ordinary integer/missing inputs.
  Record construction cost without claiming an optimization. Existing checked
  wrapper and real telemetry paths must retain their behavior.
- [x] Review the final diff and composed source, verify held-change hashes and
  commit/push separately as `c810eb0`; all 364 committed Julia files match the
  validated candidate. Fifteen final command logs and 144 Markdown links checked.
- [ ] Require all six exact-SHA CI jobs green.

This is a host metadata invariant repair. It does not establish GPU kernel
correctness, device memory capacity or a complete GPU-resident solver.

Fresh parent evidence: `fresh-parent-commands.json` retains both native commands
and log hashes. Both runtimes reproduce 185 passes / ten negative-rejection
failures. The 21 metadata cases per runtime preserve `SubString` notes and raw
keyword/full-positional `0`/`1` Boolean conversions. Invalid metadata exception
classes are recorded for comparison. No production edit preceded those runs.

The candidate moves all five checks into the full-positional inner constructor;
the keyword and three-positional forms delegate to it. New permanent tests
exercise independent counter samples and pressure clamping as well as rejection.
Both runtimes pass 1,313 permanent tests, the original 195-assertion witness,
21 metadata comparisons, 116 backend helper checks (including actual Metal),
49 concrete compiler checks and Aqua. Full package trees match the matched
parent after source-root replacement only (29 / 27 reports). Every one of 27
complete-value/allocation profiles retains zero allocation in ten paired samples.

The external compiler fixture originally required zero reports for deliberately
invalid floating counters. Its retained failures and correction are documented
in `concrete-review.md`: valid signatures remain report-free, while each invalid
signature requires the precise uncaught `ArgumentError` and runtime rejection.
No production diagnostic is suppressed.

The fixed 27-profile timing matrix passes on each runtime. The existing checked
wrapper's maximum ratio of medians / paired-median ratio is 1.046512 / 1.029693
on 1.12.7 and 1.038560 / 1.032827 on 1.13.2, below the registered 1.5 limit.
All complete outputs and zero-byte allocation samples match. Newly validated
raw keyword/full-positional forms have maximum median ratios of 3.716241 /
3.817937 and 3.175272 / 3.265602 respectively. Their new checks have a measured
cost; this correctness repair is not advertised as an optimization. All samples
remain in `timing-112.toml` and `timing-1132.toml`.

Final source review traced all telemetry callers through the unchanged checked
wrapper and pressure consumer. Both installed Julia versions' string conversion
owner delegates `convert(String, ::AbstractString)` to `String`; the captured
ordinary and SubString metadata controls retain their results. No GPU kernel,
telemetry fallback, field type or pressure formula changes in this repair.


## Exact-revision CI status

The repair is pushed as `c810eb0`. [CI run 40005434857](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40005434857)
has five passing jobs and one failed minimum-Julia job. Julia 1.13.2 passes
625,492 assertions; Julia 1.12.7 records 625,490 passes and one failed allocation
gate at `test/shock_total_flux/reuse_resources.jl:53`: 1,501,504 bytes exceed the
unchanged 1,492,110.4-byte bound. Both complete compiler trees are byte-identical
to the reviewed `ba3694c` trees (24 / 21 reports). Formatting and all MPI jobs pass.

The shock owner and resource test are byte-identical to green `ce52836`, but
the failed log lacks per-profile allocation attribution. The cause remains
unverified; no retry, relaxed gate or billing exception has been claimed.
The allocation investigation must be resolved before this revision's CI gate
can be marked complete. Full logs and hashes are retained in the persistent
`ci/40005434857/` and `shock-total-resource-ci/` evidence directories.


The documentation-only successor `d6bce16` also fails one minimum-Julia allocation
assertion (1,754,288 bytes against 1,733,265.6). Its complete compiler trees match
the same reviewed baseline. Both local runtimes pass all 12 profiles / 120 checks
with the original and actual CI coverage/bounds flags; maximum allocation ratio
is 1.014246 and every measured sample has zero compilation. The CI cause remains
unverified.

A separate test-only change records type/count/profile, bytes, compilation and
GC time when the existing gate fails. Parsed-code comparison proves that removing
this logging conditional exactly recovers the original fixture. Both runtime
lanes retain all 120 passes; replaying the exact state deliberately yields 108
passes / 12 failed allocation gates and 12 complete diagnostic records. This
instrumentation improves the next CI failure's evidence; it does not close the
allocation investigation or relax any gate.


The later documentation-only `f850b84` run also finishes with five successful
jobs and a failed minimum-Julia allocation check: 1,736,112 bytes exceed the
unchanged 1,733,300.8 limit. It records 625,490 passes / one failure. Both complete
compiler trees exactly match reviewed `ba3694c` (24/21 reports). The full archive,
log hashes, complete trees and zero-length differences are retained under
`ci/40017267108/`. This revision predates the added per-profile diagnostics.
The attribution remains open; no billing exception or retry is claimed.

The diagnostic-only successor `f1410b7` passes all six
[exact-SHA jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40027319048).
Complete logs retain 625,494 / 625,495 passes and full compiler trees equal to
the reviewed parent after checked line mappings (24 / 21 reports). There are no
failed samples in that run to attribute; the historical allocation cause remains
unresolved, and the failed original revisions remain recorded.
