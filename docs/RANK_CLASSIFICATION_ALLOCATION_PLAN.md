# Rank-classification allocation reduction

Status: pushed as `b3e2eb62e644ab7f2de970af56bbc1939202fbc8` after local
acceptance on Julia 1.12.7 and 1.13.2.
[Exact-commit GitHub CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39432448426)
passes all six jobs. Downloaded logs verify Julia 1.12.7/1.13.2 and
504,062/504,063 package assertions. Complete CI compiler trees retain 23/22
reports with no added or removed groups. The six log hashes and exact SHA
are recorded in `/tmp/hybrid-ci-b3e2eb6-20261020/verified-results.json`.

The pinned-main probe at `/tmp/hybrid-rank-allocation-review-20261020/probe.jl`
passes six value controls on Julia 1.12.7. Warm `rank_index` and
`rank_of_position` calls allocate 320/640/960 bytes in 1D/2D/3D for both Float32
and Float64 grids. The migration resource matrix independently records roughly
1.28 GB allocated while routing 1,048,576 particles on minimum Julia; that
whole-call number is not attributed solely to coordinate labels.

`rank_index` eagerly constructs `"coords[$d]"` before calling
`_rank_coordinate`. The same indexed labels appear in both pencil-coordinate
methods and both pencil-owner methods. MPI neighbor validation and the scalar
partition helper call the same validator with constant labels. This is the
complete caller set from `rg -n '_rank_coordinate' src test`.

## Pre-edit brief

- Mode: measured optimization with unchanged numerical and validation behavior.
- Requirement: remove ordinary indexed-label allocation from rank and pencil
  classification; retain public values, errors, integer conversion, periodic
  wrapping and complete input validation order.
- Baseline: the warm public probes above and the existing migration resource
  matrix. Eager labels are the source-level candidate cause; an isolated change
  and matching allocation measurements must establish the reduction.
- Canonical owner: `_rank_coordinate` in `domain_decomposition.jl`, called by
  logical rank, pencil decomposition and MPI neighbor validation.
- Compatibility: keep existing constant-label callers and public signatures.
  Supply an optional private coordinate index to the same validator; construct
  the indexed error label only after the range check fails. No duplicate
  coordinate validator or public mode is introduced.
- Strategy: change the indexed callers together, preserving their label text,
  and leave all classification arithmetic intact.
- Verification: exact ordinary/boundary values and error text, large signed and
  unsigned integer rejection, all periodic-axis combinations, independent
  integer rank oracles, pencil ownership, existing migration/pencil/MPI tests,
  zero warm allocation for ordinary classification, complete compiler-tree
  comparisons and formatting on both supported Julia versions. Measure complete
  migration allocation before claiming a whole-call improvement. Push this
  change separately from particle storage repair and require all six exact-SHA
  CI jobs.

Plausible incomplete patches include dropping integer range validation,
checking only the first coordinate, or changing the order of periodic wrapping
and conversion. Large-coordinate/error-order tests and unchanged arithmetic
syntax checks must reject those alternatives.

## First candidate and orientation follow-up

Both baseline runtimes pass 1,284 value/error assertions and fail all 15 ordinary
allocation checks. Candidate V1 preserves 46,133 checks on minimum Julia and
removes allocation from logical rank/position lookup and pencil rank lookup.
Three allocation checks remain failed: pencil owner lookup with a runtime
Symbol orientation still allocates 320 bytes. V1 is not accepted; later native
gates did not run. Its 347 formatting and 92 unchanged-arithmetic checks pass.

**VERIFIED localization:**
`/tmp/hybrid-pencil-dispatch-review-20261020/probe.jl` measures 320 bytes in the
checked pencil-owner path for each Symbol orientation, versus zero for the
corresponding `Val`. Orientation conversion, axes lookup and scalar partition
lookup separately allocate zero. `code_warntype` in the same log shows `Any`
for the axes result and split-axis indexing after the runtime `Val` conversion.
The six value controls pass. This evidence identifies a remaining dispatch cost;
it is not a failure of the lazy coordinate labels.

The follow-up preserves the same normalizer and its Symbol/Val inputs. Return
one of the three concrete orientation values for supported Symbols and share
the existing orientation-error message with the unsupported-Val handler. Keep
all classification arithmetic, public validation order and original allocation
assertions. Add Symbol/Val parity and invalid-orientation message/order checks
for bounds, ownership and orientation queries. Re-run the final candidate on
both runtimes; do not accept the three failures by raising the allocation limit.

## Final candidate evidence

V2 at `/tmp/hybrid-rank-allocation-v2-20261020` passes 46,188 focused assertions,
16 concrete compiler checks and 137 MPI/quality integration assertions on each
of Julia 1.12.7 and 1.13.2. The 52 added orientation value/error-order checks
also pass against the original code on both runtimes. Formatting passes 347
checks and the unchanged-arithmetic audit passes 93. Complete package compiler
trees match the baseline: 22/21 reports, 14/13 complete groups and no added or
removed groups. `compiler-comparison-results.json` retains both commands and
log hashes. All native logs and the frozen source manifest were checked before
promotion. Only the two classification source files, the new regression test
and its `runtests.jl` include differ from the accepted source/test baseline.

**VERIFIED allocation reduction:** final-source measurements at
`/tmp/hybrid-rank-allocation-resources-v2-20261020/summary.json` compare V2
with the retained, hash-checked pinned-main baseline logs in
`/tmp/hybrid-rank-allocation-resources-20261020`. Each command passes 175
assertions. All 35 complete output fingerprints per runtime match, with zero
compilation/recompilation in the measured calls. The seven profiles cover
Float32/64, 1D–3D, 4,096 particles per rank and a 524,288-particle 3D case.
Allocated bytes fall by 59.68–84.19% on minimum Julia and 47.04–76.16% on latest
Julia. The largest case changes from 641,606,992/440,280,400 bytes to
138,290,512 bytes on the two runtimes. All five samples per profile are retained.
These unpaired process observations establish allocation reduction, not an
elapsed-time speedup.
