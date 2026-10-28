# Distributed FFT collective construction

## Pre-edit brief

- Mode: repair configuration agreement in the existing optional FFT owner;
  establish the first CAP-003A/B prerequisite.
- Requirement: every rank supplies matching valid global dimensions and
  semantically equal expanded transforms before any Pencil communicator or
  plan is created. A local recoverable configuration/planning failure reports
  failure on all peers. Supported valid transform mathematics and layouts stay
  unchanged. No full global field is allocated by the wrapper.
- VERIFIED baseline: `constructor-mismatch-v1.jl` on real two-rank MPI accepts
  differing `(8,6,4)` / `(10,6,4)` global dimensions and FFT / BFFT directions
  on both Julia 1.12.7 and 1.13.2. Each rank records 2 passes / 2 failed rejection
  assertions. No incompatible forward transform was executed; no hang is
  claimed. Commands, per-rank parameters and source hashes remain in
  `2026-10-26/distributed-fft-contract/baseline-v1-*`.
- Canonical owner: `ext/HybridPlasmaPICPencilFFTSExt.jl`; reuse core
  `_mpi_collective_stage` from `src/parallel/mpi.jl:548` for local-only phases.
  PencilFFTs/PencilArrays continue to own topology, plan/storage and transport.
  Core fallback declarations and public exports retain their existing names.
- Contract: normalize one transform with the dependency's `expand_dims`;
  accept a three-transform tuple; validate shape, transform output dimensions
  and real/complex transitions using dependency metadata. Compare the complete
  normalized configuration by value against rank 0, including transform fields;
  do not compare hashes or display strings. Broadcast/transport failures and
  invalid/different communicators remain outside recoverable preflight.
  Custom transforms must be serializable and define value equality when their
  default equality does not describe configuration.
- Layout: the dependency's default 2D process grid uses input axes (2,3),
  intermediate(1,3), output(1,2). Reject empty rank partitions using original
  and transform-output lengths before topology construction. No replicated
  grid or custom transpose. Installed source `plans.jl::_make_intermediate_pencil`
  and `Pencils.jl::check_empty_dimension` establish these rules.
- Compatibility: matching FFT, BFFT, RFFT, inverse-real, real-to-real, identity
  and valid in-place combinations keep their dependency behavior. Equivalent
  single/expanded tuple forms agree. Invalid Int conversions become descriptive
  ArgumentError via the existing checked-integer owner. Invalid shape, transform
  and periodicity on one rank are caught inside collective local validation.
  Remote ranks receive the existing explicit collective-stage ErrorException;
  the origin retains its original exception and stack. No fallback success.
- Strategy: validate/normalize locally; collective status; broadcast rank 0
  configuration; local equality validation and collective status; construct
  topology only after agreement; wrap the dependency's noncollective plan
  construction in the same canonical stage. Shared transport helper unchanged.
- Verification: original two-rank witnesses, one-rank invalid configuration,
  valid mismatches/equivalent representations, parity and roundtrip with known
  serial transform references, DC/Nyquist/asymmetric shapes,2/4/8 real ranks on
  minimum/latest Julia; original extension/coreMPI neighbors, Aqua and concrete
  plus full compiler trees. Record constructor cost and warmed transform
  allocation/output differentials separately. No latency/scaling claim without
  the predeclared performance gate. Update docs/README/tracker and push this
  slice separately, then require all six exact-SHA CI jobs.
- Wrong patches: validating only rank 0; comparing dimensions while ignoring
  transform direction/fields; making failures into default FFT plans; applying
  rank-local checks after topology creation; comparing raw serialized padding.
  Adversarial cases and same-value expanded representations must reject these.

## Subsequent CAP-003 work

- [x] Accept the constructor slice: `00f7922` has six green exact-SHA jobs and complete compiler review, as recorded below.
- [x] Land rank-local storage and the derivative/vector/filter operator slices;
  [canonical stage and acceptance record](DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md)
  tracks `4fafe4a` and `4724947`, including their six-job/compiler acceptance.
- [ ] Compose local deposition, halos, migration, fields and hybrid timesteps.
- [ ] Add per-family restart, diagnostics, ES/EM and measured multi-node scaling.

The constructor and later operator slices have separate acceptance records.
Coupled distributed evolution, restart, cluster performance and complete
simulator scope remain open.
Primary layout contract: https://jipolanco.github.io/PencilFFTs.jl/dev/tutorial/
and https://jipolanco.github.io/PencilArrays.jl/dev/PencilArrays/.

## First complete campaign finding

The original two-rank size/direction witnesses pass on both candidate runtimes.
The broader first campaign has 257 passes / 1 failure at one rank: a heterogeneous
integer tuple recursively re-enters the generic shape method and raises
StackOverflowError. The two-rank stage was interrupted through its own launched
tool session after the same failure, with all logs/source retained. No completed
multirank pass is claimed for that stage.

The checked source keeps the inherited `_checked_fft_shape(n::Tuple)` self-call;
the NTuple overload does not match mixed BigInt/Int entries. V2 repairs the
canonical Tuple owner directly with checked per-component conversion. The
unchanged-parent dispatch witnesses retain StackOverflowError on both runtimes.
Valid heterogeneous and wider-integer tuples now run alongside the invalid
overflow case; no fatal exception catch or dimension clamping was introduced.

## Byte-transport follow-up (before editing)

The V1 concrete compiler campaign reports two possible undefined-buffer errors
per transform on both runtimes, entirely in `MPI.bcast` at dependency
`src/collective.jl:99,101`. Inspection shows the root/nonroot complementary
branches initialize that buffer; this is not a demonstrated runtime failure.
However, the wrapper already serializes its configuration inside the collective
validation stage. Calling the arbitrary-object broadcast serializes those bytes
again outside that stage and uses its Cint message-size conversion there.

Use the dependency's typed `Bcast` / `Bcast!` primitives directly: validate the
encoded byte length against the MPI Cint count range in the existing local
validation stage, broadcast that length, allocate the receiver buffer, then
broadcast the original encoded bytes. Keep one serializer and the same complete
value comparison. No analyzer filter, suppression, exception fallback, or change
to core checkpoint transport. Memory exhaustion and MPI transport failure remain
outside the recoverable local-error contract. Re-run the complete original
matrix, concrete/compiler reports and canonical Aqua settings after this change.
The V1 analysis harness stopped after its failing concrete testset; V2 must use
an outer testset and save all reports before assertion completion. Also use the
repository's exact Aqua settings (ambiguity checking enabled).

## Final constructor candidate verification

The V3 source is frozen in `candidate-v3-source.json` (394 source/harness files)
and `candidate-v3-source.tar.gz` under the evidence directory above. All frozen
hashes still match after the completed transform campaigns.

| Gate | Observed result | Evidence |
| --- | --- | --- |
| Original size/direction disagreement | Both runtimes reject collectively; 4 assertions per rank | `candidate-original-v1-*` |
| Full transform/error matrix | 1/2/4/8 ranks pass on both runtimes; 7,530 assertions per runtime across all ranks | `contract-v3-commands.json`, per-rank TOML files |
| Concrete compiler and Aqua | 4 concrete signatures have no reports; all 13 assertions pass per runtime with canonical Aqua settings | `quality-v2-*` |
| Full generic compiler trees with FFT extension loaded | All 24 / 21 report trees match the parent after only worktree-root normalization; missing report, changed frame and changed message mutations are rejected | `complete-compiler-review-v2.json` |
| Warmed transforms | 135 rank/profile pairs per runtime; 1,350 allocation comparisons and 1,350 output-hash comparisons all equal the parent | `resource-differential-review-v1.json` |
| Workflow integration | The actual new CI step passes locally on minimum Julia with two real ranks in a fresh temporary test environment | `ci-distributed-fft-v3-*` |
| Formatter | 6 assertions: three owned Julia ASTs preserved and two whole-tree idempotence passes | `format-v3-*` |
| MPI/checkpoint/extension neighbors | 47,105 assertions pass on each runtime | `neighbors-v3-*` |

The largest observed warmed transform allocation is 8,080 bytes in each direction,
identical to the parent. This is not a zero-allocation or latency claim. Constructor
allocation observations are recorded separately (nine profiles, five samples
each): the new preflight adds 4,432–4,840 bytes on Julia 1.12.7 and
3,600–4,008 bytes on Julia 1.13.2; configuration agreement adds
three status reductions and two byte-count/payload broadcasts before a plan is
returned. No constructor latency or multi-node scaling claim is made.

A subsequent independent diagnostic pass found a separate inherited defect:
`distributed_fft_roundtrip_error` calls `parent` on an in-place multi-view
container. FFT!, RFFT! and R2R! each fail with MethodError on both parent runtimes.
The source/logs remain in `roundtrip-inplace-baseline-v1-*`; its separate
[repair plan](DISTRIBUTED_FFT_DIAGNOSTIC_PLAN.md) preserves this constructor
candidate unchanged.

All local constructor acceptance gates are complete, and the repair is pushed as
`00f7922`. [Exact-revision GitHub CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40220039958)
passes all six jobs (checked 2026-10-29). The archived complete compiler review
records 797,946 / 797,947 assertions on actual Julia 1.12.7 / 1.13.2, with all
24 / 21-report forests equal after checked unchanged-source line mappings. This repair
closes the constructor prerequisite only. Subsequent operator progress belongs
to the linked canonical plan; coupled evolution, restart and cluster scaling
remain open.
