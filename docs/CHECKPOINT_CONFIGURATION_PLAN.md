# Checkpoint configuration ownership and publication

Status: configuration ownership repair verified on Julia 1.12.7 and 1.13.1;
independent serialization and publication gaps remain open.

## Configuration links to solver state

VERIFIED: captured model/shape configuration can refer to particle fields,
solver scratch or clock references. Restoration copies saved arrays into the
destination but installs deserialized configuration; its links can then point
to detached copies. The original evolution and restarted evolution diverge.
Read-only grid geometry and independent parameter arrays must remain usable.

The canonical state validator now inspects configuration and both source and
payload ownership. Its bounded object walker follows stored fields, logical
collection contents and custom wrapper metadata. Declared array parents are
visited for metadata without expanding a view's logical extent. Later logical
visits upgrade prior metadata-only visits; cycles visit each object at most
once per mode. Standard Array/Memory buffers and Dict/IdDict capacity stay opaque.
Custom arrays must accurately implement Base's parent and alias interfaces.

The owned roots include particle/field arrays, solver scratch, clock/phase
references and FFT plans. Capturing an owned FFT plan is also invalid: on both
Julia versions a public save/load probe turns its live native handle into
`C_NULL`, while the destination grid plan remains valid. Those handles were
inspected, never executed. Logs:
`/tmp/hybrid-config-ownership-grid-plan112.log` and
`/tmp/hybrid-config-ownership-grid-plan113.log` (16 checks per version).

The candidate includes forward, scaled inverse and inner inverse plan objects,
while retaining accepted read-only `kfull`, `kvec`, `midx` and `ik` references.
These checks target detachment before save/archive publication or restore
mutation. The final candidate retains Set-to-dictionary identity and registers
packed-array chunks plus Dict/IdDict backing allocations as owned, without
traversing unused capacity. Equal immutable metadata and disjoint logical views
remain accepted. Focused, neighboring, compiler and real MPI gates pass on
both Julia versions. No checkpoint schema or global Serialization method is changed.

Final acceptance includes Hybrid/CAM-CL scalar and collection layouts, both
EM APIs with fixed/mobile ions, Hall-MHD, Float32/64, 1D–3D, independent mutable
parameters, cyclic structures, disjoint views, malformed later species,
retained payload links, custom array/dict/set/string wrappers and SimpleVector.
The checks below cover these contracts, including MPI's collective preflight.

## Ownership repair verification

The final v11 ownership matrix passes **8,354 assertions on each of Julia
1.12.7 and 1.13.1**, across all 12 testsets. Source SHA256:
`006e213da683b3bc68b83318bb454eb5f74adf56969bb1c6c35f92df317b41df`.
Logs: `/tmp/hybrid-config-ownership-v11-focused112.log` and
`/tmp/hybrid-config-ownership-v11-focused113.log`.

Sixteen concrete validator JET calls per version have zero reports. Whole-package
analysis in the isolated `d2c3da7` plus IO candidate snapshot retains 24/22
reports against the unchanged ceiling of 30; later non-IO fixes are excluded
from those counts. Logs: `/tmp/hybrid-config-ownership-v11-jet112.log` and
`/tmp/hybrid-config-ownership-v11-jet113.log`.

Real two-rank MPI integration of this IO source over `8e7d5be` passes 173
assertions on rank zero and 165 on its peer on each Julia version. Coverage
includes scalar/two-species layouts, Float32/64, rejected rank-local links,
previous-generation preservation, collective incoming rejection before
mutation, independent cycles and exact restart continuation. Logs:
`/tmp/hybrid-ownership-mpi-v11-v14fixture112.log` and
`/tmp/hybrid-ownership-mpi-v11-v14fixture113.log`.
Earlier temporary harness failures used the wrong peer exception type, the
wrong Hybrid payload field path, and one directory for incompatible precisions;
they are retained as test-fixture failures, not production findings.

Measured ownership guards allocate 0 B for 18 standard configurations;
mutable-cycle cases allocate 32,544–55,792 B. A primitive parameter vector of
length 1 or 1,000,000 costs the same 42,144 B. Packed backing lengths 6 and
1,000,000 cost 44,944 and 49,104 B. No timing or speedup claim follows.
Log: `/tmp/hybrid-config-ownership-v11-resources112.log`.

The preceding v10 snapshot passes 15,100 adjacent I/O/reproducibility/MPI
assertions on each Julia version (terminal exit 0):
`/tmp/hybrid-config-ownership-v10-final112.log` and
`/tmp/hybrid-config-ownership-v10-final113.log`.
The v11 addition changes only dictionary backing registration and is covered
by 560 new assertions per version plus the complete final ownership matrices.
Final full-package integration is a separate gate after the outstanding slices.

## Separate confirmed follow-ups

- Internal parameter views: Julia's serialization can detach an Array-backed
  SubArray from its separately stored parent. Independent configuration remains
  vulnerable even after solver-owned links are rejected. Dict/backing-buffer
  relationships can detach too: both native-version public probes pass four
  witness checks (`/tmp/hybrid-config-internal-dict-backing112.log` and
  `/tmp/hybrid-config-internal-dict-backing113.log`). Plan a canonical
  configuration representation or reject unsupported internal relationships
  explicitly; verify mutation-dependent continuation across checkpoint,
  archive and MPI entry points. Do not override global serialization behavior.
- Native resources independent of the solver: arbitrary user configuration can
  contain handles that Serialization cannot reconstruct. Establish an explicit
  supported configuration contract. An already-detached legacy payload cannot
  recover the original ownership relationship from its serialized bytes.
- Serial checkpoint publication: `_save_checkpoint` writes directly to the
  destination. An injected custom serialization error changes an existing
  checkpoint from 862 to 411 bytes and leaves it unreadable. Both-version probes
  pass six witness checks each; the archive remains unchanged for that same
  serialization failure because its payload is serialized before publication.
  Logs: `/tmp/hybrid-checkpoint-publication-probe112.log` and
  `/tmp/hybrid-checkpoint-publication-probe113.log`.
  Plan same-directory staged writes and atomic publication, preserving the
  original error, cleaning only owned temporary files and retaining the previous
  checkpoint on pre-publication failure. Verify successful replacement,
  serializer/write/close/rename failures, retry, path ownership and bounded
  resources. MPI's existing generation publication is a separate owner.

The serial publication follow-up now has a [dedicated implementation plan](SERIAL_PUBLICATION_PLAN.md),
including path, permission, asynchronous snapshot and failure-cleanup contracts.

## IO-010 serialization-contract implementation brief

Reviewed before production edits; final gates remain pending.

1. **Mode:** repair canonical configuration validation.
2. **Requirement:** reject relationships and resources that cannot survive
   serialization before publication or restore mutation. Preserve independent
   mutable parameters, shared Arrays/Refs/cycles, ordinary Dict/IdDict/Set,
   standalone views and disjoint logical views. Validate incoming retained
   graphs too; already-detached legacy bytes cannot reveal their lost links.
3. **Evidence:** public checkpoint/archive probes pass 36 assertions per native
   version (`/tmp/hybrid-io010-recon-v2-{112,113}.log`). Repeated, overlapping and
   parent-linked Array-backed views, separately held mutable indices, dictionary
   backing links and native FFT plans lose relationships or handles. Scalar/Ref
   pointers become null; pointer arrays retain addresses. No resource was executed.
4. **Verified cause:** Julia's `Serialization.jl:193` nulls pointers, lines
   287–291 copy Array-backed SubArrays, and lines 378–388 rebuild dictionaries.
   The existing validator checks solver ownership, not internal reconstruction.
5. **Contract:** reject copied view storage or mutable indices exposed elsewhere,
   exposed dictionary backing storage, native FFT plans and non-null pointers,
   including nested arrays. Preserve safe value codecs such as Regex and BigInt.
   Custom serializers own their reconstruction semantics; validation does not
   establish safety for arbitrary integer-encoded handles or external callbacks.
6. **Owner:** `_validate_checkpoint_configuration` in `src/io/checkpoint.jl`,
   shared by Hybrid, CAM-CL, EM and Hall configuration capture/restore, scalar and
   collection checkpoints, archives and MPI preflight. ES has no custom graph.
   Generic `save_run` and `async_save` retain generic Serialization behavior.
7. **Compatibility:** no schema or global Serialization changes. Previously
   accepted graphs with demonstrated loss now fail explicitly. Existing mutable
   configuration/type restoration and solver-owned-storage rejection remain.
8. **Strategy:** extend the bounded walker with edge/storage ownership. Separate
   stored edges from synthetic iteration; inspect repeated views before identity
   deduplication. Preserve logical-versus-parent metadata traversal, opaque spare
   capacity, primitive fast paths and cycle bounds. Reuse overlap predicates and
   explicit safe-codec owners without serializing during validation.
9. **Verification:** both native versions; all configurable families, layouts,
   Float32/64 and 1D–3D; save/archive/incoming rollback and mutation-dependent
   continuation. Real MPI must reject a bad graph on one rank collectively while
   preserving the previous generation and destinations. Check resources,
   concrete/package JET, formatter and neighboring IO/repro/MPI tests. Disjoint
   view positives reject a blanket ban; repeated-view negatives reject early
   deduplication; unused-capacity positives reject indiscriminate reflection.

Before implementation, this contract was extended to distinct Array headers
sharing storage. A vector and `reshape(vector, 2, 2)` share storage before
serialization but detach afterward. Both-version public probes pass 42 assertions
in `/tmp/hybrid-io010-recon-v3-{112,113}.log`. Preserve repeated references to one
Array; reject overlapping storage represented by distinct serialization owners.

The initial candidate passes 33 Julia 1.12 contract assertions, including safe
numeric/regex values, null pointers, disjoint views, shared arrays, cycles,
read-only grid geometry and immutable-versus-mutable wrapper edges:
`/tmp/hybrid-io010-contract-v1-112.log`. This is intermediate evidence, not a final
verification result.
