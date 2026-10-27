# Archive payload inference repair

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: JET-002 implemented and focused gates passed on both native versions;
pushed as `6341c4b`. Full-package integration runs are in progress.

## Pre-edit brief

- **Mode:** canonical inference repair with unchanged archive behavior.
- **Acceptance:** the concrete archive loader and unchanged package compiler
  gate pass with HDF5 and Pencil extensions loaded on Julia 1.12.7 and 1.13.1;
  valid archives, checksum-before-decode ordering and existing errors remain.
- **Baseline:** full 1.12 suite at `fe519e9` reports 31 findings against ceiling
  30. A retained test environment reproduces 31 without numerical tests;
  nine concrete `load_run(String)` reports flow through unrelated HDF5 reads.
- **Cause — VERIFIED:** optimized compiler IR for `load_run(String)` retains
  repeated `container.state` reads and `IOBuffer(container.state)` as `Any`,
  despite the earlier `Vector{UInt8}` check. Binding the field once before that
  check gives `sha256(::Vector{UInt8})`, `IOBuffer`, and
  `Serializer{IOBuffer}`. The otherwise identical isolated candidate has zero
  concrete loader reports and 22 package reports on 1.12. Evidence:
  `/tmp/hybrid-jet002-payload-{baseline,v1}-ir112.log` and
  `/tmp/hybrid-jet002-payload-v1-all112.log`.
- **Contract:** validate the container, metadata, checksum type and byte-vector
  payload; hash the stored bytes before decoding them. Preserve accepted states,
  aliases, return layout, exception types and validation order.
- **Owner and surfaces:** `load_run` in `src/io/metadata.jl`; `load_archive`
  delegates to it. Metadata, reproducibility, I/O and compiler tests cover callers.
- **Compatibility:** no API, schema, dependency or serialized-format change.
  This addresses an inference defect; the reports do not prove that a valid
  archive executes the unrelated HDF5 failure paths at runtime.
- **Strategy:** retain one local payload binding through its existing validation,
  SHA-256 calculation and nested deserialization. Add concrete compiler
  regressions for String and SubString paths, with no suppression or gate change.
- **Verification:** same-environment before/after compiler reports, native
  1.12/1.13 metadata/reproducibility/I/O regressions, public string variants,
  corruption and malformed-container failures, formatter and exact final diff.
  Reject three superficially green patches: removing optional imports (retain
  both extensions), raising the report ceiling (keep 30), or decoding before
  checksum verification (retain corrupt-byte rejection checks).

## Environment controls

In fresh processes using the exact retained 1.12 test environment, base and
HDF5-only modes each have 0 concrete loader / 26 package reports; Pencil-only
has 0 / 22; all optional test imports have 9 / 31. These observations do not
attribute a runtime defect to either dependency. Logs:
`/tmp/hybrid-jet002-retained-{base,hdf5,pencils,all}112.log`.

The current stable Julia release was rechecked on 2026-10-15 against the
[official downloads page](https://julialang.org/downloads/manual-downloads/):
1.13.1. CI retains both minimum minor `1.12` and moving stable `1`.

## Verification and remaining scope

Final source SHA256 for `src/io/metadata.jl`:
`929987711d2a7f808818a940703a59bd39156aa2c0ea898219e9f6feb8697e84`.
Frozen source/test hashes: `/tmp/hybrid-jet002-final-hashes.json`.

- Final compiler test passes all three assertions on each version: package
  reports are 22 on 1.12.7 and 21 on 1.13.1 (ceiling remains 30), with zero
  concrete loader reports for String and SubString paths. The final gate loads
  HDF5, PencilArrays, PencilFFTs, Statistics and Aqua first. Logs:
  `/tmp/hybrid-jet002-final-gate{112,113}.log`; orchestrator
  `/tmp/hybrid-jet002-final-gates.py` exits 0.
- Comparing current-source 1.12 baseline and candidate report leaves removes
  exactly the nine concrete loader findings and adds none. The other 22 leaf
  messages are unchanged. Driver/log:
  `/tmp/hybrid-jet002-removed-leaves.{py,log}`. This is not a classification or
  dismissal of the remaining package reports.
- Existing metadata, reproducibility and extra-I/O files pass 234 assertions
  on each version using the identical final production bytes. Logs:
  `/tmp/hybrid-jet002-payload-v1-runtime{112,113}.log`.
- Final metadata tests pass 80 assertions on each version, including the 60
  already counted above and 20 new checks for String/SubString paths/checksums,
  exact scalar values, shared cyclic payloads, metadata and unchanged stored
  bytes. Existing checks cover malformed containers/payloads, bad checksums,
  schema mismatch, corrupted bytes before decoding and separate-process reads.
- JuliaFormatter 1.0.62 is idempotent for all three changed Julia files at
  margin 100; `/tmp/hybrid-jet002-final-format.log` exits 0. The production diff
  only binds and reuses the already validated payload; serialization, error
  messages/order, schemas, exports, dependencies and report ceiling are unchanged.

Full-package runs now use immutable `6341c4b` on both native versions, including
the later checkpoint and continuity changes. Logs:
`/tmp/hybrid-full-6341c4b-julia{112,113}.log`; initial source/environment hashes:
`/tmp/hybrid-full-6341c4b-manifest.json`. Completion remains pending.
