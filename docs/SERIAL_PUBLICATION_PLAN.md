# Atomic publication of serialized outputs

Status: IO-011 is implemented and verified on Julia 1.12.7 and 1.13.1, following
the checkpoint configuration repair `3495c64`.

## Pre-edit brief

- Mode: repair the shared serial file-publication boundary.
- Requirement: serialization, flush, close or publication failure must preserve
  the previous destination and leave no owned staging file after successful cleanup.
- Baseline: an injected serializer failure truncates a checkpoint from 862 to
  411 bytes and makes it unreadable. Six witness assertions pass on each Julia
  version in `/tmp/hybrid-checkpoint-publication-probe112.log` and `113.log`.
- Cause: `_save_checkpoint` calls `serialize(path, state)`, opening the final file
  before serialization succeeds. `save_run` buffers its state but still writes
  the final container directly; `async_save` also publishes directly.
- Contract: fully serialize, flush and close a same-directory staging file before
  renaming it over the destination. Preserve original exceptions; expose cleanup
  failures as additional errors. Delete only the staging path created by this call.
- Owners: one private writer shared by `src/io/checkpoint.jl`, `metadata.jl` and
  `hdf5_output.jl`. `archive_run` reaches it through `save_run`. MPI generation
  publication keeps its separate collective lifecycle.
- Compatibility: retain bytes, schemas, signatures and caller return values.
  `async_save` still takes its deep copy before spawning, with write errors
  delivered through the task. Existing regular-file rwx bits and write-access
  checks are retained. New staging files use mode 0600. Resolve existing and
  dangling leaf symlinks to their physical referents, including directory-link
  paths containing `..`. Publication replaces the destination inode, so other
  hardlinks retain the old file. Concurrent-writer ordering and crash durability
  are outside this contract.
- Strategy: explicit staging ownership and close/rename/cleanup, without buffering
  an additional complete raw checkpoint. Native rename must not predelete the
  destination; the lifecycle must not swallow cleanup failures.
- Verification: serializer and partial-write errors, flush/close/rename failures,
  error aggregation, closed handles, exact old bytes, fresh-file failure/retry,
  unrelated-file preservation, symlink chains/cycles, permissions, raw bytes and
  return values. Run public checkpoint/archive/async integration, neighboring I/O,
  compiler and allocation checks on Julia 1.12.7 and 1.13.1.

## Path evidence

Both native runtimes pass 14 path/permission witness checks in
`/tmp/hybrid-serial-publication-path-contract112.log` and `113.log`.
Existing serialization follows symlinks and hardlinks, preserves mode 0640 and
refuses a non-writable mode-0440 file. Rename alone replaces mode-0440 files and
leaf symlinks, so it does not provide the required path contract by itself.

The temporary physical-parent resolver passes 12 assertions on each version:
`/tmp/hybrid-serial-leaf-resolution-probe112.log` and `113.log`.
It covers relative and dangling chains, directory symlinks with `..`, cycles,
missing parents and agreement with the actual serialization destination.
These prototype checks are supplemented by the final implementation gates below.

## Final verification

Both native versions pass **23,747 I/O, reproducibility and MPI assertions**
(17m13.6s on Julia 1.12.7; 20m41.4s on 1.13.1), with terminal exit zero:
`/tmp/hybrid-publication-v3-final112.log` and `113.log`.
The 293 focused publication checks also pass separately on each version.
They cover partial writes, serializer/flush/close/rename failures, a close that
fails before closing, error aggregation, replacement staging entries, symlink
chains and cycles, directory-link `..` semantics, rwx modes, unchanged serialized
bytes, failed task propagation and caller-owned input snapshots.

All four concrete compiler paths report zero findings. Package report messages
match the prior baseline exactly: 26 on 1.12, 25 on 1.13, with the unchanged
ceiling of 30. Logs: `/tmp/hybrid-publication-v2-jet112.log` and
`/tmp/hybrid-publication-v3-jet113.log`; v2/v3 production behavior is identical.
Final docstring-only additions were checked separately. JuliaFormatter 1.0.62
is idempotent on all final source/test files and `git diff --check` passes.

Each version passes seven allocation/encoding checks for 128 B, 4 MiB and 16 MiB
payloads. Additional publication allocation is 3,808–3,824 B on 1.12 and 3,712 B
on 1.13; the samples show no payload-sized duplicate buffer. Logs:
`/tmp/hybrid-publication-v3-resources112.log` and `113.log`.
The final reviewed source/test hashes are retained in
`/tmp/hybrid-publication-v3-handoff.json`.

This evidence is from local macOS runs. Linux CI and the next immutable full
package run provide separate integration evidence. Directory write access and
exclusion of concurrent directory mutation remain requirements. Atomic rename
preserves the previous file on pre-publication failure; it does not provide
crash durability or update other hardlinks to the old inode. Configuration
serialization of internal aliases/native handles remains separate IO-010 work.
