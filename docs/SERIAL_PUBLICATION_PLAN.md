# Atomic publication of serialized outputs

Status: implementation in progress (IO-011), following the verified checkpoint
configuration repair `3495c64`.

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
These are prototype results; final implementation gates remain pending.
