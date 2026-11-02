# Physical-pencil Maxwell checkpoint

`save_mpi_checkpoint(dir, model)` and `load_mpi_checkpoint!(model, dir)` add a
source implementation of matching-rank/layout restart for
[`distributed_maxwell_vacuum`](DISTRIBUTED_MAXWELL_IMPLEMENTATION.md). They are
collective on the model grid's communicator. Every rank must call the same
operation with the same directory and immutable configuration.

Runtime/restart/liveness, compiler/JET, formatting/CI and resource acceptance
remain deferred under the current implementation-first policy. No tests or
runtime execution were added or performed for this source slice. It supplies
neither changed-layout restart nor completion of CAP-003.

## Saved state and compatibility

The dedicated `pencil_maxwell_manifest.ser` uses a checksummed envelope around a
distinct schema-1 manifest. It records the periodic spectral KDK/exact-affine,
represented-curl and matching-physical-pencil policies; precision, global n/L/dx,
speed and rank count; synchronized time/step; a unique `_pencil_maxwell_...`
generation; canonical ordered rank filenames; and each rank file's SHA-256 digest.
The manifest payload's checksum is verified before that payload is decoded.

Each rank file stores its zero-based MPI rank, logical owned ranges, actual
memory permutation and parent-array dimensions, plus six detached ordinary CPU
Array fields in parent memory order. It repeats the common configuration/phase.
Restore requires exact current ownership and permutation, even when local array
sizes alone would match. Plans, pencils, communicators, Fourier metadata and
scratch are not serialized. No RunMetadata or artificial RNG seed is required.
There is no schema migration, cross-model interchange or implicit repartition.
These are trusted Julia Serialization files; digests detect changed bytes and
do not authenticate an untrusted payload or promise cross-version portability.

## Save publication and cleanup

Save agrees mode/configuration/path and the current valid phase, then validates
storage and finite physical values collectively before snapshot allocation. The
six arrays detach before filesystem publication. Caller paths retain symlink and
`..` semantics: directory creation uses the supplied path, then ranks agree its
physical `realpath`; no prior lexical `abspath`/`normpath` changes its meaning.

Rank zero validates any superseded manifest and streams its rank-file digests,
then creates a unique generation. Generation-name encoding, receiver allocation
and decoding are guarded between collectives. Each rank atomically writes and
closes only its own new canonical file. After all writers finish successfully,
rank zero streams the completed files, builds the manifest envelope and publishes
it through the existing atomic serializer. No global field is assembled in memory.

Failure before manifest publication retains the old manifest. Cleanup removes
only the new generation's canonical rank files and its empty directory; unknown
entries are retained and reported. Generation/rank symlinks reject, and cleanup
never recurses. An original failure and cleanup failure are both retained.
After publication, failure to clean the superseded generation reports that the
new checkpoint is already usable. The directory requires a shared filesystem,
exclusive save/load access and no concurrent mutation. There is no power-loss
durability guarantee; process/transport loss and fatal resources retain the
existing MPI exclusions.

## Restore transaction

Restore validates destination storage/layout/aliases and immutable configuration,
without rejecting stale destination field/time values. All ranks validate and
agree the manifest, then each rank verifies its file digest as a stream before
decoding its detached state from that same open file. Required fields/types,
policies, phase, logical ranges, permutation, array shapes, finite values and
input/destination independence validate collectively before any live store.

Only after the last success agreement are six ordinary array copies and time/step
stores performed. Public arrays and clock Ref identities persist; fields are not
recomputed. Scratch is not restored and future stepping overwrites it. Recoverable
read/checksum/validation failure leaves physical fields and clocks unchanged on
every rank, within the stated filesystem/MPI premises. Actual continuation and
failure behavior still require the deferred runtime validation pass.
