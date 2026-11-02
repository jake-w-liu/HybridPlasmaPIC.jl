# Physical-pencil electrostatic checkpoint

`save_mpi_checkpoint(dir,model)` and `load_mpi_checkpoint!(model,dir)` add a
source implementation of matching-rank/layout restart for
[`distributed_electrostatic_pic`](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md).
They are collective on the model grid's communicator. Every rank calls the same
operation with the same directory and cached model configuration. Save requires
an initialized state; initialization resets the clock and velocity phase, so do
not initialize again after loading.

Runtime/restart/liveness, compiler, format/CI and resource acceptance remain
deferred. No tests or runtime execution were added or performed for this slice.
It does not provide changed-layout restart or completion of CAP-003.

## Saved state and exact matching

The dedicated `pencil_electrostatic_manifest.ser` contains a checksummed envelope
around a distinct schema-1 `:distributed_electrostatic_pic` manifest. It records
the shared CIC/leapfrog, original-factor exact charge routing/one density
rounding, represented odd-Nyquist-zero Poisson and matching physical-pencil/CIC
anchor policies. Precision, global n/L/dx, background charge, physical neutrality
policy/tolerance and explicit `grid_dc=:project` must match exactly. So must the
cached ordered species q/m, communicator rank count and actual rank-indexed
owned ranges. Local particle counts may change on restore.

The manifest retains initialized/primed flags, last_dt and time/step, a unique
`_pencil_electrostatic_...` generation, canonical ordered rank filenames and
their streamed SHA-256 digests. Its serialized payload checksum is verified
before that payload is decoded. Count zero has zero time, unprimed velocity and
zero history; an advanced count has positive time, primed velocity and positive
resolved history. Variable timesteps are supported. Time is not reconstructed
from count times timestep.

Each rank file repeats configuration/phase and its zero-based rank, logical
ranges, actual memory permutation and parent-array dimensions. It stores three
ordinary CPU Array E components and rho in parent memory order, plus ordered
local x/v/weight/ID/tag/q/m records. Empty populations are retained. Fields and
particles restore verbatim; there is no Poisson solve, recentering, migration,
charge repair or ID renumbering during loading. The physical/grid charge,
compensator, volume and deposition ledgers are derived from these saved values
and configuration; no separate mutable ledger state exists.

Plans, communicators, Fourier metadata, derived partition/lookup caches, ghost
maps, particle candidates and FFT scratch are not serialized. A matching model
supplies its own valid layout and scratch. There is no cross-model interchange,
schema migration, repartition or cross-version portability promise. These are
trusted Julia Serialization files: digests detect changed bytes, not malicious
payloads.

## Save publication and directory reuse

Save collectively validates current storage, values, bound factors, phase,
particle ownership and original physical neutrality before detaching the local
snapshot. It requires the usual consistent carried field/source state; manual
source changes require the model's explicit field refresh before saving.

An existing manifest must match the current cached configuration, including
q/m. After an allowed scalar-rebinding refresh changes q/m, use a new directory.
An incompatible existing manifest is rejected before a new generation is
created, retaining its valid checkpoint. Changing only rank-local particle
counts does not change the cached configuration.

Directory creation uses the supplied path, then ranks agree its physical
`realpath`; symlink/`..` semantics are retained. Rank zero validates the previous
manifest and streams every old rank-file digest before creating a unique
generation. Generation-name encoding, receive allocation and decoding are
guarded between collectives. Every rank atomically writes and closes only its
own new canonical file. After all writers succeed, rank zero streams the
completed files and atomically publishes the checksummed manifest envelope.
No global field or population is assembled in memory.

Recoverable failure before manifest publication preserves the old manifest.
Cleanup removes only known canonical files and their empty generation directory;
unknown entries remain and are reported. Generation/rank symlinks reject, and
cleanup never recurses. Original and cleanup exceptions are both retained. A
superseded cleanup failure after publication reports the new checkpoint usable.
The directory requires a shared filesystem, exclusive save/load access and no
concurrent mutation. No power-loss durability is promised; MPI/process loss and
fatal resource errors retain the existing collective-stage exclusions.

## Restore transaction and retained identities

Restore validates destination ordinary storage, aliases, actual layout and
cached immutable policy. Old particle q/m/x/v/weights, field values, clocks and
phase flags may be stale or invalid values; current phase agreement is not used.
Saved q/m must match the cached model binding, rather than stale particle q/m.
Species count and order stay fixed.

All ranks agree the checksummed manifest. Each rank streams its file digest
before decoding the detached state from that same open file. Exact layout and
memory permutation are required even when array sizes alone would match. Typed
phase, finite fields/particles, nonnegative weights, positive masses, source
ownership and incoming/live/private storage independence validate before every
replacement leaf is allocated. Original physical neutrality is then checked
collectively using incoming populations.

After the final all-rank success, publication only rebinds complete same-type
particle components, assigns typed q/m, copies fixed ordinary field arrays and
assigns typed phase/clock Refs. It performs no public resizing, late validation,
communication or field reconstruction. ParticleSet owners, model field arrays
and phase Ref identities persist; component-array references retained before
loading do not follow replacements. Scratch resets through later normal staging.

Recoverable preparation/read/checksum/configuration/validation failure leaves
public state unchanged on every rank within these premises. Ordinary publication
excludes protected memory, custom setters, concurrent mutation and fatal resource
failure. Actual continuation, allocation, resource and collective failure behavior
remain unexecuted under the current source-only instruction.
