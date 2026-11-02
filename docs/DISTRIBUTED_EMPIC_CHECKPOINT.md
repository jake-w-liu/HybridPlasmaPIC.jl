# Physical-pencil electromagnetic PIC checkpoint

`save_mpi_checkpoint(dir,model)` and `load_mpi_checkpoint!(model,dir)` provide
matching-rank/layout persistence for
[`distributed_empic`](DISTRIBUTED_EMPIC_IMPLEMENTATION.md). Both operations are
collective on the model grid's communicator. Every rank supplies the same
directory and cached configuration. Save requires an initialized boundary;
loading restores that boundary without another initialization.

This is a source implementation. Module loading, continuation, collective
failure/liveness, allocation and resource acceptance remain unexecuted. No tests
or runtime execution were added or performed for this slice. Changed-layout
restart and completion of EM CAP-003 remain outside it.

## Physical state and exact matching

The dedicated `pencil_empic_manifest.ser` wraps a checksummed schema-1
`:distributed_empic` manifest. Its explicit policies identify single-rate
nonrelativistic, all-mobile CIC/Boris evolution, carried magnetic/velocity phase,
original-factor exact charge routing with one density rounding, represented
odd-Nyquist-zero derivatives, and corrected midpoint current that retains raw
true DC and zeros other fully derivative-null current modes. These identify the
implemented family, rather than selectable checkpoint options.

Precision, global n/L/dx, speed c, background charge, physical neutrality
policy/tolerance, explicit `grid_dc=:project`, ordered cached species q/m,
communicator rank count and actual rank-indexed owned ranges must match exactly.
Each rank also matches its logical ranges, actual memory permutation and parent
dimensions. Matching dimensions alone do not establish a matching layout.
Species count/order stay fixed; local particle counts may change, including zero.

Each rank retains eleven independent ordinary CPU arrays in parent memory order:
three components each of E, B and J, plus rho_n and rho_np1. It also retains every
ordered local x/v/weight/ID/tag/q/m record. Empty populations remain explicit.
The manifest and rank records carry initialized/primed flags, the sole last_dt
history, time and step. At count zero, time/history are stored zero, velocity and
magnetic field are unprimed, J is stored zero and the density endpoints are
identical. Advanced boundaries have positive finite time, primed state and
positive resolved history strictly below the model's derived CFL. Variable
timesteps are supported; time is not reconstructed from step times last_dt.

E and positions are at integer time. Initially B and velocities are full-time;
after a positive step they have the nominal offset -last_dt/2. J is the completed
last corrected midpoint current, and rho_n/rho_np1 are that interval's preceding
and new endpoints. Loading copies all eleven arrays verbatim. It does not infer
old current/density from current particles, redeposit, project, recenter, migrate,
renumber IDs or prime the mover. No exact floating Gauss/continuity identity or
zero rounded-grid DC is imposed as a checkpoint gate. Original physical
neutrality is checked using incoming populations and the bound policy.

The local count cache is derived from fully prepared incoming populations; it is
not serialized or accepted as independent source state. Charge/DC/volume,
energy, momentum and residual diagnostics derive from saved physical inputs and
configuration. Plans, communicators, layout partition/lookup caches, field and
particle candidates, Bc/Jraw/curls, source spectrum, midpoint/ghost and FFT scratch
are omitted. The matching destination supplies valid storage for later staging.

Unknown keys, versions, tags and policy values reject, as do serial EM, ES and
Maxwell payloads. There is no schema migration, repartition or cross-version
portability promise. These are trusted Julia Serialization files: SHA-256
detects changed bytes, rather than malicious payloads.

## Save publication and directory reuse

Save collectively validates current storage, values, cached factors/counts,
phase, particle ownership and original physical neutrality before detaching the
local snapshot. Unsupported external population mutation requires the model's
explicit reinitialization before ready operations.

An existing directory must match the current cached configuration, including
q/m/c and policies. After reinitialization changes a cached binding, use a new
directory. An incompatible manifest rejects before a new generation is created,
preserving its previous checkpoint. Local count changes alone do not change the
cached configuration.

Ranks create the supplied directory and then agree its physical `realpath`,
retaining symlink/`..` path semantics. Rank zero validates the old manifest and
streams every old rank digest before creating a unique `_pencil_empic_...`
generation. Generation-name encoding, receive allocation and decoding are
guarded between collectives. Each rank atomically writes and closes its own
canonical file. Only after all writers succeed does rank zero stream the new
file digests and atomically publish the checksummed manifest. Its payload digest
is verified before decoding. No global field or population is assembled.

Recoverable pre-publication failure preserves the previous manifest. Cleanup
removes only canonical known files and their empty generation directory;
generation/rank symlinks reject, unknown entries remain and are reported, and
cleanup never recurses. Original and cleanup exceptions are both retained. A
superseded cleanup failure reports that the new checkpoint was published.
Operations require a shared filesystem, exclusive save/load access and no
concurrent model mutation. No power-loss durability is promised. MPI/process
loss and fatal resource errors retain the existing collective-stage exclusions.

## Restore preparation and publication

Restore requires valid ordinary destination storage, aliases, actual layout and
cached immutable configuration. Old particle q/m/x/v/weights, physical fields,
phase/clock values, local lengths and count-cache entries may be stale. It does
not call current ready/phase/count agreement. Incoming q/m matches the cached
binding, rather than replaceable destination particle q/m.

All ranks agree the verified manifest. Each rank streams its file digest before
decoding from that same open file, then checks exact schema/configuration/rank,
layout, typed phase/history, initial current/endpoints, finite fields/particles,
source ownership and storage independence. Incoming buffers cannot alias one
another, public/private particles, fields, spectra, metadata, caches or Fourier
storage. Every same-type replacement leaf, typed scalar row and derived count
row is allocated before the final collective physical-neutrality check.

After all-rank preparation succeeds, publication only rebinds complete particle
components, assigns typed q/m, copies eleven fixed field arrays and prepared
scalar/count rows, and assigns typed phase/clock Refs. There is no public
resize/append, late validation, communication or physical reconstruction after
the first store. ParticleSet owners, field arrays, cache Vectors and phase Refs
retain their identities; old component-array references do not follow particle
replacements. Scratch is reused through normal later staging.

Within these premises, recoverable read/checksum/configuration/preparation
failure leaves public particles, fields, caches and phase unchanged on every
rank. Ordinary publication excludes custom setters, protected memory, concurrent
mutation and fatal resource failure. Its runtime continuation and failure
behavior remain unexecuted under the source-only instruction.
