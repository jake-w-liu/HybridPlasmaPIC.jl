# Serial Yee vacuum checkpoint and archive

`save_checkpoint(path, model::MaxwellVacuum)` and
`load_checkpoint!(model, path)` source-implement field-only restart for the
existing CPU Float32/64 Yee model in one to three dimensions. Periodic, PEC and
mixed axes retain their actual staggered component shapes. No tests or runtime
execution were added for this slice; restart, compiler/format/CI and resource
acceptance remain deferred during the implementation-first pass.

The distinct `:yee_maxwell_vacuum` schema 1 records the synchronized native/
normwise Yee KDK policy, precision/dimension, cell counts, lengths, stored dx,
periodicity, speed, six physical arrays and completed time/step. All grid and
speed values must match the destination. It does not reuse a particle-coupled
closed-Yee or distributed-pencil schema. No previous timestep is required by
this synchronized variable-step model; no clock is reconstructed as count*dt.
The next duration and any prescribed 1D current impulse remain caller inputs.
No waveform, particles, RNG, plan, staging or derived CFL is serialized.

Save validates storage independence, shapes, finite fields, PEC constraints,
speed and typed nonnegative clock before existing atomic serialization. The
ordinary checkpoint file has the existing trusted Julia Serialization format.
`archive_run(path, model; rng_seed, ...)` wraps the same state in the existing
schema-4 checksummed run container. The seed is provenance only. Default
normalization describes epsilon0=1, mu0=1/c²; the boundary description lists the
actual periodic/PEC axes. Explicit metadata overrides remain available.

Restore checks destination physical/private array shapes and independence,
allowing old nonfinite field values, invalid PEC entries or stale clock values
to be replaced. The full incoming policy, exact configuration, typed clock,
ordinary array shapes/precision, finite values, pairwise independence and PEC
constraints validate before publication. All incoming leaves detach together
and revalidate before six same-shape copies and typed clock stores. This also
handles an in-memory payload whose arrays permute live destination fields.
Physical arrays and clock Ref identities remain; fields are copied verbatim
without projection, reconstruction or boundary clamping. Scratch is not loaded
and the next ordinary step overwrites it.

The model and directory require exclusive access during save/load. The existing
atomic writer publishes only after successful serialization, flush and close;
there is no crash-durability or universal cross-version binary promise. Ordinary
copy guarantees exclude fatal resource exhaustion, OS protection and concurrent
external mutation. Runtime continuation and rejection behavior are unverified
until the deferred validation pass. Absorbing layers, coupled particle boundary
models, backends and broader CAP-002 integration remain separate work.
