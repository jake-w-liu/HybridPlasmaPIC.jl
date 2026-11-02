# OpenHallHybrid checkpoint and archive

This is the distinct schema 1 persistence contract for the CPU D2/3
`OpenHallHybrid` numerical reservoir family. Source implementation does not
establish runtime, stochastic continuation, compiler or resource acceptance.
Exact integrated remote CI must close for every push; no local Julia or device
workload was run for this candidate.

Save only an initialized completed boundary:

```julia
save_checkpoint("hall.chk", model, species)
load_checkpoint!(model, species, "hall.chk")

archive_run("hall.run", model, species; rng_seed = provenance_seed)
load_archive!(model, species, "hall.run")
archive = load_archive("hall.run")  # inspect metadata/state without restoring
```

The supplied ordered population objects must be the model's bound owners.
`rng_seed` records caller provenance. Continuation uses the saved full Xoshiro or
MersenneTwister state, including its caches and supported metadata, and never
reseeds from that provenance value.

The raw checkpoint is atomically serialized and has no additional checksum.
The archive uses the existing run container schema 4 with a SHA-256 checksum of the
serialized state bytes. `load_archive!` verifies that container through `load_run`
and then invokes the same family restore as `load_checkpoint!`. Inspection via
`load_archive` retains the existing nonmutating API. File permissions, symlink
handling and same-directory atomic replacement follow the shared publication
owner; a rejected save leaves the previous published file intact.

The payload saves all of the following:

- Exact fixed geometry, physics, namespace/charge/mass bindings, CIC/SBP/derivative,
  SAT/source/gather/moment/phase/report policies, all controller/resource limits
  and converted tolerances.
- Complete ordered particles: positions, all three synchronous velocities,
  weights, UInt64 provenance IDs, UInt32 tags, original charges and masses.
- Every public `HybridFields` role, including density, ion velocity, B/E/J,
  pressure, pressure gradient, reciprocal density, zero lapJ, the empty scalar-CGL
  roles and the inactive floor counter.
- Actual stored initialized/time/step/last_dt and detached last scalar report.
- Exact source elapsed duration, fractional count carries, next and initial-next
  provenance indices, initial particle totals, last and cumulative per-face
  birth/exit aggregates and cumulative physical/work/energy ledgers.
- Full stock RNG state, its algorithm and the same-Julia/runtime-platform stamp.

Private FFT/projection plans, trial solutions, capacities, gather buffers and
controller guesses are omitted. The expected population-count cache is derived
from the incoming particles. The controller restarts its trial policy from its
configured first level rather than from a saved warm trial.

Restore requires matching fixed configuration, ordered species count,
namespace/original q/m binding, stock RNG algorithm and runtime/platform. Legacy
shock, periodic hybrid, other family/policy and different configuration payloads
are rejected; there is no implicit migration. The complete saved carried fields
are authoritative. Validation checks their finite/floor/zero/direct-constraint
requirements without projecting or otherwise changing incoming B.

Old destination physical values, public q/m, counts, clock/report and RNG
metadata values are replaceable. The destination must still have its fixed
ordinary storage, bound particle objects, geometry, cache shapes, independent
live/private roles and valid stock copy destinations. No view/custom-array
substitution or concurrent mutation is supported by this family.

Before copying incoming data, restore checks original physical/source/RNG
storage independence and rejects opaque RNG metadata linked into incoming or
destination-owned state. It then detaches the whole payload and validates it
again. All particle replacements, fields, source data, full RNG copy source and
work capacities are prepared before the first public store. The canonical
commit retains model, ParticleSet, field-array and RNG object identities; particle
leaves may be replaced as allowed by the model's publication policy.

Restore does not sample births, advance fields, redeposit moments, recompute
Ohm fields, normalize density, project B, kick, recenter or prime velocities.
The model is synchronous integer-time state, including physical loaded
velocities before the first completed step; it has no periodic leapfrog half-phase to reconstruct. No
`time=step*last_dt` equality is imposed. Exact source duration is saved separately
from the rounded public clock.

The saved report and algebraic source identities document the accepted map;
schema validation does not replay its field trials or prove past physical
accuracy. This numerical reservoir boundary is not transparent Hall outflow,
material/sheath physics, a stability theorem or a claim of complete open-boundary
energy conservation. Expected validation/preparation failure leaves public
physical/source/RNG/clock/report state unchanged; private work may change.
