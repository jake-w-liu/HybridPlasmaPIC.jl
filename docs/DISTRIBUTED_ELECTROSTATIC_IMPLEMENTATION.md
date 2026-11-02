# Rank-owned periodic electrostatic PIC

The source implements a 3D CPU Float32/64 electrostatic collection on actual
physical PencilFFTs pencils: owned density/E, rank-local particles, complete CIC
stencils, direct migration and a common leapfrog phase. This is the first coupled
CAP-003 slice. Runtime, numerical, compiler, format and resource acceptance remain
deferred during the implementation-first pass; no tests or executions were added.
Matching-layout checkpoint is the next source slice, not supplied by this one.

Load PencilArrays/PencilFFTs with HybridPlasmaPIC, construct a
`distributed_fourier_grid`, then call
`distributed_electrostatic_pic(g, species; grid_dc=:project, ...)` on every rank.
The species collection itself is an ordinary Vector. Every rank supplies the same
number/order of species, with the same q/m; each local population may be empty. Particle leaves must be independent ordinary Vectors of
the matching precision with finite positions/velocities, nonnegative finite
weights, finite signed charge and positive finite mass. The mover retains original
q/m factors. Weights denote particle number; epsilon0 is normalized to one.

## State, phase and ownership

Call `init_espic!(model,species)` with physical full-time velocities. It wraps and
redistributes private copies, forms density/E, then publishes together, resetting
time/step and priming history. The factory alone does not initialize a field.
`step_espic!(model,species,dt)` gathers old E, calls the shared electrostatic Boris
prime/kick/drift owner, migrates to the new owner and solves the new source.
After stepping, velocities are stored at the half level; diagnostics reports its
negative last_dt/2 offset. Variable steps retain the preceding interval convention
of the serial owner. This is not an exact discrete energy conservation law or an
unconditional timestep stability guarantee.

Owner tables come from actual logical pencil ranges and topology rank mappings,
with identical/congruent communicator ordering. A particle belongs to the rank
owning its wrapped CIC base node. Its eight terms retain repeated nodes on tiny
grids; transport deduplicates only requested field nodes. Remote positive-one-cell
face/corner stencil values and exact charge contributions are routed only to their
owners. Owned fields are accessed with global logical views, not raw parent
indices. No global field/population is gathered or replicated.

Migration can cross any number of domains in one step. It retains every position,
velocity, weight, ID and tag; no particle is created/deleted. Population owner
objects remain stable, but their component arrays can be replaced. Retained
references to old component arrays do not follow migration. Caller-provided IDs
are not silently renumbered or certified globally unique.

After manually changing positions, weights or factors, call `poisson_E!(model)`
to wrap/redistribute and refresh density/E without changing the existing clock or
velocity half level. Added particles must already have velocities at that level.
Use initialization to replace/reorder owners and reset phase. The constructor
species count is fixed. Ordinary stepping requires canonical periodic positions,
correct local ownership and the bound factors used for its carried field.

## Source and neutrality policy

The three existing canonical owners are shared with serial code: periodic cell
coordinates, original-factor exact moment contributions and one vector gather
stencil row. Charge contributions cross ranks exactly before normalization by the
product of original stored dx values. Each complete cell density, including the
uniform background, rounds once. This deliberately defines a new density map;
serial normal-path density bits are not promised. Exact arithmetic and byte
serialization have unmeasured cost. Native gather and the canonical pusher retain
their existing operation order given equal inputs. Finite materialized gathers
are required; no complete original-field gather range recovery is claimed.

`neutrality=:require` applies the serial physical-charge check and declared
`charge_rtol` (default 64eps(T)) to original particle/background charges.
`:neutralize` instead declares a uniform physical-policy compensator. Neither
changes weights. The required **explicit** `grid_dc=:project` independently chooses
the canonical Poisson `:neutralize` policy for the represented deposited density.
Standalone Poisson's strict exact-zero `:require` admission is unchanged.

The exact diagnostic ledgers separately report particle charge, physical
background charge and neutrality residual/compensator; deposited-grid charge and
its DC compensator; declared box volume and stored grid volume; background volume
discrepancy and particle-deposition discrepancy. These explain why physical
neutrality does not imply an exactly zero rounded grid sum. Physical-policy and
grid compensators describe the same projection at different representations;
they are not two added backgrounds. No hidden exact charge equality is asserted.

The [canonical Poisson field](DISTRIBUTED_POISSON_FIELD.md) sets derivative-null
DC/Nyquist combinations to zero. Uniform DC projection does not solve every
unresolved Nyquist charge mode. Existing FFT scale, coefficient, finite-output,
conditioning and range boundaries remain; no full Gauss guarantee for such modes.

## Collective publication and diagnostics

Nearest rounding and gradual underflow are required. Converted dt is nonnegative
and either true stored zero or at least 2floatmin(T). A zero step validates then
leaves physical state and scratch unchanged. Nonzero steps require finite clock
advancement. All ranks must invoke model APIs in the same order. The grid and
communicator remain caller/dependency owned; concurrent use is unsupported.

Private particle copies, gathered fields, migrated complete replacements and all
four density/E fields are prepared before publication. Every fallible local
transport preparation/assembly and arithmetic stage agrees before the next
communication. Final publication only rebinds complete same-type particle leaves,
copies fixed ordinary field arrays and assigns typed phase/clocks. Recoverable
prepublication failure preserves physical state on every rank; private work may
change. MPI loss and fatal resource failures retain existing exclusions. Aliases
must be exposed by the supplied ordinary storage and Base.mightalias contracts.

`electrostatic_diagnostics(model,species)` exchanges exact scalar charge, count,
kinetic and field-energy contributions. Energies round once globally and the
field is counted once; an unrepresentable total returns Inf. No field is gathered.
The source provides no global momentum diagnostic, changed-layout restart,
collisions, subcycling, hybrid/EM composition, GPU path or measured cluster scaling.
