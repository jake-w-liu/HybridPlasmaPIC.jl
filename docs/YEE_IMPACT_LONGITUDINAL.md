# Event-local longitudinal reflection source candidate

`YeeImpactLongitudinalPIC` is a distinct nonrelativistic CPU Float32/64 model.
The source implements zero contacts or one certified common contact cluster per
external step, using at most two synchronized local kick-drift-kick chunks.
Exact-revision GitHub runtime, convergence, rollback/restart, compiler, formatting
and complete-step resource acceptance remain pending. This is an additional
bounded capability, not completion of all bounded PIC or a material sheath model.

```julia
state = YeeImpactLongitudinalPIC(grid;
    shape = YeeTriangularCloud(0.125), max_root_refinements = 128)
init!(state, particles)
step!(state, particles, dt)
diagnostics = yee_impact_longitudinal_diagnostics(state, particles)
```

Use a closed 1D YeeGrid and explicitly select a fixed physical cloud half-width.
Particles may be one ParticleSet or an ordered one-based vector. The first
successful positive converted dt fixes the external step. The default root cap is
`precision(T)`; an explicit positive integer override changes only the finite
bisection budget. It is not a tolerance and does not make an inconclusive
physical certificate successful. Nearest hardware rounding, gradual underflow,
zero transverse velocities/E and all B are required. Original width/background/
left-plate Float32/64 types and bits are preserved; other real source keywords
resolve once to grid precision.

Initialization uses the existing complete original-source Gauss expression and
isolated free plates, resets histories and binds ordered particles, source values,
counts, provenance and component identities. Public completed centers must be
strict interior. Reinitialize explicitly after source/configuration changes.
Reflection acts on the cloud center while retaining particle-bound endpoint
tails; no charge is absorbed or mirrored.

Both walls and all quadratic branches are checked, including interior extrema.
The method completes the event-duration final electric half-kick before reversing
its incoming velocity, then advances the exact rational remainder. Internal
durations/halves, gathers and q/m retain their original rational values through
canonical Boris. Actual constrained endpoints deposit the canonical current.
For `x_free=x+s*u`, the signed snap is `r_snap=wall-x_free`; the predictor residual
is separately `r_pred=r_snap+s*(u-u_star)`. Snap/history/source-sensitivity gates
use actual r_snap, while event-time conditioning uses r_pred. A four-ulp bound
on `|x|+|s*u|` is an explicit endpoint policy, not an unreported nudge.

The [conditional carried-error contract](BOUNDED_EVENT_COVERAGE_CONTRACT.md)
retains physical hit-time intervals, position/velocity/force uncertainty and
the signed field-invariant defect through reflection and the enlarged remainder.
A fresh zero-error tube from the numerical contact state is insufficient.
Only a just-certified selected center may privately depart its own wall; its
inward velocity must be strict, and only the departure root at zero is ignored.
Later same-wall returns and opposite/other contacts remain checked. Grazing,
terminal/public contact, unresolved ordering, a second cluster or exhausted
refinement reject the whole external candidate. These conservative gates do not
establish a general stability or entire-run initial-trajectory accuracy theorem.

`last_work` is once-rounded complete signed work from the actual accepted local
fields/currents. `last_impulse` is a once-rounded net diagnostic of the stored
local currents, not either actual stage current. `last_energy_defect` rounds
complete rational DeltaK minus exact local work sum; it differs from the old
event-free longitudinal rounded-kinetic diagnostic. `last_physical_energy_change`
rounds complete physical Delta(K+U). Neither is forced to zero. Read-only energy/
Gauss can overflow; stored energy/work histories must be finite before publish.

Per-wall last/cumulative counts include zero-weight tracers. Exact signed wall
impulse `2*m*w*v_incoming` is delivered to each wall; represented electric-kick
impulse is separate. Fixed elastic wall charge and work are zero by law, without
an artificial heat reservoir. This ledger does not include conductor/background
stress or claim conserved plasma+wall momentum. Exact scalar histories report
snap, predictor/physical time bounds, source-invariant and terminal state bounds,
and a snap-energy sensitivity bound; they are not conservation corrections.

All physical/clock/history work stages before publication. Ordinary writable
CPU strided storage and public identities are preserved; scratch can change on
error. Special-store ancestry and forbidden aliases reject. OS protection
changes, arbitrary custom setter failure and concurrent mutation/stepping remain
outside that ordinary-write guarantee.

Distinct checkpoint tag `:yee_impact_longitudinal_pic`, schema1, saves completed
phase, original configuration bits, verbatim fields, ordered particles and every
exposed ledger/history. Required versioned policies identify the numerical map,
endpoint/coverage/cluster scope and zero wall charge/work law; max-root budget
must match. Scalar/collection layouts and other model schemas cannot mix.
Restore detaches all incoming values, rebuilds both private stages and reserves
Vector capacity before publication; matched strided destinations retain their
shape. No field projection or reconstruction occurs. Archives use the existing
schema4 metadata/checksum container and atomic publication owners.

The existing longitudinal manual suite includes compact independent physical-
volume/current/work, bounded travel-time refinement, mixed selected/nonselected
particles, carried remainder clearance, common-time rejection, branch,
late-failure/retry, view/empty/restart and concrete compiler controls. They remain
unrun for this candidate. Absorption/emission, circuits, multiple clusters, transverse walls,
relativity, GPU/MPI and broader physical/cost acceptance remain separate work.
