# Synchronous implicit electromagnetic collections

`ImplicitEMPICStepper` is a distinct coupled source implementation. It solves
particle velocities and both Maxwell fields together. Numerical, compiler,
restart and resource acceptance remain pending remote verification; source
implementation is not an exact energy or stability theorem.

## Supported boundary and API

The first family is periodic CPU Float32/64 D1–3 with all three velocity/field
components, ordinary Vector-backed ordered populations and CIC. Every population
is mobile. Charges may have either sign, masses are positive and weights are
nonnegative. Exact physical total particle charge must be zero; no fixed or
neutralizing background is supplied. Typed empty populations/collections give
the vacuum limit. Every grid length is odd, avoiding non-DC derivative-null
charge modes. Existing explicit EMPIC even-grid behavior is unchanged.

```julia
state = ImplicitEMPICStepper(g, CIC(), species; c = light_speed)
# Fill state.fields.E/B and synchronous particle values.
init!(state, species)
step!(state, species, dt)
diag = implicit_em_diagnostics(state, species)
save_checkpoint("state.chk", state, species)
load_checkpoint!(state, species, "state.chk")
archive_run("state.run", state, species; rng_seed = provenance_seed)
load_archive!(state, species, "state.run")
```

This is nonrelativistic motion. There is no speed clipping, collision/source/wall
map, subcycling, species schedule, orbit averaging or device/distributed path.
The normalization is epsilon0=1 with field energy `(E²+c²B²)/2`; `c` must resolve
positive and finite at grid precision. It is not an arbitrary SI-units adapter.

The model binds population identity/order and original q/m at construction.
Initialization atomically binds counts, weights, IDs and tags; those remain
fixed until reinit or restore. Replacing/reordering objects requires a new
model. External physical edits are the caller's responsibility: they must form
a consistent synchronous boundary, including fields/source handling. This
family supplies no charge-injection current or external-work ledger. Ready
checks catch source-binding and charge/Gauss inconsistencies; they do not replay
or certify the history of arbitrary external changes.

Positive `dt` uses the canonical resolved timestep floor and must advance the
stored clock. Original nonzero input converting to zero rejects; original zero
is a validated no-op. All physical particles and E/B are integer-time state.
A different later positive dt requires no leapfrog phase conversion. Variable
steps do not imply a proven accuracy guarantee. Last J/Jraw are interval
midpoint data; their reported offset is nominal, not an exact timestamp obtained
by subtracting dt from the rounded clock.

## Actual coupled map

Every trial starts from immutable old snapshots and one velocity/E/B unknown.
It forms the represented average velocity, lifted half/full coordinates and
periodic endpoints. Mean E/B are gathered at the lifted midpoint. The guarded
canonical Boris map retains original q/m and full stored dt; gathered samples
are finite rounded T operands. Original geometry-through-Lorentz completion is
not promised by this family.

Endpoint charge and raw midpoint current are redeposited inside every residual.
The canonical longitudinal spectral correction enforces represented continuity,
retaining true DC. The same trial mean fields enter both the particle residual
and coupled midpoint Ampere/Faraday residuals. No field/current is held to an
explicit predictor. There is no projection after the nonlinear solve.

The existing bounded Newton/forward-Jv/restarted-GMRES/Armijo algorithm is shared
with the synchronous hybrid family. Defaults are nonlinear_rtol=sqrt(eps(T)),
nonlinear_atol=64eps(T), max_newton=20, max_linear=200, memory=20,
max_backtrack=12 and max_jv_shrink=8. These are bounded policy, not convergence
claims. Three block scales stay fixed for the whole solve. Actual action and
fresh final residual gates are required; owned trial-domain failures alone allow
bounded probe/line shrink. Unexpected canonical errors rethrow. Failure does
not retry a different timestep.

Constraint defaults are constraint_rtol=256eps(T), constraint_atol=0. Positive
tolerances must survive conversion. Actual continuity/Gauss/divB defects are
normalized by fixed positive original-operand reference scales and gated
separately. Initialization corrects only longitudinal E and rejects inadmissible
seeded B. Rounded charge totals, FFT/store errors and source/current constraints
are measured; there is no bitwise Gauss-storage claim. Required finite curls,
FFT/current/gather intermediates are an explicit range boundary, even when a
more general completed expression might fit.

## Work and diagnostics

The detached `ImplicitEMSolveReport` saves actual residuals/scales/tolerance,
iteration/evaluation/backtrack counters, constraints and exact signed work.
Canonical exact kinetic, field and momentum owners supply endpoint differences.
It reports four pieces whose signed sum is the physical energy change:

- Particle map defect: deltaK minus actual particle electric work.
- Rounded midpoint/gather/deposition pairing defect: particle work minus raw
  grid-current work.
- Continuity correction work: raw minus corrected grid-current work.
- Field equation/mean/FFT defect: deltaU plus corrected grid-current work.

The continuity correction is a spatial change to current. Its work can persist
as dt tends to zero at fixed mesh. A tighter nonlinear tolerance need not remove
it. No term is forced to zero or called conservation. This differs from the
compatible interpolation used for the work theorem in
[Angus et al. (2023), equations9–11/18–19](https://debog.github.io/Files/2023_Angus_EtAl_JCP.pdf).

Diagnostics round complete exact reservoirs for display and retain true overflow
as Inf. Work defects are not computed by subtracting separately rounded Inf
energies. No exact momentum theorem or matched-error performance claim is made.

## Publication and persistence

One instance/population/grid permits one nonreentrant caller; shared FFT work
requires exclusive borrowing. Every input/storage/alias check and field,
particle, report and capacity calculation precedes one public commit. Expected
failure preserves all public physical/source-binding/clock/report state; private
work may change. Model, ParticleSet and field-array identities remain; particle
leaves may be replaced. Custom setters, OS protection changes and concurrent
mutation are outside the ordinary-storage publication contract.

Distinct schema1 `:implicit_empic_midpoint` saves all E/B/J/Jraw/rho_n/rho_np1,
ordered particles/provenance, fixed configuration/policies, actual synchronous
clock/history and detached report. Count0 requires zero currents, equal charge
snapshots, zero history and no report. Advanced history/report must agree and
last dt must not exceed the completed stored time; variable dt does not require
time to equal step count times last dt.
No nonlinear/FFT plans, iterates or warm guesses are saved. There is no RNG;
archive seed is provenance. Raw checkpoint publication is atomic without a
separate checksum; the existing run container4 provides the archive checksum.

Restore requires the same fixed grid/precision/qm/shape/solver/constraint policy
and population count/order. Destination values/counts/qm/clock/report/source
bindings may be stale; fixed storage and owner identities remain required.
Original incoming aliases are checked before whole-graph detachment, then
validation/preparation runs again. Counts/source bindings and all solver capacity
are prepared before one commit. Incoming fields/current/rho are copied verbatim;
private charge/constraint validation does not reconstruct public fields. No
priming, recentering, projection or nonlinear re-solve occurs. Other explicit or
hybrid payloads reject with no implicit migration.

All actual convergence/order/dispersion, failure/reuse/restart, compiler/format
and matched physical-error cost acceptance remains open. Larger stable-looking
steps or the vacuum linear CN response do not establish nonlinear stability or
resolved kinetic physics. Exact-head remote CI must close for every push; no
local Julia, tests, compiler, formatter or numerical/device workload was run.
