# Rank-owned isothermal kinetic-ion hybrid

`distributed_hybrid(g,species,physics)` source-implements periodic 3D CPU
Float32/64 kinetic ions and massless isothermal electrons on actual PencilFFTs
physical pencils. It shares the canonical ion prime/Boris/joint midpoint drift,
Ohm scalar rows, Faraday sign and RK4 affine owners. It does not wrap a serial
`HybridStepper`, construct a substitute `FourierGrid`, or replicate global fields
or populations. This bounded model is one CAP-003 family, not CAP-003 completion.

The current implementation-first pass adds and runs no tests. Package loading,
numerical/convergence, actual-rank rollback/liveness, compiler/JET, formatting/CI,
restart continuation and complete-operation memory/performance acceptance remain
deferred. The [complete matching-layout checkpoint source](DISTRIBUTED_HYBRID_CHECKPOINT.md)
is implemented; its runtime acceptance is also deferred.

## Construction and state

Load PencilArrays/PencilFFTs and HybridPlasmaPIC, construct a
`distributed_fourier_grid`, then collectively call
`distributed_hybrid(g,species,HybridModel(IsothermalElectrons(Te);nfloor))`.
Te has the grid's Float32/64 type and is finite/nonnegative. The original model
η, ηH and de2 must be stored zero. Original nfloor and its converted grid value
must be finite and positive. Nonzero coefficients are not allowed to disappear
through conversion. The original ion q/m factors feed the canonical mover.

Every rank supplies the same nonempty ordered Vector of species slots and q/m,
including locally empty populations. Leaves are independent ordinary Vectors:
positive finite q/m, finite x/v and nonnegative finite weights. All populations
are mobile, nonrelativistic and CIC. At every moment stage, represented deposited
n must be at least nfloor. Insufficient coverage rejects; no density clipping or
hidden background is used. There is no |v|<c gate in this nonrelativistic model.

`model.fields` exposes B/E/n/ui/pe/J/gradp/ninv, eighteen owned real pencil arrays.
Seed finite B with `PencilArrays.global_view(model.fields.B[c])`, then call
`init!(model,species)` with physical full-time velocities. Init privately wraps,
migrates, forms moments and Ohm fields, then publishes fields/particles/reset
phase together. Old E and other derived fields may be stale. Reinit never infers
or undoes a previous velocity half-phase. Initial B should be solenoidal; no
implicit projection or filter is applied, and diagnostics report represented
magnetic divergence.

Particle ownership is the actual physical-pencil owner of the wrapped CIC base
node. The existing partition lookup, eight-term stencil, routed exact moments,
gather replies and arbitrary-distance migration are reused. Repeated stencil
nodes on tiny axes retain all contributions. Local counts can change through
migration. Owner objects remain stable; component leaves may be replaced, so
references to old leaves do not track migration. IDs/tags persist without ID
renumbering or a uniqueness certificate.

Ready step, diagnostics and save enforce cached owner order/q/m/local lengths.
Explicit init can bind new counts and factors with reset phase. Same-count
replacement and source/velocity/weight/field edits cannot all be detected and
are unsupported between ready operations. The count cache is not a content or
identity certificate.

Calls are collective and nonconcurrent, nearest rounding and gradual underflow
only. Grid metadata is read-only. The communicator/FFT plan remains caller-owned.
Base alias reporting is the supported overlap contract; fields, particle leaves,
bindings and private work must be independent and disjoint from FFT workspace.
There is no raw-pointer overlap framework or device/replicated fallback.

## Normalization and numerical map

The model always uses positive-ion **collection** semantics, even one species:

```text
n  = sum(q*w*S(x))/stored_cell_volume
Ji = sum(q*w*v*S(x))/stored_cell_volume
ui = Ji/n
J  = curl(B)
pe = Te*n
E  = -ui×B + (J×B)*ninv - grad(pe)*ninv, ninv=round_T(1/n)
dB = -curl(E)
```

The deposit combines exact original q/weight/stencil terms across ranks before
volume division and one T rounding. This differs from serial intermediate
summation, so bitwise serial/distributed trajectories are not promised. The
pressure product, derivative and stored reciprocal remain in canonical order;
there is no `Te*grad(log n)` or product-rule rewrite. Normalization/clipping
owners are reused only after strict density validation proves clipping inactive.

At an unstepped ready boundary x/B/E and v are full-time, step/time/history are
zero, and primed is false. At an advanced boundary x/B/E and exposed predictor
moments are nominal integer-time; v has offset `-last_dt/2`. The last accepted
interval controls the next kick. Clock need not equal step times current dt.
B has no Maxwell half-phase or retained CAMCL pair.

`step!(model,species,dt;NB=1)` uses positive integer NB magnetic subcycles:

1. Validate ready state, original/configuration agreement, converted nonnegative
   dt, resolved nonzero duration, ownership and finite advancing clock. Stored
   zero dt validates and leaves private work/phase unchanged, without consuming
   initial priming.
2. Gather old carried E/integer B to detached ions. On first use back velocity
   up by the canonical h/2 Euler prime. Boris spans h initially and
   `last_dt/2+h/2` thereafter. The shared joint kick/coordinate admission records
   midpoint and endpoint positions, then wraps both periodically.
3. Route midpoint q-weighted n/Ji, normalize, and prepare pe/gradp/ninv once.
   These moments stay frozen throughout every RK stage and subcycle.
4. Perform the same four-stage RK4 map for represented `hb=h/NB`, using existing
   original-h affine completion owners. Each local affine stage is validated
   collectively before the following RHS communication.
5. Form preliminary endpoint E from actual half-time velocities and Bnew.
   Gather that same preliminary E/B for every species, predict v at integer
   time into private buffers, redeposit predicted moments and recompute carried
   E. Actual stored v remains at its half-time level.
6. Migrate endpoints and prepare all replacements/count bindings/final fields
   before one final all-rank barrier. Then publish all eighteen fields, complete
   particle leaves, bindings and phase with ordinary stores only.

Finite raw moments, normalized velocity, pressure, reciprocal, gradient, J/E,
Faraday rates and completed RK stages must be representable. A completed local
Ohm E can be finite while its original product source needs special Faraday
handling. The serial original-source recovery is tied to its own Fourier/dyadic
owner. This distributed family **rejects a canonical source-risk event before
curl(E) in every differentiated RK RHS**. It does not silently discard that
route or copy serial FFT physics. Init/carried/predictor E may use canonical
completed scalar rows when their physical outputs are finite. Full actual-pencil
original-source Faraday overflow/underflow recovery remains open. Likewise a
nonfinite materialized curl/rate rejects at its boundary, without asserting that
a later duration-scaled complete expression is impossible.

All output candidates remain private until complete success. Recoverable local
errors agree before the next collective and preserve public fields, particles,
counts and clocks. Private work may change. Communicator/transport failures,
fatal-resource failures and concurrent mutation are outside that transaction.

`recommended_dt(model;NB,safety,d_i,Ω_ci)` shares the existing magnetic RK4
advisory and exact/bounded fallback on the actual spacings. It does not establish
nonlinear hybrid, particle or ion-acoustic stability, and is not an automatic
step limiter. NB may change: RK stages have no retained subcycle phase.

## Diagnostics and persistence

`hybrid_diagnostics(model,species)` reports original-factor count, charge,
momentum and kinetic energy; magnetic energy; `Te*integral(n*log(n))`; their
complete global total; density extrema and represented divB. No E² reservoir is
included. It reports physical ion charge and represented electron charge
separately, with deposition/rounding and stored-cell-volume/declared-volume
ledgers. There is no Poisson solve, physical-neutrality rejection, background or
DC neutralizer: electron density is the represented deposited positive charge.
No exact discrete energy, Gauss or continuity guarantee is implied.

Exact kinetic/magnetic partials are combined before final rounding. Electron
energy keeps each rank's canonical log-coefficient dictionary local. Only exact
rational lower/upper endpoints cross ranks per common precision. Global any-odd
and log(2)-coefficient bookkeeping identifies the zero-log case; positive endpoint
odd coefficients preserve the canonical noncancellation premise. This is not an
endpoint-difference or mixed-sign-odd evaluator. Total energy adds the exact
mechanical offset before interval rounding. Complete diagnostic outputs outside
T may report Inf through canonical rounding, without changing physical state.

Local dictionary size can be O(owned cells), exact scalar bit length and precision
rounds vary, and cost is unmeasured. No full dictionary/field is globally gathered.
The diagnostic retains nominal particle `-last_dt/2` and integer-time field/moment
offsets; it does not time-center ion kinetic energy. Private diagnostic scratch
may change, public state does not.

A separate matching-layout checkpoint retains all eighteen public fields,
complete ordered particles and initialized/primed/time/step/last_dt verbatim.
Carried E and predictor ui cannot be rebuilt by depositing half-time velocities.
Counts are derived on restore; no RNG, RK stage, midpoint, gather, plan,
communicator or layout lookup is serialized. Stale destination physical/qm/count/
clock values are replaceable after representation/configuration checks. Restore
performs no moment/Ohm/FFT reconstruction. See the
[checkpoint guide](DISTRIBUTED_HYBRID_CHECKPOINT.md);
numerical/compiler/rollback/resource acceptance remains deferred.

The first slice excludes CAMCL, pressure evolution/transport/CGL, resistivity,
conduction, electron inertia, collisions/sources/RNG, relativistic ions, shape or
species schedules, boundaries, accelerators and rank-layout changes. Those and
the remaining distributed original-source range work keep CAP-003 open.
