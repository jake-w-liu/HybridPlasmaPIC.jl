# Synchronous open Hall hybrid

`OpenHallHybrid` is a distinct CPU Float32/Float64 D2/3 numerical reservoir
model: mobile positive ions, CIC, massless isothermal electrons, both x faces
open, and periodic transverse axes. It uses the actual mixed SBP/Fourier
operators and the restrictive perpendicular space `Bx=Jx=0` at both faces.
Bx is stored zero; differentiated Jx/divergence satisfy configured represented
residual gates, not a claim of exact native operator identities.
It supplies full nonlinear Hall and pressure induction; it does not turn the
existing reduced `PerpShock` into a different method.

This source candidate still needs exact integrated GitHub numerical, compiler,
formatting, continuation and resource acceptance. No local Julia workload was
run. The family requires both the advancing model and its distinct
[checkpoint/archive owner](OPEN_HALL_HYBRID_CHECKPOINT.md). Full CAP009 remains
open: this is not transparent Hall outflow, a continuum well-posedness or
stability theorem, material/sheath physics, or complete open-boundary energy
conservation.

## Construction and phase

Construct actual `MixedShockOperators(n,L)`, grid-precision
`HybridModel(IsothermalElectrons(Te); η=0, ηH=0, de2=0)` parameters with original
`η=ηH=de2=0`, and ordinary Vector-backed `ParticleSet` populations.
Pass `sources[s]=(left,right)` using same-precision `OpenHallReservoir` records.
Each source binds density, inward normal drift/spread, two tangential
drifts/spreads, positive macro weight and UInt32 tag. `targets` binds two
constant tangential B targets; `sat_speeds` binds nonnegative face coefficients.
The default seed is a complete rational convex interpolation between targets,
with Bx zero. A supplied seed is copied and must pass direct constraints at init;
it is not silently projected.

Bind distinct 16-bit species namespaces (default ordered 0,1,...) and nonzero
unique 48-bit particle indices. Use `assign_global_particle_ids!` explicitly to
prepare a load; the model never rewrites IDs implicitly. Charges/masses and
population object order are fixed. The expected count cache updates only after
prepared init, source/transport or restore publication. External ready-state
mutation is unsupported; same-count replacement cannot be detected.

Required controls are `rng` (Xoshiro or MersenneTwister), `max_levels`,
`max_births`, `max_particles` and `max_draws`. The supplied RNG becomes model-
owned; do not draw from it externally. Birth/particle limits include all predicted
births before exits, not only retained survivors. No limit automatically grows
or truncates the physical batch. Defaults are `field_rtol=sqrt(eps(T))`,
`field_atol=0`, projection/state/work relative tolerances `100eps(T)`, absolute
ones 0, `maxiter=200`, and `log_bits=precision(T)+16` (maximum 4096). Their adequacy
and cost are unmeasured.

Fill positions/velocities/weights/IDs, then call `init!(model,species)` once.
All x, v, B, E and moments are synchronous integer-time values. Initialization
consumes no draws and performs no velocity backing or priming. To reset, construct
a fresh model. Advance with `step!(model,species,dt)` or `step!(model,dt)`.
Original nonnegative finite dt must survive conversion to T. True signed zero
requires a ready finite/domain-valid state and leaves all physical/source/RNG/
report and private transform/solver work unchanged. It skips differentiated B
constraint checks; it does not certify newly mutated fields.

## Completed external interval

1. Hold old represented n/ui fixed and evolve B with the nonlinear mixed Hall
   RHS, explicit tangential SAT and constrained projected RK4 controller.
2. Complete E from the accepted B and those same frozen moments. Gather the
   complete original CIC coefficients at old positions and apply the existing
   original-q,m Boris velocity map for the whole converted duration.
3. Transport with fixed post-kick velocity, absorbing the first x exit and
   wrapping transverse displacement by one complete rational mod. Create
   reservoir births during this transport and drift each only for its remaining
   age. Newborns do not receive the preceding full kick.
4. Deposit strict endpoint moments and complete pressure/current/Ohm fields.
   Only after every candidate/source/RNG/ledger/report check succeeds do all
   public particles, fields, source state and clocks publish together.

The step is an explicitly synchronous first-order Lie map. No numerical order,
reflection or stochastic-stream acceptance follows from the source design.
Exact first-exit equality counts as an exit, including a through-flight birth.
Stable retained-old order precedes newborns, which use species/left/right/birth
order. Birth and exit both count for through-flight particles and consume their
reserved IDs. The last legal index is permitted; exhausted IDs are not reused.

For node g, charge and charge-weighted momentum use the exact represented CIC
stencil divided by actual `H[ix]*prod(dx_transverse)`, then the existing represented
ion normalization. Require finite n/momentum/ui and `n>=nfloor`; there is no
clipping. Gather uses the same unrounded per-axis coefficients. Their unrounded
weighted-transpose relation does not imply exact native work, continuity or
partition of unity after rounding. Diagnostics report represented-grid versus
physical particle charge and stored-norm versus declared geometry separately.

## Field controller, work and limits

Every new-family RK affine row uses the canonical exact original-duration
expression and one final T round (`:exact_original_duration_rk4_v1`). Ordinary
legacy RK/gather/drift bits need not match this new map; old methods are unchanged.
Intermediate and final magnetic candidates are explicitly projected. Each
segment compares a full map with two half maps using an exact stored-value
infinity difference; accept the two-half result without extrapolation or a
claimed global-error bound. Failed owned trials restart the entire duration at
the next bounded dyadic level. No partial external duration publishes.

At each attempted level, exact segment sums and exact halving must preserve h;
the actual half-map duration passes `_require_resolved_dt`, requiring nonzero
half>=2floatmin(T). This is a declared admission policy, not a claim that exact
arithmetic cannot operate below it. Only own error/nonfinite affine trials and
named projection failures trigger retry. Canonical source/materialized-range,
storage and other errors propagate. Source-risk Faraday differentiation remains
explicitly unsupported; finite completed local E does not settle that source.
See [mixed operator qualifications](MIXED_HALL_OPERATORS.md).

A level budget L permits at most `2^(L+1)-1` trial segments and 12 RHS/rate/field
projections per segment, each with bounded projection iterations and fresh gates.
`max_draws` counts stock uniform requests, including canonical normal/negative-
drift rejection routes. Exhaustion rejects the whole private stream candidate;
there is no truncated distribution or reseeding. Element counts are bounded;
stock-normal recursion can use O(max_draws) stack and exact scalar byte sizes
depend on values/history. No measured speed, allocation or wall-time claim is made.

`open_hall_diagnostics` returns detached exact source/particle scalars, finite-
rounding field/SAT/projection/RK work defects and fixed-log-bit energy intervals. Signed endpoint logarithm coefficients merge
before bounding the energy change, retaining exact cancellation of shared terms.
Last and cumulative per-species/face birth/death mass, charge, momentum and
kinetic exchanges check retained records; kick energy-minus-work defects are
reported separately.
The electron boundary enthalpy is an endpoint trapezoid estimate using the
complete represented `ui_x-Jx*ninv` source. The whole-energy residual includes
unresolved electron-reservoir/coupling/spatial/time terms; it is not silently
identified with zero or with numerical error alone. Stage projection changes are
reported separately and are not added twice to endpoint magnetic exchanges.
