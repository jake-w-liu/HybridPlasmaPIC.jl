# Closed Yee particle charge, current and path work

`YeeChargeWorkspace`, `yee_particle_charge!` and `yee_longitudinal_impulse!`
provide the CPU Float32/Float64 charge/current prerequisite of CAP-002B on a
closed one-dimensional PEC interval. They do not select or advance a particle
trajectory. Coupled bounded electromagnetic PIC remains incomplete.

`yee_longitudinal_voltage(Ex,xold,xnew,g)` supplies the read-only scalar path
work prerequisite on the same closed interval. It does not choose a force or
advance particles.

An optional `YeeTriangularCloud(a)` shape now supplies static finite physical-width
charge/current/path work and a stationary work-adjoint field. Its source review is complete; essential regressions and exact-revision GitHub
acceptance remain pending.
It is a separate declared cloud model; point/CIC remains the default.

For N cells, nodal charge uses the CIC hats at the exact stored-input coordinate
`s=N*x/L`, including both wall nodes. It is extensive charge, without division
by dx. The longitudinal impulse integrates the oriented path overlap of the
same hats between caller-supplied endpoints, using the original species charge
and particle weight. It has N entries and requires no additional dt or dx factor.
In exact arithmetic, `Qnew-Qold+[I[1],diff(I)...,-I[end]]=0`. Independently rounded
charge and impulse have a normwise continuity residual, not bitwise continuity.

Accept one `ParticleSet`, a tuple or a one-based vector of species. Abstract
collection/component declarations and strided views remain supported; consumed
coordinates and weights have grid precision. Charge and weights must be finite,
weights nonnegative and endpoints in `[0,L]`. Mass, velocity and provenance
values are not operands, although all component shapes/storage are checked.
Outputs are ordinary vectors of grid precision. The output and both workspace
arrays must be independent of one another and all inputs; read-only inputs may
share storage. One workspace belongs to one call at a time.

All populations, endpoints, shapes and aliases are checked before scratch writes.
The default point/CIC native accumulation has bound `gamma*S`, where
`gamma=(P+10)*u/(1-(P+10)*u)`,
`u=2^-53`, P is particle count and S is the sum of absolute contributions to the
entry. Final conversion adds `eps(T)/2*(abs(exact)+gamma*S)` plus half a subnormal
quantum. Native arithmetic requires nearest hardware rounding, gradual underflow
and the count/coefficient/product/range checks. Otherwise the same kernel replays
exact rational original operands, including zero, order and geometry decisions.
Final exact conversion is nearest-even independently of hardware context. A complete
entry whose canonical rounding is nonfinite raises `ArgumentError` before output
publication. Scratch contents can change on arithmetic rejection and remain
reusable; invalid-input preflight preserves scratch too.

## Longitudinal path voltage

Let `C[p,e]` be the oriented normalized cell overlap already used by the impulse,
and let `V[p]=dx*sum(Ex[e]*C[p,e])`. The exact contraction identity is
`sum(q[p]*weight[p]*V[p])=dx*sum(Ex[e]*I[e])`. Supplying the exact temporal
mean `(Eold+Enew)/2` makes this the longitudinal part of the existing field
source-work expression. Constructing that mean, rounding the impulse and
rounding each voltage introduce their own errors; no bitwise equality between
independently rounded APIs is promised.

The Maxwell work metric uses stored `dx`, while endpoint geometry uses exact
`N*x/L`. For nonzero physical displacement, the derived work-adjoint average
is `V/(xnew-xold)`, including `N*dx/L`. A constant field therefore gives
`V=(N*dx/L)*Ex*(xnew-xold)` before rounding. The factor differs from one when
stored `L/N` is rounded; the operator preserves this existing discrete metric.
It includes longitudinal DC without imposing a fixed voltage or free plate
charge. Wall endpoints remain geometric endpoints, not absorption events.

The field is a CPU strided vector of grid precision with N finite entries.
Both scalar endpoints have that precision and lie in `[0,L]`. Validation
precedes calculation, including for a zero path. The complete expression uses
exact stored operands and rounds once through the shared nearest-even IEEE
owner, independently of hardware rounding/flush modes. A nonfinite final
rounding raises `ArgumentError`; inputs remain unchanged on every outcome.
A zero path returns zero voltage and does not select a stationary-node force.

Charge/impulse and voltage share the canonical CIC path stencil. Path voltage
scans N field entries for validation, then uses O(crossed cells) rational
operations and scalar temporaries, without a new workspace. A native fast path,
allocation and latency acceptance remain open.
This work contraction does not resolve the retained first-order CIC point-force
or discontinuous node-at-rest failures. A coupled force/pusher needs its own
shape, time-centering, stability and wall contracts.

The charge/impulse workspace payload is `(8+sizeof(T))*(N+1)` bytes plus object headers. Charge costs
O(N+P); impulse costs O(N+P+crossed cells), at worst O(N+P*N). Exact replay owns
O(N) temporary rational entries. No universal cold allocation or latency bound
is claimed.

The private native/exact stencil signatures express the admitted IEEE scalar
domain without converting positions or restricting public containers. Earlier
isolated physical/resource checks and the failed broad-input compiler diagnostics
remain retained in the external audit archive. New essential regressions and
concrete compiler checks run through GitHub CI; changed-source acceptance is
pending until those jobs pass.

## Optional physical-width triangular cloud

```julia
cloud = YeeTriangularCloud(0.05) # physical half-width, fixed under refinement
yee_particle_charge!(charge, particles, grid, workspace; shape=cloud)
yee_longitudinal_impulse!(impulse, particles, xold, xnew, grid, workspace; shape=cloud)
voltage = yee_longitudinal_voltage(Ex, xold[1], xnew[1], grid; shape=cloud)
stationary = yee_longitudinal_field(Ex, xold[1], grid; shape=cloud)
```

The scalar calls take one coordinate, whereas the impulse takes one endpoint
vector per supplied population. This example's scalar coordinates refer to a
single `ParticleSet`. The shape retains a finite positive Float32/Float64 width
at its original precision even if the grid has another precision. It is a
standalone Yee shape, with no periodic deposit/gather or GPU implementation.
The stationary helper requires explicit cloud selection and returns a field,
without q/m, acceleration or a trajectory.

The density is `S_a(t)=max(a-|t|,0)/a²`. Convolve the closed point/CIC hats with
this density, extending the left/right endpoint hat constantly outside its
nearest wall. Outside cloud tails remain **particle-bound endpoint support**;
they are not absorbed particles, free plate charge or an injection ledger.
All centers in `[0,L]`, coincident/overlapping clouds and every mesh node are
admitted. This static convention does not define a finite-radius mechanical
wall, reflection or absorption.

For exact `h=L/N`, let `P_a` be the density's cumulative and `R_a'=P_a`. Each edge
has one canonical cumulative and derivative:

```
A_e(x) = (R_a(x-e*h)-R_a(x-(e+1)*h))/h
g_e(x) = A'_e(x) = integral_cell_e S_a(r-x) dr/h
C(x) = [1-A_0, A_0-A_1, ..., A_last]
I_e = sum(q*w*(A_e(xnew)-A_e(xold)))
V = dx*sum(Ex*(A(xnew)-A(xold)))
G_E(x) = dx*sum(Ex*g(x))
```

`C` is nonnegative and sums to one. Exact continuity and work contraction follow
from these same cumulative differences. `g` is C1 through mesh nodes, so the
stationary field has one continuous value there. At a zero-field rest equilibrium,
the linear coupling matrix is `K_ef=dx*sum(q²*w*g_e*g_f/m)`, symmetric positive
semidefinite for positive masses and nonnegative weights. This is a spatial
coupling property, not nonlinear/time-step stability or an accepted mover.

Interior DC response retains `eta=N*dx/L`. Near a wall it is
`eta*Ex*integral_0^L S_a(r-x)dr`; at either wall with `a<=L` it is half interior
DC. The continuum `Ex=r` cloud force integral at the left wall is `a/6` when
`a<=L`; for `a>L` it is `L²/(2a)-L³/(3a²)`.
The discrete helper uses the actual edge samples and its declared work metric;
it does not return the unchanged point-wall target zero. Holding a fixed is
essential: the source-derived cell-average force bound is O(h²) with a constant
depending on a. Taking `a=O(h)` does not establish fixed-macrocharge point-force
convergence. Point-sheet/contact and point-wall rejection evidence stays open.

The initial cloud tier evaluates complete stored-input rational polynomials and
whole-entry sums, then rounds once using the shared nearest-even IEEE owner.
It includes original q and weight cancellation and is independent of hardware
rounding/FTZ. Invalid input preserves output and scratch; unrepresentable required
output rejects before publication, with reusable scratch. No point native error
bound is reused. Exact rational support bounds omit edges with identically
zero contributions. With k_p retained edges for particle p, charge/impulse use
O(N+sum(k_p)) rational operations, with worst case O(N+P*N). Fixed physical
width can still span an increasing number of cells under refinement. Scalar
field/voltage validate all N entries and contract only the retained support.
The workspace payload is unchanged. Native arithmetic, allocation and latency
acceptance remain open; no measured speedup is claimed.

The compact new regressions use independent polynomial volume/path quadrature,
whole original-source cancellation, mixed width/grid precision, broad/view
storage, continuity, node/overlap/wall/DC response, stationary Gram pairing,
fixed-width cell-average force bounds and rejected publication/context controls.
They have not run locally. GitHub runtime, formatting and full compiler review
must pass on the integrated revision before this prerequisite is accepted.

## Event-free longitudinal coupled source candidate

`YeeLongitudinalPIC` adds a synchronized CPU Float32/64 longitudinal
**electrostatic** model. Source implementation and essential tests are present;
static-cloud and coupled exact-revision GitHub runtime, compiler and formatting
acceptance remain pending. This does not accept full bounded electromagnetic
PIC, point-sheet/contact behavior or wall reflection.

```julia
g = YeeGrid((32,), (1.0,))
particles = ParticleSet{1,Float64}(1; q=-1, m=1)
particles.x[1] .= 0.5
particles.weight .= 0.1
state = YeeLongitudinalPIC(g; shape=YeeTriangularCloud(0.1),
    background_charge=0.1, left_plate_charge=0.02)
init!(state, particles)
step!(state, particles, 0.001)
diagnostics = yee_longitudinal_diagnostics(state, particles)
save_checkpoint("longitudinal.bin", state, particles)
```

One `ParticleSet` or an ordered one-based vector is admitted, including empty
populations and empty collections. All populations move. Charges are finite and
signed, masses positive and weights finite/nonnegative. Zero-weight tracers move
without contributing source or weighted kinetic energy. Positions/velocities
have grid precision and independent writable ordinary or strided CPU storage. Centers
must lie strictly inside `(0,L)`; cloud tails may overlap the walls. All
transverse velocities, transverse E and all B must have stored zero values.
Private particle candidates own ordinary vector storage through the shared PIC
staging owner; public particle objects, views and component identities remain bound.
All particle leaves, including frozen weights and provenance, must support
ordinary Julia writes and expose aliases: checkpoint restore writes every leaf.
Shared preflight rejects `Base.Experimental.Const` and `AtomicMemory` ancestry
through standard views, reshapes and reinterpretations. Writable views and
ordinary reinterpretations retain their admission. This is a bounded storage
check, not a universal writability trait; OS protection changes, arbitrary custom
setter failures and concurrent mutation are outside the rejection guarantee.
Initialization, stepping and restart require nearest hardware rounding and
gradual underflow; other contexts reject before staging. Static cloud APIs keep
their independent exact-arithmetic context contract.

The constructor requires explicit cloud selection. Original Float32/64 width,
background density `b` and left isolated free-plate charge `sigma` retain their
own precision and bits; other real background/plate arguments resolve once to
grid precision. With exact `h=L/N`, initialization evaluates each complete
original-source expression and rounds once:

```
E_e = sigma + b*(e+1/2)*h + sum(q*w*(1-A_e(x)))
right_free_plate = -sigma-b*L-sum(q*w)
```

The derived right free plate is retained exactly as `Rational{BigInt}`, even
when it is outside grid floating-point range. Background endpoint node charges
are `b*h/2`; interior nodes have `b*h`. Free plates and particle-bound endpoint
tails remain distinct. Initialization reuses the exact support-limited charge
kernel and an exact nodal prefix before final field rounding; nodal charge need
not fit the grid scalar range. It uses O(N+sum(k_p)) rational operations and O(N)
transient exact storage, with worst-case O(N+P*N). These are source-level costs,
not measured allocation or latency acceptance. No DC removal, voltage condition
or Gauss projection is applied. Required initialized fields must round finite.
Successful `init!` resets time, step and last-step diagnostics and binds particle objects, component
identities, q/m, weights, IDs and tags. Changing that configuration requires
explicit reinitialization. Particle positions, longitudinal velocities and Ex
are evolving state; callers remain responsible for any deliberate external
changes to them.

`step!` uses the existing guarded original-factor B=0 kick, checked drift,
canonical cloud impulse and sourced Yee row arithmetic:

```
u = v + (dt/2)*(q/m)*G(E,x)
y = x + dt*u
I = sum(q*w*(A(y)-A(x)))
Enew = E-I
vnew = u + (dt/2)*(q/m)*G(Enew,y)
```

Each kick retains the complete rational stationary gather, original q/m and
original resolved full-dt half through the canonical Boris expression. Only
the final velocity rounds to grid precision, nearest-even. A gathered field
need not fit grid precision if the completed kick is finite. This source
completion changes the numerical map; ordinary native velocity bits are not
preserved. Static public field/voltage APIs retain their once-rounded finite
output contract. Runtime/compiler/format acceptance remains pending.

The first successful step freezes the positive resolved dt. Each kick retains
its exact original half duration, including an unrepresentable half. x, v and E
remain at one integer-time phase. No Maxwell
clock or transverse CFL is introduced. Before staging, exact original-operand
arithmetic enforces both strict scoped guards:

```
dt^2*stored_dx*min(1/h,1/a)/h * sum(q^2*w/m) < 4
Qsum = sum(abs(q*w)); Ebound = max(abs(E)) + Qsum
Bparticle = abs(q/m)*(stored_dx/h)*Ebound
abs(v)*dt + Bparticle*dt^2/2 < min(x,L-x)
```

The first guard bounds zero-field rest modes, not general nonlinear stability
or accuracy. The second is a conservative exclusion of a first physical wall
exit, not a temporal error bound. Rounded drift endpoints must also stay
strictly interior. Risk/contact rejection is atomic; there is no reflection,
clipping, wrapping or endpoint nudge. All particle/field/clock/diagnostic
candidates are completed before public writes, retaining public array identities.
Scratch may change on rejection and is reusable.

`em_field_energy(state)` uses `stored_dx*sum(Ex^2)/2`. The coupled diagnostic
uses the existing kinetic and field-energy owners and computes each complete
signed original-source Gauss residual before one final rounding, including
background and both free plates. Read-only energy/residual outputs can be
infinite when their complete value exceeds scalar range. Both old/new kinetic
diagnostic evaluations and the stored last-step impulse, source work and energy
defect must be finite for a step to publish.
`last_work[]` is positive for transfer out of the field and uses the old/new
fields directly. `last_energy_defect[] = Knew-Kold-last_work[]` uses the existing
rounded kinetic-energy diagnostics. Exact-arithmetic Verlet generally has a
nonzero physical `DeltaH=DeltaK-work`; in floating point, field-source staging
and diagnostic rounding also separate this defect from measured `DeltaH`.
Neither exact energy nor bitwise Gauss conservation is claimed.

Distinct checkpoint schema 2 saves initialized completed boundaries through
`save_checkpoint` and `archive_run`: original configuration types/bits, closed
grid including stored dx, policies, ordered phase space/provenance/q/m/weights,
current fields, clock/fixed-step phase, exact derived right plate and all exposed
last-step diagnostics. Restore validates and detaches the complete payload;
current Ex is restored verbatim, without projection or reconstruction. Scalar
and collection payloads cannot be interchanged. Matching configuration is
required; saved particle values replace destinations and rebind their existing
component identities. Ordinary vectors may resize; strided destinations must
already have matching counts. Scratch is rebuilt. Archives reuse the existing
metadata/checksum owner and checkpoints reuse atomic file publication.
Both coupled schemas require `force_policy=:original_geometry_boris_nearest_even`
and integer `force_policy_version=1`. Schema 1 belongs to the former
rounded-gather/native map and is rejected before restore mutation; no silent
migration or legacy mode is provided. The archive container remains schema 4:
`load_archive` checks integrity, while model restoration checks numerical policy.

The `longitudinal-pic` implementation-essentials suite includes independent
physical-volume initializer/time-map references, continuity/energy checks,
nonzero oscillator temporal refinement, exact admission boundaries, source and
alias failures, late diagnostic rollback, mixed shape/grid widths, and restart
schema/continuation cases. These are pending GitHub execution; no local Julia,
formatter, compiler, build or profiling result is claimed. Complete-step resource
costs and broader nonlinear/long-time physical validation remain separate gates.

## Event-free 1D3V electromagnetic source candidate

`YeeElectromagneticPIC` adds both transverse polarizations and all three velocity
components on the closed 1D PEC grid. Source implementation and compact essential
tests are present; exact-revision GitHub runtime, compiler and formatting acceptance
remain pending. This finite-width event-free scope does not resolve retained
point-sheet/contact/wall-event failures or provide broader nonlinear stability.

```julia
state = YeeElectromagneticPIC(YeeGrid((32,), (1.0,));
    shape=YeeTriangularCloud(0.1), c=1.0,
    background_charge=0.1, left_plate_charge=0.02)
# Set finite transverse E/B, with Ey/Ez endpoints and all Bx stored zero.
init!(state, particles)
step!(state, particles, 0.001)
diagnostics = yee_electromagnetic_diagnostics(state, particles)
```

The population, source/provenance/storage binding, original width/b/sigma,
free-plate, strict-interior, nearest/gradual and fixed-after-first-success dt
contracts follow the longitudinal model. c resolves once to grid precision.
All public x/v/E/B share integer time. Ex/By/Bz use N cell/edge entries;
Ey/Ez/Bx use N+1 nodes. All Bx must be stored zero, consistent with the 1D
solenoidal cavity and homogeneous normal-B wall condition. Ey/Ez endpoint zeros
are required, never projected. There is no normal guide flux, wall event,
periodic/multi-D, relativistic, MPI/device, circuit or injection mode.

With exact h=L/N and stored delta=dx, retain the canonical A/g and nodal C:

```
G1(F,x) = delta*sum(F_e*g_e(x))
G0(F,x) = sum(F_i*C_i(x))
GE = (G1(Ex), G0(Ey), G0(Ez)); GB = (0, G1(By), G1(Bz))
Cbar_i(x,y) = integral_x^y C_i(s) ds/(y-x), or C_i(x) for a stored zero path
Ix_e = sum(q*w*(A_e(y)-A_e(x)))
Iy/z_i = sum(q*w*dt*u_y/z*Cbar_i)/delta, for interior nodes
Iy/z endpoints = 0
```

G0 has no extra delta/h factor or wall renormalization. The derivative identity
`d_x G0(E)=G1(Ddelta E)` follows from the same cumulative. Transverse current
endpoints are absent PEC Dirichlet degrees of freedom: their E is zero, so
masking them preserves the exact work contraction without removing endpoint
particle-bound charge. A private quartic antiderivative of the original
triangular density supplies exact path averages. Original q/weight/duration/
velocity/metric factors and complete signed current sums round only at final
entry publication. Required final impulses must be finite. GE/GB remain complete
rational contractions through canonical Boris with original q/m and exact
resolved full-dt half. Only each half-kick's final three velocity components
round to grid precision, nearest-even; an unrepresentable gathered component is
admitted if that completed velocity is finite. The first rounded velocity still
owns drift/current and the second kick. Original factors are prepared once per
step/population, with no per-particle/grid matrices. Ordinary native bits are
not preserved. This source completion awaits GitHub runtime/compiler/format
acceptance; rational-kick latency, allocations and peak memory are unmeasured.
Static point/CIC and existing cloud A/g arithmetic remain unchanged.

The explicit map is a guarded full-vector Boris half-kick at old x/E/B,
a checked straight drift using its velocity u, compatible currents, one canonical
B-half/E-full-minus-current/B-half field candidate, then a Boris half-kick at
new x/E/B. Exact original half duration is retained when stored dt/2 is unresolved.
The ideal map is self-adjoint and consistent, giving formal second order away
from events; floating rounding and runtime order remain separate verification.
Zero-transverse states reduce algebraically to longitudinal Verlet for steps
admitted by both models; the EM admission restrictions remain.

A step must be strictly below the existing downward-rounded Maxwell CFL. It
also enforces the existing longitudinal rest bound and the combined transverse
zero-B-rest guard `dt²*(4*c²/delta²+sum(q²*w/m)/delta)<4`. This follows from the
actual composed rest operator `c²*Ddelta'*Ddelta+sum(q²*w*Cint*Cint'/(m*delta))`.
Neither independent-axis guards nor this scoped rest bound prove general
nonlinear or finite-guide-field stability, positivity or accuracy.

The continuous declared cloud semidiscretization has the work-adjoint physical
energy `H=K+delta*sum(E²+c²B²)/2`. It bounds the gathered electric norm by
`E*=max(1,delta/h)*R`, where exact rational R encloses `sqrt(2H/delta)` upward.
The strict full-speed tube
`dt*sum(abs(v))+dt²*abs(q/m)*E*/2 < min(x,L-x)` excludes its first physical exit,
including transverse rotation and zero-weight tracers. It uses original rational
energy values, not an overflowing rounded energy. Separately, every rounded
straight-drift endpoint must remain interior. This tube does not enclose a true
continuum Maxwell trajectory or certify temporal accuracy. No clipping or nudge
repairs rejected paths.

`em_field_energy` returns physical field energy. `yee_modified_energy(state,dt)`
uses the existing fixed-step modified field quadratic; it is distinct from
physical energy. Diagnostics expose both plus synchronized kinetic energy and
complete signed Gauss residuals with plates. `last_work` is positive out of the
field; `last_modified_energy_defect` is rounded-diagnostic DeltaK-work, and
`last_physical_energy_change` is one complete original-input physical H difference
between stored states. None is forced to zero. Read-only diagnostic overflow may
return infinity; stepping requires finite old/new rounded kinetic diagnostics,
source work and stored last-step values before commit.
No exact fully discrete energy or bitwise Gauss conservation is claimed.

All candidates, late diagnostics and clocks stage before public writes. For admitted
writable host storage, rejection preserves particles, fields, histories, array
identities and clock; scratch may change. One instance must not be mutated/stepped concurrently. The distinct EM
schema 2 saves original config types/bits, resolved c, policies/versions, ordered
particles/provenance, current E/B verbatim, synchronized fixed phase and histories.
Restore validates/detaches everything before publication, with ordinary particle
arrays resizable from empty storage and matched strided views retaining identity.
Scalar/collection layouts and the longitudinal model's schema remain distinct.
The required shared force policy and schema-1 rejection are described above;
the EM spatial gather/current/PEC policies keep version 1.

The `yee-electromagnetic-pic` essentials suite includes independent physical
quadrature/current and matrix time-map/rest references, nonzero temporal/energy
refinement, original-factor/half-duration and late-error controls, zero-transverse
reduction, alias/context rejection and complete checkpoint/archive continuation.
Compact independent original-gather underflow/overflow, signed cancellation,
magnetic overflow rescue and final-half rejection/retry controls cover the new
force policy in both coupled models.
These are source tests awaiting GitHub. No runtime, accuracy, allocation, latency
or peak-memory acceptance is claimed. Wall histories/events, additional boundary
models, broader physical validation and parallel/device paths remain open.

### Complete-step cost observation

The existing manual `yee-electromagnetic-pic` essentials suite runs
`benchmark/yee_coupled_steps.jl` before its unchanged regression assertions and
retains the output in `implementation-essentials.log`. The script observes 13
small scalar profiles: Float32/64 longitudinal and EM B=0/B!=0 ordinary states,
plus the existing compact original-force range tracers. Three fresh matched
samples each measure two complete public steps after one unmeasured step fixes
the phase. Initialization, checkpoint-state checks, explicit GC and logging are
outside timing; all production validation/current/field/diagnostic work stays in.

Raw Base `@timed` seconds, allocation traffic, GC/compilation counters and actual
Julia/source/CPU/thread metadata provide a baseline for future paired work.
Allocation traffic is not live or peak memory. The observation has no resource
threshold, numerical-map change or acceptance claim; its first remote results
remain pending. Larger workloads, collection/view/backend scaling and peak
memory remain unmeasured. Regression and compiler assertions retain their gates.

## Planned event-local reflection

Method design reviewed 2026-11-01; carried-error design reviewed 2026-11-02.
**The wall model has a source candidate; runtime acceptance remains pending.**
[`YeeImpactLongitudinalPIC`](YEE_IMPACT_LONGITUDINAL.md) is a distinct longitudinal CPU Float32/64 model with
at most one common contact cluster and two synchronized local KDK
(kick-drift-kick) chunks per external step. Existing event-free models keep their methods and restart tags.
Reflection acts on the cloud center; the current particle-bound cloud tails,
isolated free plates and fixed source counts remain. This does not supply a
material sheath, mirrored cloud, absorption, emission or transverse wall model.

Complete the duration-s KDK, including its final electric half-kick, before
reversing the incoming velocity. Then advance the exact remainder `tau-s` from
that synchronized state. Internal durations and halves retain original rational
values; conversion to a stored timestep or a dummy Maxwell duration is not
allowed. Reflecting the initial whole-step half-kick has a first-order impact
velocity error for an off-midpoint hit under constant acceleration: the exact
one-hit velocity is `-v+a*(tau-2*tc)`, while that split gives `-v` when both endpoint
accelerations equal a. This algebra rejects that split as the proposed map; it
is not a measured convergence result. Event-local force centering has a
conditional second-order route for isolated transversal impacts away from the representation floor;
independent convergence evidence remains required.

For every particle, examine both walls and every monotone branch of its
quadratic predictor, including interior extrema. Endpoint-in-box checks alone
miss exit-and-return trajectories. Exact rational brackets and finite refinement
limits must identify the earliest outgoing branch. Every member of a common
cluster must independently pass endpoint, timing, source-sensitivity and actual
completed incoming-velocity gates; overlapping time intervals alone do not prove
simultaneity. Grazing, unresolved ordering, terminal contacts, a second cluster
or exhausted refinement reject the complete external candidate.

At a represented contact, constrain the private endpoint to the exact wall and
report the correction `r = wall - (x+s*u)`. Its proposed bound is four ulps of
`|x|+|s*u|`, evaluated without first rounding that scale. If `u_star=v+s*a0/2`,
the exact predictor residual is `r+s*(u-u_star)`; `r` alone does not bound the time
error. Positive lower incoming-speed bounds and position/velocity/source-force
sensitivity bounds must also admit the common time. This is an explicit
constrained endpoint policy, not an unreported coordinate nudge.

The canonical endpoint-current owner must deposit the actual constrained
positions. Keep each accepted local current/field work contraction exact until
the signed external sum is rounded once. Sum net-current history separately;
it is not either stage current. Reflection leaves position, charge and field
unchanged, transfers momentum from the completed incoming velocity, and has
zero physical wall work and charge transfer. Track represented electric-kick
momentum, wall impulse, all-node Gauss, physical energy change and generally
nonzero particle-work defect without adding artificial wall heat. Stage the
whole external operation, histories, clock and ledgers before public stores.
A distinct versioned checkpoint must identify the method, endpoint policy,
cluster scope and complete wall/history state.

The first-contact physical tube can be derived from the existing triangular
cloud and field-current owners. After an approximate event, a fresh tube from
the numerical reflected state only bounds flow initialized at that state.
**Original external-start trajectory coverage has candidate source enforcement;
exact-revision runtime and independent physical acceptance remain pending.**
Carried position, velocity, acceleration, field-invariant and event-time errors
must cover time alignment, reflection and the enlarged remaining horizon before
this candidate can claim no hidden second physical contact. Refining a predictor
root does not remove physical truncation error. This proof is a prerequisite
for the stronger wall contract, not a runtime result or an impossibility claim.

- [x] Specify the one-cluster map, rational durations, represented-contact
  policy and local physical-tube bounds; independently review their limits.
- [x] Review the conditional [carried-error contract](BOUNDED_EVENT_COVERAGE_CONTRACT.md)
  through the event and remainder, retaining state, invariant and timing errors.
  Its original-start coverage premises have candidate source enforcement;
  exact-revision remote and independent physical gates remain open.
- [x] Source extraction of the canonical duration-free exact longitudinal row
  and unrounded work value, preserving existing methods by literal expansion.
  Direct runtime/inference/format acceptance remains pending on GitHub.
- [x] Extract the canonical raw charge/current sum and its existing candidate
  conversion so the event stages retain Istar without lifting rounded current.
  Existing public/current routing remains; GitHub acceptance is pending.
- [x] Source candidate for staged reflection, carried certificates, diagnostics
  and distinct restart policies; exact-revision GitHub acceptance remains open.
- [ ] On GitHub, run essential independent force-driven impact refinement,
  branch/cluster failures, complete rollback and restart controls. Preserve
  existing event-free and compiler assertions; measure complete-step cost.
- [x] Separate source candidate for one isolated passive full-cloud capture,
  including remaining charge-to-plate current, dynamic plates and work ledgers;
  [scope and pending acceptance](#one-isolated-passive-capture-source-candidate).

## One isolated passive capture source candidate

[`YeeAbsorbingLongitudinalPIC`](YEE_ABSORBING_LONGITUDINAL.md) implements a distinct
nonrelativistic CPU Float32/64 closed 1D finite-cloud source candidate: zero or one
uniquely isolated certified-passive instantaneous capture per external step.
Complete the rational-duration KDK final half-kick with original forces before
retiring the entire cloud, including particle-bound endpoint tails. Its additional
`I_death,* = q*w*(Bside-A(wall))` current transfers full signed charge to the
selected isolated plate through the canonical source row. Both current free
plates and separate captured kinetic/death-work/received-energy ledgers are exact.
A negative received total or inconclusive physical-passivity lower bound rejects.

Absorption has a separate carry across the finite force jump. It aligns survivor
uncertainty at the true predeath hit, applies the same death map, recomputes the
postdeath gathers and checks both walls over the enlarged remainder. Smooth-epoch
bounds are not extended across the jump. The [reflection contract](BOUNDED_EVENT_COVERAGE_CONTRACT.md)
retains its unchanged fixed-source/velocity-reversal premises and departure rule;
those premises do not cover source removal.

Stable compaction requires ordinary Vector leaves and preserves particle objects,
array identities and surviving provenance, including empty ordered populations.
Distinct schema-1 checkpoint/archive state includes dynamic plates, initial source
records and capture/certificate histories; no implicit cross-model migration is
provided. See the model guide for binding, transaction and restore limits.

GitHub numerical/convergence, rollback/restart, JET/compiler, formatting/CI and
complete-step resource acceptance are deferred during the current
implementation-first pass; existing authored controls remain unrun. General
material/sheath/circuit/emission laws, multiple captures, transverse walls,
relativistic/device/distributed execution and coupled absorbing field layers
remain open. The separate [scalar-conductivity field layer](YEE_CONDUCTIVE_LAYER.md)
now has source evolution, exact budget and restart; attaching it to particle
charge/current requires its own sink/injection derivation.
This source candidate does not complete CAP-002B or establish whole-run stability.

## Implementation tracker

The optional static-cloud API, surface-tail model and canonical owners above
were recorded before production edits and reviewed on 2026-11-01. The source
implementation preserves point defaults; exact-revision runtime/compiler/format
acceptance remains a separate gate.

- [x] Closed extensive charge and compatible longitudinal endpoint impulse.
- [x] Canonical longitudinal path voltage and exact discrete work contraction.
- [x] Optional explicitly named finite physical-width triangular-cloud static
  charge/current/voltage and stationary work-adjoint field source implementation.
- [ ] Accept the optional cloud's independent essential regressions, full compiler
  review and formatting on the exact integrated GitHub revision.
- [ ] GitHub validation of the path-voltage implementation and unchanged impulse.
- [ ] Choose a force/current pairing that passes retained point-sheet force,
  isolated-plate/DC, work and semidiscrete stability controls.
- [x] Source candidate for exact PEC-masked transverse path currents and a
  synchronized event-free 1D3V Boris/Yee method.
- [ ] Accept the EM candidate on the exact integrated GitHub revision.
- [x] Source completion of original-geometry-through-Lorentz semantics for the two
  event-free coupled models, retaining exact gathers/factors through canonical Boris and
  rounding only each final half-kick velocity; schema-2 force policy is explicit.
- [ ] Accept the original-force map, schema transition and whole-step resource
  behavior on the exact integrated GitHub revision.
- [x] Extract private duration-free exact source-row and unrounded rational
  source-work owners; existing native, exact and public expressions are preserved
  by helper expansion. The reflection source candidate consumes these owners.
- [ ] Accept direct signed-zero/range/work and unchanged-caller controls on the
  exact integrated GitHub revision, including their event consumer.
- [x] Source candidate for one-cluster longitudinal reflection with branch and
  carried-error checks, atomic failure and wall momentum/charge/work ledgers.
- [ ] Accept independent impact convergence, branch/common-time rejection and
  rollback/restart on GitHub; broader repeated-impact dynamics remain open.
- [x] Source candidate for event-free longitudinal Verlet coupling, energy/Gauss
  diagnostics, uniform background and exact isolated free plates.
- [ ] Accept the coupled candidate on the exact integrated GitHub revision.
- [x] Source candidate for synchronized EM physical/modified-energy and Gauss
  diagnostics with distinct complete checkpoint/archive.
- [ ] Accept EM diagnostics/restart; transverse wall/plate events remain open.
- [x] Distinct source checkpoint/archive schema for initialized event-free
  longitudinal fields, ordered particles, fixed phase and exposed diagnostics.
- [x] Distinct source checkpoint/archive for completed longitudinal reflection,
  root budget, wall/electric ledgers and carried error histories.
- [x] Source candidate for one uniquely isolated certified-passive capture with
  full death current, dynamic exact plates, separate ledgers, jump-aware survivor
  carry, stable Vector compaction and distinct checkpoint/archive.
- [ ] Accept event-free longitudinal, EM, reflection and capture restart controls,
  plus independent capture passivity/force-jump/rollback and full compiler checks
  on the exact integrated GitHub revision; measure complete-step resources.
- [ ] General material absorption/injection/circuits, multiple captures,
  relativistic coupling, higher dimensions, GPU/MPI and broader physical validation.

The coupled source candidates require their exact-revision GitHub acceptance;
the reflecting utility and prescribed-source field step retain separate contracts.
Retained force and impact-convergence failures remain open.
