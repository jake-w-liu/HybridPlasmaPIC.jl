# Coupled error control for periodic EMPIC

`EMPIC` and planar `EMPIC1D` accept an explicit `error_control` proposal on an
ordered persistent species schedule. This CPU Float32/Float64 policy evolves
physical positions, velocities, E and B at the same full-time phase, using a
projected Euler–Heun pair and the canonical charge/current/Maxwell operators.
It measures a coupled embedded defect. Scientific order, true-error guarantees,
effectivity across support cuts and resource acceptance require independent
remote results; source implementation alone does not establish them.

## Select the policy

Use `species_dt`, `event_budget` and `error_control` together. The opt-in supports
CIC/TSC; NGP is rejected at this constructor selection. Existing NGP ordinary,
dividing and persistent schedules remain available with their existing phases.
`species_stride` cannot be combined with this policy. No old phase is migrated
or inferred by a setter or restart.

The proposal has exactly these keys, including exactly five named tolerance
roles in the displayed order. There are no hidden common-unit tolerances:

```julia
T = Float64
proposal = (
    pair=:projected_euler_heun_v1, metric=:component_max_v1,
    domain=:cic_tsc_continuous_v1,
    atol=(x=T(1e-5), v=T(1e-5), E=T(1e-5), B=T(1e-5), clock=T(1e-12)),
    rtol=(x=T(1e-3), v=T(1e-3), E=T(1e-3), B=T(1e-3), clock=zero(T)),
    h_min=T(1/65536), h_max=T(1/128), max_trials=1024,
    max_force_evals=100_000, max_deposits=1_000_000,
    max_fft_calls=100_000, max_report_cells=2_000_000,
)
```

This is an unexecuted API illustration; its tolerances, budget and time bounds
must be chosen for the user's units/problem. Atols must remain finite positive
after conversion; rtols are finite in `[0,1)`. Time bounds are positive resolved
T intervals, with `h_min<=h_max`; each budget is a positive Int-sized integer.
The policy requires nearest rounding and gradual underflow.

```julia
using HybridPlasmaPIC
g = FourierGrid((8,), (one(T),))
a = ParticleSet{1,T}(1; q=T(-2), m=T(3))
b = ParticleSet{1,T}(1; q=one(T), m=T(5))
a.x[1][1], a.weight[1], a.v[1][1] = T(7/32), T(1/2), T(1/16)
b.x[1][1], b.weight[1], b.v[1][1] = T(21/32), one(T), T(-1/32)
species = [a,b]
st = EMPIC(g,species; c=one(T),shape=CIC(),
    species_dt=T[1/8,3/16],event_budget=1024,error_control=proposal)
fill!(st.E[2],T(1/32))
init_empic!(st,species)
step_empic!(st,species,T(1/64))
diag = electromagnetic_diagnostics(st,species)
```

Particles must have finite loaded positions and velocities, finite nonnegative
weights, original finite charges and positive masses. Ordinary injective writable
CPU storage and ordered fixed source inventory retain their collection meanings.
Fixed rows require zero velocity and `species_dt=0`; they deposit charge and
stay fixed. Empty rows and zero-weight mobile markers are admitted. Mobile
markers participate in force and error/domain checks regardless of weight.

## The actual coupled pair

`init_empic!` receives full-time x/v/E/B, prepares the canonical source/Gauss
state and starts a zero-time controller/calendar. It performs no velocity or
magnetic half prime. `phase_policy=:full_time_coupled`; physical velocity and
magnetic diagnostic offsets are zero. The planar model keeps its existing
polarization, including stored `v_z=0`.

For each trial, choose an exact source-calendar gap `delta` bounded by the next
mandatory edge and proposed interval. Materialize one resolved actual `h::T`;
BOTH branches consume `H=Q(h)`. Gather/deposit at every mobile origin, then form
Euler x/v/B. Deposit Euler charge and correct its current through the canonical
finite-interval continuity owner before completing Euler E. The actual stored
Euler fields and particles supply the second force/current stage. Construct
Heun from the original state and both stages, deposit its own charge, correct
its actual averaged current and complete its own E/B update. Publish only an
admitted high branch. Euler and Heun normally have different endpoint charge;
they share the source endpoint, actual h, source inventory and physical phase.

Nonrelativistic force retains original q/m in exact scaled increments.
Relativistic force uses physical velocity:

```text
Delta v = H*q/(m*gamma) * [E + v × B - v*(v·E)/c²].
```

Completed Euler/Heun outputs use the shared exact rational/radical combination
and rounding owners. No pre-rounded q/m, unscaled overflow-prone acceleration,
rounded half-h or separately rounded stage increment replaces that formula.
Every materialized relativistic stage/output must satisfy the shared supported
velocity rule: its rounded normalized squared speed stays strictly below one,
with an exact stored-speed check on the admitted near-boundary side. Exact
subluminality alone does not widen that rule. No velocity clipping or
renormalization is supplied.

On effective Fourier range modes, let `G=-i*k/|k|²` and `P=I-G*D`, using the
actual stored `kvec`, including its per-axis Nyquist zeros. The consumed current
is `P*Jraw-G*(rho_end-rho_start)/H`. The longitudinal high-minus-low field
therefore includes `G*(rho_high-rho_low)`. A raw-current difference alone is
not the coupled embedded defect. True DC raw current remains physical; non-DC
all-zero-effective-k current is suppressed. Mixed Nyquist/range modes retain
their valid components. Initial homogeneous excluded E/B modes are retained.
Gauss is a range-mode constraint; physical DC charge and excluded density modes
remain separately observed channels.

## Norm, cuts, clocks and retries

The component maximum includes all mobile spatial positions/all three velocities
and every represented E/B component. Velocity and field scales are
`atol+rtol*max(abs(start),abs(low),abs(high))`; position scale is
`atol_x+rtol_x*L_d`, with the actual unique lifted high-minus-low displacement.
There is no Richardson divisor or assumed convergence exponent.

CIC integer cuts, TSC half-integer cuts and physical periodic seams are admitted.
All stages use their actual canonical rounded stencil, including knot ties;
there is no artificial smoothing or single-piece rejection. Actual branch lifts,
BigInt winding, ideal/represented endpoint difference, sampled stencils and
compact directed cut ranges are retained. Ideal and represented branch travel
must each stay below `L/4`, and branch separation below `L/2`, to give an
unambiguous torus metric. These sampled straight stage paths are not an exact
orbit event history. The physical seam remains explicit when represented
`L` differs from `n*dx`.

The controller separately observes actual-h minus calendar-gap, public rounded
clock increment minus actual h, and their sum under the clock role tolerance.
Source time advances by delta; public time follows actual repeated additions.
Neither clock replaces the other. Nominal species caps/deadlines and exact
base-partition edges remain mandatory. They control future nominal decisions,
not whether a mobile force is sampled: every trial has both force stages even
during a pending nominal interval. Right-edge decisions stay pending until the
next positive call. Zero duration validates and consumes nothing.

Accept only norm<=1 plus finite/source/constraint/speed/travel/role/phase/budget
admission. On an observed retryable rejection discard both branches and retry
half the gap from the identical accepted origin. After success, recommend at
most `2h` when norm<=1/4, otherwise h, within the explicit bounds. These are
controller rules, not an order theorem. Unresolved/below-hmin mandatory gaps or
exhausted budgets fail visibly. There is no promise that every tolerance/macro
endpoint is reachable. Canonical source/FFT/programming errors propagate; they
are not converted by broad exception catches into successful retries.

## Receipts, budgets and persistence

`last_error_control_report` retains its historical control/event proposal,
accepted `intervals`, bounded rejection/deadline receipts, actual origin/low/high
particles/fields/source/current/derivative operands, lifts, norm and work/cost
observations. `error_control` and `controller_next_h` describe future choices;
`nominal_H`, `last_deadline` and `next_deadline` describe the nominal calendar.
`admission_control` retains the policy that admits the current physical state
and report, separately from the future `error_control` proposal.
The full external `step_empic!` increments public step once.

Signed field, particle and coupling defects sum to the actual reported total
energy change. Consumed field work uses actual high current and mean endpoint
E; particle work uses actual gathered electric samples. Generic raw work has
all three current components; planar raw work is explicitly transverse Jy-only.
Relativistic kinetic energy is the canonical reported rounded physical moment,
not an exact irrational number. Rejected trials spend computational budget and
contribute no accepted physical work. Conservation, passivity and stability are
not inferred from these accounting identities.

Budgets count all attempted force evaluations, scalar deposition contributions,
reserved primitive FFT-plan calls including range replays, trials and retained
report scalar cells. Actual admission and final-verification replay costs are
separate receipts and count against each positive operation's work budget.
Valid historical reports use their own retained admission policy for observation
and capture; a smaller future budget can require raising it before advancing.
FFT reservations bound primitive plan calls in the reviewed canonical
dependency; dependency changes require rechecking them. They are not measured
CPU time or peak bytes. Full-operation
allocation/peak-memory/performance acceptance remains separate.

Existing family setters admit copied future proposals. Pair/metric/domain are
immutable; changing tolerances/caps does not retime pending deadlines or rewrite
old report/history. Each old report is validated under its own historical
proposal. A new owner-local guard covers init/step/setter/diagnostics/save/load.
All fallible work, including final cumulative ledger, complete report and shadow,
precedes ordinary publication. A late failure discards the entire external
operation, including previously accepted private subintervals. Independent states
sharing a grid/FFT scratch retain the existing concurrency restriction.

Completed-boundary raw/archive schema 4/policy 3 stores full-time primaries, source,
provenance, exact spectral-operator signature, clocks/calendar, future and
historical controller state, complete bounded report and cumulative ledger.
Original graph/live links are checked before copying. Restore supports changed
ordinary destination counts and stale old numerical values, then publishes the
saved phase verbatim without solving, reprojecting or re-priming. Different
operator bits or old schema 1/2/3 phases do not silently migrate. The canonical
atomic serializer/checksummed archive container remains the IO owner.
Uninitialized count-zero constructor seeds are also retained as uninitialized,
with finite data and empty zero-time history; their source/Gauss preparation
starts only at explicit `init_empic!`. Capture materializes admitted public
strided views into ordinary saved vectors; matching-count view destinations keep
their existing arrays on restore.

Accepted reports undergo exact typed native numerical replay. Reproducible
admission therefore uses the retained native arithmetic/FFT environment and
matching replay results. The operator signature is saved; native plans, wisdom
and library environment are not serialized. Essentials construct fresh target
grids/plans within their runtime. Broader planning/wisdom/backend portability is
unexecuted.

Full NGP jump handling, a true-error guarantee, knot/seam effectivity and
convergence, asynchronous multirate force reduction, source-growth/collision/
boundary/backend intersections and measured resource acceptance remain open.
Minimum/latest supported Julia acceptance runs on GitHub; no local numerical
workload or unexecuted result is presented as a passing check by this guide.
