# Synchronous coupled implicit hybrid

`ImplicitHybridStepper` advances every ion velocity and the magnetic field in one
bounded nonlinear solve. This is a separate integer-time state. Load full-time
particle velocities, seed `st.fields.B`, and call `init!(st, species)` before
`step!(st, species, dt)`. `init!` wraps positions and computes the initial moments
and Ohm fields privately before replacing the public state.

The admitted model is periodic Fourier/CIC in 1D3V, 2D3V or 3D3V, CPU Float32 or
Float64 with nearest rounding and gradual underflow, and an ordinary `Vector` of
ordinary Vector-backed populations. Ion charges and masses must be finite and
positive; weights finite and nonnegative. Every population is mobile. The
isothermal electron temperature uses the grid precision and is finite and
nonnegative; the density floor must remain positive at that precision.
`η = ηH = de2 = 0`. There are no sources, collisions, evolving electron
pressure, subcycles, transport, device or distributed execution in this family.

The stored normalization is the existing hybrid ion/Alfvén convention:
`n = sum(q_s*n_s)`, `ui = sum(q_s*n_s*u_s)/n`, `J = curl(B)`, `pe = Te*n`, and
`E = -ui × B + (J × B)/n - grad(pe)/n`. Magnetic energy has coefficient 1/2.
Species use the canonical Boris map with their original charge and mass. No
Maxwell speed or collision permittivity/count conversion is introduced.

Constructor species count/order, owner object identities and cached charge/mass
are fixed. Ready calls reject changed charge/mass or local particle counts.
Construct a new stepper to change the binding. Reinitialization resets the clock
and solve history while retaining that binding. Same-count external physical
mutation cannot be detected; initialize explicitly at the intended common time
level. Use a separate state/grid for concurrent work. The busy guard rejects
reentry under this nonconcurrent-access premise; it is not a thread-safe lock.

## Represented coupled equations

The unknown contains all new particle velocities `v1` and every component of
`B1`, in population/component/particle order followed by component/grid order.
For each residual evaluation, form complete stored-input means
`vbar = RN((v0+v1)/2)` and `Bbar = RN((B0+B1)/2)`. Construct both physical lifts
`ellmid = x0 + (h/2)*vbar` and `ell1 = x0 + h*vbar` with the canonical guarded
coordinate operation. Both must remain finite and representable. Wrap `ellmid`
for the common deposition/gather stencil; publish wrapped `ell1` only after
success. Wrapped endpoints are never averaged, and no minimum-image path is
chosen. Multicell lifts are admitted; there is no orbit average.

Deposit midpoint charge density and raw charge momentum `Mbar` with the canonical
collection moments and `vbar`. Require every actual midpoint density to be finite
and at least `nfloor` before normalization. The same strict condition applies to
initial and completed endpoint moments. Denominator clipping is therefore
inactive. Canonical Ohm uses `nbar`, `ubar`, `Bbar` and private zero inactive
`lapJ` inputs. Gather its finite `Emid` and `Bbar` at the same midpoint positions.
Let `vstar` be the guarded canonical Boris map of `v0` under these fields for `h`,
and `Fmid = faraday_rhs!(Emid) = -curl(Emid)`.

Freeze `Sv = max(1,maxabs(v0))` and `SB = max(1,maxabs(B0))` for the attempt.
The reference 1 is the normalized velocity/field unit. With `Q` denoting the
existing exact rational conversion, define

```
Gv = RN_T((Q(v1) - Q(vstar)) / Q(Sv))
GB = RN_T((Q(B1) - Q(B0) - Q(h)*Q(Fmid)) / Q(SB))
rho = max(maxabs(Gv), maxabs(GB))
tau = RN_T(Q(atol) + Q(rtol)*Q(rho0))
```

The residual arithmetic retains the original stored operands and rounds once,
including ordinary cases. It does not undo upstream FFT, scatter, normalization,
gather or Boris rounding. It has no measured cost claim. Scales remain fixed;
component errors are never divided by a population or grid count. `rho0` is a
fresh residual at `(v0,B0)`; nonzero `rho0` requires positive finite `tau`.

The velocity equation is equality to the represented Boris map. A small residual
is not an exact-arithmetic Lorentz equation, physical-energy certificate,
unconditional stability guarantee or accuracy estimate for a multicell orbit.
The represented odd Fourier derivatives retain their canonical zero-Nyquist
policy. No new divergence projection or exact discrete continuity law is added.

## Bounded Newton and GMRES policy

The constructor controls and defaults are:

| Keyword | Default | Meaning |
|---|---:|---|
| `nonlinear_rtol` | `sqrt(eps(T))` | Relative part of fixed dimensionless outer tolerance; below 1 |
| `nonlinear_atol` | `64eps(T)` | Absolute part of outer tolerance |
| `max_newton` | 20 | Maximum accepted Newton updates; zero allows initial-only acceptance |
| `max_linear` | 200 | Maximum Arnoldi iterations per GMRES call |
| `memory` | 20 | Restart basis limit, capped at the actual system dimension |
| `max_backtrack` | 12 | Maximum halvings after alpha=1 |
| `max_jv_shrink` | 8 | Maximum halvings of an inadmissible difference probe |

These dimensionless defaults are engineering policy without calibration or a
convergence guarantee. Tolerances must be finite, nonnegative and not both zero;
positive supplied tolerances must survive conversion. Linear/memory limits are
positive Int-sized integers. The other limits may be zero.

Directions are scaled by the fixed block scales. The concrete forward-difference
action is `A(p) = [G(z + epsilon*S*p) - G(z)]/epsilon`, initially using
`epsilon = sqrt(eps(T))*max(1,maxabs(z/S))/norm2(p)`. Physical perturbations and
finite-difference components retain their complete stored operands. `A(0)` is
explicitly zero. A nonzero direction must change at least one represented
physical component. Only the owned expected trial-domain error permits bounded
probe shrink. Unresolved, nonfinite or exhausted probes fail visibly.

Restarted Krylov GMRES uses identity preconditioning, reorthogonalization and an
explicit zero correction. The zero-correction residual is checked first. The
fixed inexact target is `norm2(Gbase)/2`; the internal estimate uses half that
target, which must remain positive. After GMRES, a fresh application of the
actual defined finite-difference action must satisfy
`norm2(A(delta)+Gbase) <= norm2(Gbase)/2`. Krylov status is not acceptance. This
checks the defined action, not its approximation error to a mathematical
Jacobian. CIC stencil transitions provide no differentiability guarantee.

Backtracking tries `alpha = 1, 1/2, ...` within its limit. A fresh admissible trial
must meet `rho_trial <= (1-alpha/4)*rho_base` or already satisfy `tau`. An unchanged
represented iterate fails. There is no hidden timestep retry or explicit
fallback. After apparent convergence, the complete residual is reconstructed
again from the old snapshots and final unknowns. Only then are endpoint moments,
endpoint Ohm fields and the report prepared; endpoint density/range or canonical
report-derivative failure can still reject the attempt.

`ImplicitHybridConvergenceError` records the reason, most recently completed
attempted block norms (an incomplete block is Inf), tolerance, iteration/call
counts and captured expected-domain cause when available. Initialization failures
before any residual attempt have zero solve counters/metrics. The internal callback
counts completed Arnoldi iterations; users cannot install a callback or
preconditioner. Only owned `_ImplicitHybridTrialError` is caught for shrink or
translated into a typed convergence failure. Other canonical exceptions and
cancellation propagate. Counter overflow also propagates before publication.

## Physical phase, publication and diagnostics

Accepted positions, velocities, B, E and moments are all integer-time values.
There is no priming or half-velocity conversion. Count zero has stored-zero time
and timestep history and no successful report. An advanced state has positive
finite time and positive resolved `last_dt`; time advances by actual stored
additions, never `count*h`. Positive `dt` requires a resolved half step and strict
stored clock progress. An originally nonzero timestep that converts to stored
zero is rejected, including a tiny negative value that becomes -0. An originally zero timestep validates the
ready state and returns without staging or changing the report.

All calculations use private state. Same-type particle replacement leaves,
endpoint fields, counts, clock and report are prepared before the first public
store. Publication preserves the public `ParticleSet` objects and field arrays;
particle leaf arrays are replaced. No resize, validation, deposition, solve or
report calculation follows the first store. Recoverable failure leaves public
state and the previous report unchanged; private work may change. Fatal resource
failure, protected memory and concurrent mutation are outside this boundary.

`implicit_hybrid_diagnostics(st, species)` returns the synchronous energy and
momentum budgets, clock, zero velocity-time offset and last accepted report.
The report retains both scaled residual blocks, scales, costs, canonical endpoint
electron free energies, pressure work and a density-continuity defect. With the
actual canonical `ue = ui - J/n`, represented midpoint `gradp` and raw `Mbar`,

```
pressure_work = RN_T(-h * sum(ue · gradp_mid) * cell_volume)
continuity_defect = RN_T(maxabs(Q(n1)-Q(n0)+Q(h)*Q(div(Mbar))))
free_energy_defect = RN_T(Q(F1_T)-Q(F0_T)+pressure_work_exact)
```

Work and the continuity combination accumulate stored factors exactly before
final rounding. The standalone canonical `div(Mbar)` must complete first; this
family does not fuse its range with a later h factor. `F0_T` and `F1_T` come from
the canonical isothermal energy owner. NaN endpoints reject; if either endpoint
is infinite, the report uses `:unresolved_free_energy_range` and
`free_energy_defect = nothing`. Finite endpoints use `:represented_reservoir`;
the resulting work/defect may themselves exceed T's range. The separate canonical
[complete-budget source](HYBRID_DIAGNOSTIC_RANGE_PLAN.md#current-complete-budget-source-2026-11-02)
retains original reservoir expressions for range-sensitive totals; it does not
change this report's represented F0_T/F1_T boundary or its unresolved policy.

This is a represented-reservoir diagnostic, not an exact electron-work theorem.
The isothermal free-energy secant leaves continuity and nonlinear chain-rule/work
defects, in addition to the actual normalization, quadrature, FFT, moment and
Boris rounding. The [Stanier–Chacón–Chen reference](https://arxiv.org/abs/1803.07158)
advances electron pressure and uses its own compatible discretization and orbit
averages. Its energy/momentum theorem is not transferred to this spectral/CIC,
single-midpoint represented-map method.

## Checkpoint and archive

`save_checkpoint(path, st, species)`, `load_checkpoint!(st, species, path)` and
`archive_run(path, st, species; rng_seed, ...)` use distinct schema 1
`:implicit_hybrid_midpoint`. The existing atomic writer and schema-4 archive
checksum remain the file owners. The seed is provenance; this family has no RNG.

The exact payload keys are `integrator`, `schema`, `configuration`, `phase`,
`initialized`, `species`, `fields`, `floor_count`, `time`, `step`, `last_dt` and
`last_report`. Configuration includes D/T/grid n/L/dx, the exact fixed normalized
model/CIC/charge-mass binding, all resolved controls, and residual/quadrature/
phase/moment/derivative/density/field/norm/solver policy tags. Unknown keys,
other families/versions and mismatched fixed policy are rejected.

Only initialized completed states can be saved. Ordered x/v/weight/ID/tag/q/m
and every public `_hybrid_field_state` array are saved verbatim: n/ui/B/E/J/pe/
gradp/ninv/lapJ and the correctly empty scalar CGL roles. Public inactive lapJ
values are exposed state and survive restoration; residual Ohm inputs are
separate private zero arrays. Floor count must be zero. The typed last report
must match clock/history and fixed-control metadata bounds. Saved reports are
history metadata, not proof of the unavailable preceding solve.

Restore validates destination storage, owner identity and cached immutable
binding while allowing stale destination physical values, particle q/m, counts,
clock, initialization and report. Incoming ordinary arrays must be finite,
independent and detached from destination public/private/grid storage; incoming
strict density and complete phase/report policy validate. All incoming state is
detached again, complete same-type replacement leaves and counts are prepared,
and solver capacity is reserved for the incoming local counts before any public
store. The species count/order and cached q/m must match; local particle counts
may change. No initialization, redeposition, Ohm recomputation, solve, projection
or leapfrog migration is performed. No midpoint/initial guess/Jv/Krylov/FFT work,
plans or other scratch is persisted.

This implementation has received source review only. No Julia, compiler, parser,
numerical, convergence, temporal-order, dispersion, work-error, rollback/restart,
resource or performance acceptance was executed in this implementation phase.
