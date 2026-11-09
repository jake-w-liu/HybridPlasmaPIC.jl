# Coupled EMPIC projected-pair contract V3

Design only against HybridPlasmaPIC123c4a23e59923f02ce96150c34ce32ccae72442.
This amendment completes the implementable first controlled stage of inherited
CONTRACT_V2.md; the original and its source freeze remain untouched. No production
approval, execution, order/effectivity theorem or full CAP008 closure is implied.
SOURCE_PINS.json pins32 repository files and the actual SpectralOperators v0.1.9
source required by Project.toml. Existing Nothing/dividing/persistent half phases
and schema1/2/3 keep their meanings. The first new opt-in is deliberately bounded
to CIC/TSC finite-interval evolution across support knots; NGP jump/error
control stays OPEN and is rejected at constructor selection for this new policy.

## 1. Actual spectral model and the limiting constrained equations

VERIFIED owners: full_pic.jl:2412–2472, electrostatic.jl:123–244; dependency
periodic_fourier.jl:75–115. Use the ACTUAL stored effective vector
k=(g.kvec[d][I[d]])_d, not continuum2pi*m/L or kfull. Nyquist is zeroed per even
axis. Let R contain modes with at least one nonzero effective component, Z the
non-DC all-zero-k combinations, and0 the true first Cartesian/DC index. Embed
D spatial components in3; absent gradient components are zero.

On R define D_k e=i*k·e, G_k r=-i*k*r/|k|², L_k=G_k D_k, P_k=I-L_k.
Then D_k G_k=1. On DC set G=0,P=I; on Z set G=P=0. Curl C uses this same
effective k (i*k×), absent derivatives zero. Thus DC raw current is physical,
ALL components of Z current are zero, and mixed Nyquist/nonzero-k modes keep
ordinary transverse components. Curl is zero on DC/Z. Initial Z E/B values are
constant ideal homogeneous modes, not projected away by current correction.
Gauss applies to R: D E=R rho. Physical DC charge and excluded Nyquist density
remain observed source channels; do not claim a full-grid physical Gauss law.

Let rho(X) be the ideal real-polynomial counterpart of the SAME tensor stencil,
fixed physical source/neutralizer and uniform charge-reconciliation owner.
Reconciliation retains the prescribed original source target; its ideal uniform
shift is X-independent by partition of unity. Actual stored rho includes stencil,
product, sum and final reconciliation rounding. Jraw(X,V) is the canonical
original-q*w*physical-v deposit at full-time stored positions of ALL mobile rows;
fixed rows contribute charge only. Define dotrho=D_X rho(X)*V_D on a smooth piece.
The finite-interval correction is, on R,

    C_h(J,Delta rho) = P J - G(Delta rho)/h,
    D C_h = -Delta rho/h.

This sign follows directly from canonical source=-delta/divisor in
_spectral_project_mode. On DC C_h=Jraw; on Z C_h=0. Consequently the limiting
constrained ODE is

    X'=V_D,  V'=F(X,V,E,B),  B'=-C E,
    J*=P Jraw-G dotrho,
    E'=c² C B-P Jraw+G dotrho.

The longitudinal current is NOT the longitudinal raw velocity deposit. Reduced
primaries are X,V,E_free,B plus fixed initial longitudinal residual, with
E_L=G rho(X)+constant. DC E'=-Jraw_DC; Z E'=0. This is the precise target on
smooth pieces. For NGP, rho jumps at support cuts: dotrho then has a measure
part and no ordinary differentiable ODE is asserted at those jumps.

## 2. One fully coupled projected Euler–Heun trial

Use one exact calendar gap delta>0, actual resolved h::T from the existing
_empic_event_interval(T,delta), and H=Q(h). BOTH branches consume H, never delta
in one branch and h in the other. Origin E0/B0/V0 are full-time physical inputs.
From origin gather ALL mobile forces F0 and deposit J0 at stored positions.
The equations below are exact-arithmetic definitions; all implemented field rows
consume canonical derivative parts and current correction, not copied operators.

    XE=X0+H*V0_D, VE=V0+H*F0, BE=B0-H*C E0,
    rhoE=rho(XE), JE_used=C_h(J0,rhoE-rho0),
    EE=E0+H*c²*C B0-H*JE_used.

Wrap/materialize XE with canonical represented-lift ownership. Materialize VE,
EE,BE; require finite values, supported relativistic C1 and branch/lift admission.
These ACTUAL represented Euler values supply stage1 gathers F1 and current J1.
Then form, from the SAME origin (not from a mutated Euler state),

    XH=X0+H*(V0_D+VE_D)/2, VH=V0+H*(F0+F1)/2,
    BH=B0-H*C(E0+EE)/2,
    Jbar_raw=once_round_T((Q(J0)+Q(J1))/2),
    rhoH=rho(XH), JH_used=C_h(Jbar_raw,rhoH-rho0),
    EH=E0+H*c²*C(B0+BE)/2-H*JH_used.

There is no B half prime/center, velocity half prime, Boris finite-map substitute,
second old macro call or extrapolated primary. Generic branches retain all3 raw
and consumed currents. Planar stages retain only actual raw Jy plus physical
Jx_DC; canonical _esirkepov_Jx! supplies consumed Jx from EACH branch's charge
endpoints and its actual DC (stage0 for Euler, once-rounded stage average for
Heun). Nyquist suppression of Jy remains canonical. Do not compute/report an
unconsumed full raw Jx, or solve/re-prime Ex again after consuming this Jx.

Important V2 correction: Euler and Heun normally have DIFFERENT represented
positions and deposited rho. They agree in original source inventory, source
calendar endpoint, actual h and physical full-time phases, not spatial charge.
The pair is therefore compared at a common physical phase, with its own complete
endpoint source on each branch. No old pending half-phase subtraction is used.

## 3. Projection consistency and the actual defect

For this algebra only remove ALL materialization/FFT/current-average rounding.
In ideal exact arithmetic on a smooth admitted piece, the equations in section2
are Euler/Heun on the reduced X,V,E_free,B system followed by the constraint lift
E_L=G rho(X). In particular,

    XH-XE = H*(VE-V0)_D/2,
    VH-VE = H*(F1-F0)/2,
    BH-BE = -H*C(EE-E0)/2,
    EH-EE = H*c²*C(BE-B0)/2 - H*P(J1-J0)/2 + G(rhoH-rhoE).

The last term is load-bearing. A raw-current difference alone misses it.
Expanding rho(X0+H V0) includes H²*D²rho[V0,V0]/2 already in BOTH lifted branches;
it cancels from their difference. The leading longitudinal pair term is instead
H²*G*D rho*F0_D/2. Conditional Taylor algebra on a sufficiently smooth fixed
piece therefore has the expected reduced embedded consistency; it is not an
ordinary full-state Euler defect and supplies no true-error or global-order
statement for the implemented projected/rounded/cut-crossing method.

Actual observations are computed directly from Q(high)-Q(low), with no Richardson
divisor or assumed convergence exponent. Preserve these perturbations explicitly:
represented lift minus ideal drift; stored-source/stencil/reconciliation error;
native FFT/curl/projection error; current-average rounding; complete field-row
rounding; and supported-output materialization. They are not silently set to0.
Actual source/current/stage-curl operands are retained. Exact-Q row receipts use

    etaE_b=Q(E_b)-Q(E0)-H*c²*recorded_curlB_b+H*Q(Jused_b),
    etaB_b=Q(B_b)-Q(B0)+H*recorded_curlE_b.

Here recorded derivatives are exact sums of the actual T parts and Q scales
returned by canonical range owners. These eta identities are checkable without
calling an ideal DFT. Analytically Jused_b=C_h^ideal(Jraw_b,Delta rho_b)+xi_b;
xi_b includes native spectral perturbations. Then the ideal difference identity
has additional -H*(xi_H-xi_E), eta and derivative-part errors. The first stage
does NOT claim to compute exact ideal-DFT error xi or exact ideal-stencil error.
It reports actual continuity/Gauss observations and the real pair difference;
independent essential references must quantify numerical floors remotely.

A native field-space constraint observation projects a DETACHED copy with the
canonical divergence projection and compares it to the actual endpoint, in E
units; analogously observe B's range longitudinal component in B units. Require
finite observations within the corresponding explicit field tolerance. These
copies are observers, never repaired/published fields. DC/Z source channels and
native continuity residual are reported separately; zero charge increment must
not create a fictitious small relative continuity claim.

## 4. Canonical scaled force and strict physical-v admission

Propose two scalar methods in particles/kinematics_range.jl; no old DAG changes:

    _em_lorentz_scaled_parts(v,E,B,factor,H::Q,c,Val(rel))
        -> (coefficient::NTuple{3,Q},gap::Q)
    _em_lorentz_pair_output(v0,p0,p1,c,Val(rel)) -> NTuple{3,T}

p1=nothing is Euler; otherwise compose the2 increments with exact1/2. NR parts
reuse _euler_kick_parts on original Q operands and original _boris_exact_factor,
returning its scaled kick parts4:6. Relativistic physical-v acceleration is

    F = (q/m)/gamma * [E+v×B-v*(v·E)/c²],
    g=1-|v|²/c², coefficient=H*(q/m)*[...], increment=coefficient*sqrt(g).

Derivation: p=m*gamma*v, p'=q(E+v×B), v·p'=q*v·E gives the displayed inverse
mass-tensor formula. Existing relativistic Boris/priming are finite momentum
maps and cannot replace this RHS. Preserve original q/m even if rounded qm is
Inf/0, and never materialize an unscaled acceleration, T(H*q/m), T(h/2) or T
stage increment before composition. Each full output is rounded only once.

New source trace improves the independent force review's fixed256 proposal:
reuse existing _combine_relativistic_radical! followed by
_finish_relativistic_energy at extrapolated_leapfrog.jl:2001–2080 for
Q(v0)+a0*sqrt(g0)+a1*sqrt(g1). FIRST canonicalize both weighted terms with the
combine owner: it extracts rational roots and merges square classes. This is
load-bearing; the finish owner assumes those moment-caller preconditions and
its single-term canonical path is not a standalone rational-root tie resolver.
It already rounds canonical radical sums for momentum as well as energy, using
exact isqrt enclosures/cancellation. At most2 terms occur here. No new radical
implementation or fixed256 rounding assumption is needed. If the canonicalized
net increment is exactly zero, return the original stored component bits.

Both stage inputs and EACH materialized Euler/high velocity require canonical
_supported_relativistic_output in addition to positive exact gap. Rounded beta²>=1
still rejects even when the original exact vector is subluminal. No clipping,
normalization, speed-domain widening or invalid predictor excused by a valid
high endpoint. Empty/fixed rows retain canonical semantics; q0/weight0 MOBILE
markers participate in all force/domain/error checks. Planar omitted Ez/By and
stored vz remain canonical zeros. New force helpers are shared particle owners,
not a wrapper's copied pusher. Their remote exact-radical/C1 tests are required.

Canonical field caller proposal in full_pic.jl:
_empic_full_time_field_update!(out,origin_fields,derivative_stages,weights,h).
weights are exact (1,) for Euler or(1/2,1/2) for Heun. Each stage packet owns
its actual curl/derivative arrays and optional normalized Q-scaled parts once;
recipes reference packets, not overlapping duplicate array roles. Combine those
parts with the exact weights and reuse _em_magnetic_update! with preceding=h,
following=h (one full interval) and the existing complete electric scalar rows.
Planar Ex consumes actual Jx, Ey consumes minus Dx Bz and Jy, and Bz consumes
Dx Ey. Never materialize a rounded mean derivative or field before this row.
The additive caller does not edit any old centering/priming/DAG or field equation.


## 5. Usable CIC/TSC crossings and explicit NGP selection boundary

Use the actual rounded s_T=mod(x_T,L_T)/dx_T stencil owner. NGP cuts are k+1/2
with ties UP and discontinuous values; CIC cuts are integers with continuous
values/derivative jumps; TSC cuts are half-integers with ties-even, continuous
value/first derivative and second-derivative jumps. Q(x)/Q(dx) is a geometric
guard, not a replacement deposit/gather coordinate. Physical seams are additional
cuts when Q(L)!=n*Q(dx); a continuum torus smoothness argument cannot hide them.

The first closed policy is domain=:cic_tsc_continuous_v1. The constructor rejects
NGP with nonnothing error_control BEFORE allocating an authoritative controlled
owner; old NGP Nothing/dividing/persistent schedules remain available unchanged.
This is an explicit new-policy restriction, not a hidden fallback or mid-run
stall. Full NGP jump/force-side/continuity/work/restart and effectivity acceptance
remain required follow-through. No generic NGP error-control claim is made.

CIC/TSC trials DO cross their support knots. Evaluate every origin/Euler/high
source and force with the actual canonical stencil, including integer CIC and
half-integer TSC ties; no fake smoothing, one-sided override or zero-duration
event. Charge endpoints pass through the same finite-interval correction and
Maxwell work receipts across a knot. The source/stencil and FFT perturbations
remain in the observed branch difference. Smooth Taylor reasoning in section3
applies only when its smooth-domain hypotheses hold; nonsmooth/crossing trials
remain measured embedded observations with no inferred order/effectivity theorem.

For every mobile marker/axis retain actual wrapped origin/Euler/high values,
BigInt winding, represented lifted branch endpoints, ideal endpoints and exact
rounding defects, and the canonical base/weights at all3 sampled states. Require
both ideal AND represented branch displacements <L/4 and high-low lifted
difference <L/2, so the signed torus metric has a unique unaliased branch. This
geometric travel gate is independent of support-knot crossing and is shrinkable
from ordinary loaded knot positions; it does not require a single cell/piece.
Use compact checked integer cut ranges/directions (not one array per knot)
for the straight origin-to-branch represented segments, plus explicit seam/cut
labels. These are discrete stage-path receipts, not an exact-orbit event history.
Actual stage stencils supply the consumed fields and source. When Q(L)!=n*Q(dx),
record the canonical seam as a represented geometry cut; do not replace dx by
L/n or claim smooth periodic matching there.

A support-knot crossing is NEVER alone a rejection cause. Finite/domain/C1,
travel, norm, field/source, clock or work-budget gates may reject transactionally
and halve from the unchanged origin. In particular an exact CIC/TSC knot origin
with finite motion is admitted under the same continuous-stencil policy. The
NGP negative-drift/ties-up counterexample in FORCE_DOMAIN_REVIEW.md explains why
that historical open-piece proposal was rejected; it is not this runtime policy.

## 6. Closed controller, clocks, calendar and failure

Constructor-only owner selection: error_control=p requires species_dt/event_budget
and forbids species_stride. p has EXACT keys

    (pair=:projected_euler_heun_v1,metric=:component_max_v1,
     domain=:cic_tsc_continuous_v1,atol=(x,v,E,B,clock),rtol=(x,v,E,B,clock),
     h_min,h_max,max_trials,max_force_evals,max_deposits,max_fft_calls,max_report_cells).

Every new-owner operation requires nearest rounding and gradual underflow.
Role atols are finite positive T after checked original conversion; role rtols
finite nonnegative T<1. No shared tolerance with hidden units. h_min/h_max are
resolved positive T with h_min<=h_max; all budgets positive fitting Int. One
trial tests all mobile3V, D positions and every represented E/B component. Scale
v/E/B by atol+rtol*max(abs(start),abs(low),abs(high)); scale position error by
atol_x+rtol_x*L_d, using the UNIQUE actual lifted high-low difference guaranteed
by the actual lift/travel admission. No cut-dependent absolute-coordinate denominator.
Fixed/provenance/source leaves match exactly and do not dilute the maximum.

At source tau choose delta=min(Q(next_h),next mandatory edge-tau), h=interval(T,delta).
Mandatory edges are the old exact base partition and pending nominal H/a/b
calendar, never merged/clipped/retimed. Nonresolved or below-hmin edges reject.
Use independent clock defects Q(h)-delta and Q(t_next)-Q(t)-Q(h), plus their
sum; clock norm is their maximum absolute value divided by the declared clock
scale (atol_clock+rtol_clock*max(abs(t),abs(t_next),abs(tau),abs(tau+delta))).
Pair phases/source edges agree exactly. Source tau advances by delta; public
clock uses _validated_step_clock with actual h, not source time replacement.
Typed trial rejection is limited to explicitly observed norm, finite output, C1,
travel and constraint gates. Do not catch all ArgumentError or match error text.
Canonical source/FFT/programming failures propagate as external-operation
failures with unchanged public state unless a separately approved typed
materialization-status seam makes that precise failure retryable.

Accept only norm<=1, field/source/branch/C1/role/phase/work gates all valid.
Reject BOTH branches otherwise. Retry delta/2 from identical accepted origin;
never consume a nominal event or advance source/public state on rejection.
After success recommend min(h_max,2h) only for norm<=1/4, otherwise h, clamped to
[h_min,h_max] by exact cap conversion. Halving/doubling and threshold1/4 are
explicit conservative algorithm policy, not an asserted error exponent. Stop
with a specific hmin/resolution/domain/budget failure. No guarantee every
requested tolerance or macro endpoint is reachable is made.

Nominal due-row decisions use actual accepted full-time samples and canonical
cap/rotation/electric/travel limits, but cause NO velocity kick/prime.
Extract the existing interval-choice calculation into one canonical helper; the
old persistent decision wraps it with its unchanged half-duration receipt, while
the new owner consumes only the interval choice. Do not duplicate limit policy.
Candidate due decisions are replayed/discarded with a rejected trial; no source
event becomes authoritative until an actual interval commits privately. New labels
are nominal_H/last_deadline/next_deadline, not old last_kick/half timestamps.
Right-edge decisions remain pending until the next positive call; zero dt
validates and consumes nothing. Every mobile row force-steps at BOTH stages,
even while its nominal interval is pending. No multirate force saving is claimed.
Reuse _empic_event_next with the new phase's time/interval/last/next calendar
interface and the same captured macro origin; do not duplicate the selector.
The existing event_budget bounds accepted intervals plus deadline decisions in
this new policy; max_trials additionally bounds ALL attempted pair trials.
Public macro step increments ONCE after all private accepted intervals complete.
Physical velocity/B offsets are0; last_sub_dt remains last actual h metadata.

Future control/cap updates use existing family schedule setters only on this new
concrete owner: admit/copy a complete proposal before stores; pair/metric/domain
are immutable, pending nominal H/a/b and old reports remain unchanged. Adjust
future next_h to the new bounds without rewriting its preceding acceptance.
Each old receipt carries its OWN historical proposal; validators never judge it
against a later tolerance/budget. Old owners reject nonnothing error_control;
no setter or restore silently enables/migrates a phase.

## 7. Work, budgets, transaction and complete persistence

Per accepted interval retain actual origin/Euler/high particle stage operands,
E/B, rho0/rhoE/rhoH, stage raw currents/DC, raw effective and consumed branch
currents, actual native derivative parts/scales, lift/rounding receipts, norm
components, source/public clocks, deadline decisions and historical proposal.
Retain only bounded scalar failure summaries for rejected trials, not their
physical work/impulses. Preflight checked report-cell growth BEFORE retaining
arrays. Last external report is bounded; cumulative ledgers are scalar receipts.

For published high branch, Wused=H*dV*sum(Jused*(E0+EH)/2). Generic Wraw uses its
actual complete raw effective current with the same E endpoints. Planar raw work
is Jy-only and labeled transverse; never invent full raw-longitudinal work.
Particle electric quadrature Wp=H*sum(q*w*(v0·Ep0+vE·EpE))/2 uses actual gathered
samples, original factors and physical v. Report signed

    Dfield=DeltaU+Wused, Dparticle=DeltaK_reported-Wp,
    Dcoupling=Wp-Wused, Dtotal=DeltaU+DeltaK_reported.

Their exact-Q sum/ledger is checked, not forced zero. U uses canonical exact
represented-field energy. NR K is canonical rational energy; rel K_reported is
Q of the finite canonical physical relativistic moment output, explicitly a
reported rounded measurement, not an exact irrational-energy claim. A required
nonfinite measurement rejects visibly. Persist initial reported energy and all
cumulative defects; Energy_reported-Energy_initial=Dtotal_cumulative. No passivity/conservation or
stability theorem is inferred. Rejected attempted costs never enter this ledger.

Count all attempted particle force evaluations, scalar particle-deposit
contributions (including DC), retries and retained report scalar cells. Fourier
budget is conservative reserved primitive-plan calls: pinned v0.1.9 derivative2,
curl4 inD1/6 inD2/3; generic continuity correction7; planar Jx2; divergence D+1.
Thread a CALL-LOCAL optional meter through canonical normalized derivative and
continuity replay owners, charging every actual replay transform before dispatch.
Do not charge accepted intervals only or hide band replays. The reservation is
an algorithm work bound for pinned owners, not measured FFT/CPU time; reverify
counts on dependency changes. Scratch/peak-byte costs remain measured acceptance.

Extend canonical collection staging with independent low/high candidates and a
committed full-time shadow; inspect ALL original/public/private/report/FFT roles.
Snapshot the admitted canonical n/L/dx/midx/kfull/ik/kvec operator bits as fixed
configuration; compare against the actual grid and save that signature. Validate
coherence through the dependency's canonical _fourier_metadata (no FFT plans),
never a copied wavevector formula. Restore rejects a changed symbol signature
rather than reinterpreting saved physical fields on it.
A new owner-local AtomicBool guards init/step/setter/diagnostics/capture/restore;
old owners retain their concurrency semantics. All branches, final cumulative
ledger, replayable report/controller/calendar and detached shadow validate before
ordinary stores. Any late failure discards the entire external operation,
including its private accepted subintervals. Guards release in finally. Shared FFT scratch may change. Ordinary writable/injective
CPU storage and independent-state/grid concurrency limits remain unchanged.

New same-family schema4/policy3 tag phase_policy=:full_time_coupled stores physical
primaries/source/provenance, generic raw or planar actual field inventory, exact operator signature, complete
full-time phase, nominal calendar, future proposal/next_h, last historical proposal/
acceptance and deterministic counters, last complete bounded report, cumulative
energy/work/clock observations. Count-zero initializes from full-time data only;
NO schema1/2/3 advanced or zero-state conversion is silently supplied. Save only
completed boundaries. Original graph and live links are checked BEFORE copying;
strict finite fields/own-source/phase/report/accounting/schema validations precede
verbatim restore and prospective-shadow publication. Stale old numerical values
and changed ordinary destination counts follow canonical restore admission, not
an implicit reinit/reprojection. Archive container remains canonical.

## 8. Exact proposed ownership and essential remote gates

NEW exclusive implementation paths after approval: src/models/empic_coupled_error_control.jl
(concrete controller/pair/calendar/receipt owner), src/diagnostics/empic_coupled_error_control.jl,
src/io/empic_coupled_error_control_checkpoint.jl, docs/EMPIC_COUPLED_ERROR_CONTROL.md,
test/test_empic_coupled_error_control.jl. No copied Maxwell/current/source solver.

Existing seams requiring explicit approval: particles/kinematics_range.jl adds the
2 scalar scaled-force/output methods; models/empic_collection.jl adds constructor
keyword/owner selection, init/advance/diagnostic/setter/schedule-array hooks and
canonical interval-choice extraction with literal legacy arithmetic preserved;
models/full_pic.jl adds ONLY optional local Fourier-meter propagation through
canonical derivative/current/continuity helpers and a weighted full-time field
caller using their original complete scalar update owners. IO collection dispatcher
adds strict schema4 routing and rejects reserved fields in old schemas. Root owns
module includes/exports/shared writer/global docs/workflow/full registrations.
Initialization uses canonical charge/Gauss preparation, zero-time calendar and
zero-time work/controller counters, with NO prime; initial range-divergence B
and E constraint observations must meet the explicit B/E role tolerances.
Reuse canonical _finish_relativistic_energy; do not edit/duplicate its radical
math. Reuse the root-adopted shared represented-lift API after periodic integration;
it is not yet present in the pinned123c4a2 and must be re-pinned before edits.

Compact independent REMOTE gates: actual effective-k transverse/longitudinal/DC/
pure/mixed-Nyquist current and endpoint Gauss reference; interior CIC/TSC projected
Euler/high formulas and the G(rhoH-rhoE) term; varying vacuum wave plus charged
weight0 marker during a nominal pending interval; nonzero-weight fully coupled
refinement/effectivity at common full-time endpoints; original q/m range and
exact radical parallel/perpendicular/near-C1 outputs; actual exact/neighboring support
cuts and seam crossings with consumed branch-stencil receipts; explicit NGP
constructor restriction; late hmin/C1/travel/budget/constraint failure atomic
retry; zero/fixed/empty/non-dyadic/right-edge clocks;
changed-count/stale-destination raw/archive continuation after retries and future
proposal updates; malformed historical-vs-future config/source/alias/schema
rejection; unchanged old Nothing/dividing/persistent and half-phase suites.

No local Julia/parser/compiler/formatter/numerical/proxy/test/device/profile work,
CI, commit, push or process control. Algebra and actual source map support this
contract; numerical/scientific/effectivity/resource acceptance remain unexecuted.
Full NGP jump/error-guarantee handling and CIC/TSC crossing effectivity,
asynchronous error-controlled
multirate force reduction and source/collision/boundary/backend intersections
remain required follow-through, not completed by this bounded opt-in stage.
