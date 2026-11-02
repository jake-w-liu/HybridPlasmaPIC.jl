# CAP-008 specialized EMPIC1D species implementation

Updated 2026-11-01. **S1 and S2 sources are implemented and reviewed;
exact-revision GitHub acceptance remains pending.** The [audit tracker](DEEP_DEBUG_TRACKER.md) owns
publication and CI status. This plan expands the specialized-model and multirate
items in [CAP-008](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-008--arbitrary-kinetic-species-in-esem).
The [generic EMPIC collection contract](EMPIC_SPECIES_IMPLEMENTATION.md) remains
the source of its existing scope and pending acceptance.

The concrete S2 schedule/phase/persistence implementation is integrated after
root and independent source review. The API below is source-present; numerical,
compiler, formatter and resource acceptance require its exact-revision remote
checks. The original synchronous S1 schema is retained below as migration history.

## Scope and source owners

`src/models/full_pic.jl` carries a specialized `EMPIC1D` with `Ex/Ey/Bz`,
longitudinal continuity current, midpoint transverse current, staggered magnetic
and velocity phases, and legacy electron/optional coarse-ion stepping. Its
collection overload now uses those specialized owners. Generic `EMPIC{1}` already carries both
transverse polarizations; it is not an adapter for this implementation.

S1 adds actual specialized CPU Float32/64 collections with arbitrary finite
signed charge, positive finite mass, mobile/fixed populations, explicit
background/neutrality, NGP/CIC/TSC and optional relativistic motion. The reduced
polarization requires the invariant planar subspace: every original stored
`v_z` is exactly zero, including rejecting nonzero stored subnormals. Fixed
populations require zero velocity. The legacy scalar API retains its accepted
inputs and numerical ordering. Shared field energy is counted once; particle
momentum is reported as particle momentum.

Reuse the canonical collection storage, charge, owned-stage, binding and phase
owners in `electrostatic.jl`/`empic_collection.jl`, original-factor Boris/prime
owners, specialized current/Maxwell blocks, and collection checkpoint/archive
owners. Extract shared legacy blocks without changing their expressions or
order. Reject known unsupported ordinary-store ancestry before staging can
hide public storage failures. Stage all specialized transaction arrays,
particles and clocks for the complete external operation before publication.
Private scratch may change on rejection; public destination identities remain.

## Implementation and acceptance tracker

| Stage | Status | Required result |
| --- | --- | --- |
| S1 source design | Reviewed | Genuine specialized state; planar scope; synchronous fast phases; complete publication and distinct restart contract |
| S1 construction and orchestration | Source implemented | Ordered populations gather the same old fields, then sum currents before one specialized field update per fast interval |
| S1 diagnostics and persistence | Source implemented | Physical moments and verbatim field restore, complete graph detachment, saved-mode validation, identity/capacity handling and layout rejection |
| S1 source review and essential checks | Source reviewed; runtime pending | Legacy helper-expansion comparison and independent caller review; essential current/phase and failure/retry/restore controls registered; existing assertions retained |
| S1 GitHub acceptance | Pending | Exact revision on minimum `1.12` and moving stable `1`; full neighboring tests, complete compiler review and formatter |
| S1 resource acceptance | Pending | Remote complete init/step/restart cost versus specialized legacy path, actual hardware/runtime/threads, allocation traffic and peak memory distinguished |
| S2 schedule/phase/schema design | Reviewed | Dividing strides, independent previous intervals, atomic future-configuration updates, nominal offsets and explicit S1 migration |
| S2 construction and advance | Source implemented and reviewed | Due-event kicks with every-fast drift/current, whole-operation history staging and atomic future schedule setter; generic/scalar/all-fast source reductions retained |
| S2 persistence and essential checks | Source implemented and reviewed; runtime pending | Schema 2/policy 1 and tagged S1 migration; independent interior-event, legacy, pair/unequal-mass refinement, rejection/retry and restart controls registered |
| S2 GitHub/resource acceptance | Pending | Exact revision on minimum/latest Julia; full neighboring tests, complete unsuppressed compiler review, formatter and remote complete-operation cost/peak-memory evidence |
| Broader CAP-008 integration | Open | Collisions, boundaries, partition ownership, device residency and matched physical/cost validation |

All Julia tests, compiler checks, formatter, builds and profiling run on GitHub
or suitable remote hardware. Local work is source/caller review and lightweight
file checks. Add compact essential regressions through the existing test
framework and `electromagnetic-species` selector; do not add a test campaign or
replace neighboring checks with mocks.

## S1 synchronous baseline and migration history

Every mobile population kicks and drifts on every field interval `h=dt/n_sub`.
Total old/new charge includes fixed populations and the explicit background.
The longitudinal DC current uses the complete original mobile `q*w*v_x/L`
sum; represented nonzero modes use the existing continuity solve. Sum midpoint
`Jy` before the existing Nyquist treatment, including whole-source range replay
where the canonical current owner requires it. Do not reject an overflowing
partial species sum when its complete physical sum can cancel.

First advance primes each mobile velocity and B once. A changed `dt` or `n_sub`
retains the previous fast interval: the kick spans previous-half plus new-half,
while drift/current use the new interval. Store full and fast histories
separately. A zero step validates inputs without consuming the prime or phase.
S1 `n_sub>1` is synchronous; legacy ions remain coarse movers. Compare legacy
and collection APIs on their common planar `n_sub=1` inputs.

S1 used a distinct schema-1 `:empic1d_collection` payload with required
`polarization=:ex_ey_bz`, ordered particles, policy/mobility, initialization,
physical `Ex/Ey/Bz/Jx/Jy/rho_n/rho_np1` and full/fast clocks. Do not persist
gathers, centered B or other scratch. Protect actual scalar init/step/save
entry points from collection states before mutation. Validate saved particles
against saved c/mode and planar constraints, detach all incoming leaves, reserve
ordinary capacity, and restore fields verbatim. Legacy and generic collection
layouts must reject cross-loading before publication.

Essential checks cover independent current/phase references, variable-step
continuation, fixed/neutral/typed-empty populations, stored nonplanar rejection,
one nontrivial relativistic case, late failure/retry, and checkpoint/archive
ownership and layout boundaries. Source comparisons are not numerical proof;
all runtime results remain pending until reviewed at the integrated revision.

## S2 dividing per-species multirate schedule

For external fast count `M=n_sub`, give each mobile species an integer stride
`r_s` dividing M. Kick at fast indices divisible by its stride, but drift every
mobile species and deposit total charge and midpoint current on every fast
interval. Update fields once after all species contributions. `r_s=1` recovers
S1; electron stride 1 and ion stride M provide the legacy comparison target.

The constructor keyword is `species_stride=nothing`: default mobile
rows use 1 and fixed rows use 0. Explicit rows must be one-based integers fitting
Int, with Bool rejected. Mobile strides divide M; fixed rows are exactly 0.
Typed-empty and zero-particle mobile populations retain their declared schedules.
The `set_empic1d_schedule!` setter validates and copies M/strides together
before publication. It works at completed boundaries, including step zero, and
preserves prior phase histories. Direct `n_sub` edits are checked against current
strides before init/step/diagnostics/save mutation; zero steps still validate.

Define `b_s=M÷r_s`, then the represented coarse interval `H_s=dt/b_s`.
Use zero-based events `j%r_s==0`, or `(k-1)%r_s==0` in a one-based fast loop.
Prime each mobile velocity by its own H_s/2 and B by fast h/2. The first ever
kick uses virtual previous H equal to new H; later external first kicks consume
saved previous H/2 plus new H/2. Interior events then use new halves. Gather and
kick only due species, but drift and deposit current from every mobile species
at every h. Preserve one common field update and the existing numerical owners.

Each specialized collection owns strides, authoritative `last_species_dt`
and private history candidates. Generic EMPIC keeps synchronous behavior by
default and now shares this owner through an explicit dividing-schedule opt-in;
see the [generic contract](EMPIC_SPECIES_IMPLEMENTATION.md).
Stage histories with the entire external operation and publish only after all
candidate checks. At step zero histories are stored zero; afterward each mobile
history is finite/resolved and lies between saved last_sub_dt and last_dt, with
fixed history zero. Do not reconstruct previous durations from current settings:
the schedule may have changed after the last advance. Currents/rho still describe
the last fast interval, so continuity diagnostics use last_sub_dt.

Specialized S2 writes use schema 2 with required Int `schedule_policy=1`,
strides and per-species histories; archive container 4 stays. Only validated
specialized schema-1 `:empic1d_collection`/`:ex_ey_bz` state migrates to all-fast:
advanced mobile history comes from saved last_sub_dt, fixed/step-zero history
is zero. Saved current n_sub may differ from the configuration used for its last
advance. Reject reserved S2 fields in schema 1 and unknown versions/policies.
Saved M/strides replace destination configuration together. Normalize once,
validate/detach the complete incoming graph, reserve capacity, then restore
verbatim fields and phase. Destination requested settings being replaced do not
govern incoming schedule validation. Generic/scalar layouts remain distinct.

Velocity offsets `-last_species_dt/2` and magnetic offset `-last_sub_dt/2` are
nominal declared stage durations. H_s need not equal r_s*h exactly, and repeated
rounded h additions need not equal one full-dt addition. Preserve that existing
arithmetic; do not claim exact timestamps relative to the rounded public clock.

Gate all-fast identity, actual legacy coarse-ion comparison, variable intervals,
history rollback/retry, explicit migration and restart after changing a schedule
before another advance. A compact independent interior-event oracle and positive
pair/unequal-mass refinement signals must expose held currents and phase mistakes;
observed orders apply only to those problems. Continuity for the split path does
not establish general multirate accuracy or nonlinear stability. Nondividing,
asynchronous/adaptive schedules and mid-step checkpoints require a later phase
design. Extending specialized EMPIC1D to a second polarization also requires
additional fields/current, diagnostics and persistence; removing the planar
check is insufficient.

The spectral derivative's Nyquist null, finite-precision FFT conditioning and
documented range boundaries remain. No exact discrete total-energy theorem,
universal stability certificate, global optimization or CAP-008 completion is
claimed by either source stage.
