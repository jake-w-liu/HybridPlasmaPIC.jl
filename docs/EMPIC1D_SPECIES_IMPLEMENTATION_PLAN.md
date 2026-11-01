# CAP-008 specialized EMPIC1D species implementation

Updated 2026-11-01. **S1 is in implementation; no new solver behavior or runtime
acceptance is claimed here.** The [audit tracker](DEEP_DEBUG_TRACKER.md) owns
publication and CI status. This plan expands the specialized-model and multirate
items in [CAP-008](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-008--arbitrary-kinetic-species-in-esem).
The [generic EMPIC collection contract](EMPIC_SPECIES_IMPLEMENTATION.md) remains
the source of its existing scope and pending acceptance.

## Scope and source owners

`src/models/full_pic.jl` carries a specialized `EMPIC1D` with `Ex/Ey/Bz`,
longitudinal continuity current, midpoint transverse current, staggered magnetic
and velocity phases, and legacy electron/optional coarse-ion stepping. Its
collection overload is missing. Generic `EMPIC{1}` already carries both
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
| S1 construction and orchestration | In implementation | Ordered populations gather the same old fields, then sum currents before one specialized field update per fast interval |
| S1 diagnostics and persistence | In implementation | Physical moments and verbatim field restore, complete graph detachment, saved-mode validation, identity/capacity handling and layout rejection |
| S1 source review and essential checks | Pending | Legacy helper-expansion comparison; independent current/phase and failure/retry/restore controls; existing assertions retained |
| S1 GitHub acceptance | Pending | Exact revision on minimum `1.12` and moving stable `1`; full neighboring tests, complete compiler review and formatter |
| S1 resource acceptance | Pending | Remote complete init/step/restart cost versus specialized legacy path, actual hardware/runtime/threads, allocation traffic and peak memory distinguished |
| S2 per-species schedules | Planned separately | Dividing integer strides, independent species phase/history, explicit restart version/migration and physical refinement |
| Broader CAP-008 integration | Open | Collisions, boundaries, partition ownership, device residency and matched physical/cost validation |

All Julia tests, compiler checks, formatter, builds and profiling run on GitHub
or suitable remote hardware. Local work is source/caller review and lightweight
file checks. Add compact essential regressions through the existing test
framework and `electromagnetic-species` selector; do not add a test campaign or
replace neighboring checks with mocks.

## S1 synchronous time and current contract

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

Use a distinct schema-1 `:empic1d_collection` payload with required
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

Define the represented coarse interval as `H_s=dt/(M/r_s)`. Its kick spans the
previous and new species half intervals; B centering retains its independent
fast history. Persist each species' previous interval and schedule semantics,
and report its own diagnostic velocity offset. Allow schedule changes only at
completed external boundaries where all blocks close. Zero steps preserve all
histories. Before coding, record a schema version and explicit S1 migration or
strict rejection policy; never infer species phase from fast dt alone.

Gate all-fast identity, actual legacy coarse-ion comparison, variable intervals,
restart and independent pair/unequal-mass refinement. Continuity for the split
path does not establish multirate accuracy or nonlinear stability. Nondividing,
asynchronous/adaptive schedules and mid-step checkpoints require a later phase
design. Extending specialized EMPIC1D to a second polarization also requires
additional fields/current, diagnostics and persistence; removing the planar
check is insufficient.

The spectral derivative's Nyquist null, finite-precision FFT conditioning and
documented range boundaries remain. No exact discrete total-energy theorem,
universal stability certificate, global optimization or CAP-008 completion is
claimed by either source stage.
