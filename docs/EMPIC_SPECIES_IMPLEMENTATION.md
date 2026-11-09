# CAP-008 periodic electromagnetic species collections

Persistent nondividing/adaptive species scheduling is a source candidate in
`empic_collection.jl` and schema3 of `empic_collection_checkpoint.jl`; see
[persistent schedules](#persistent-nondividing-species-schedules). Its compact
acceptance source has not run. The historical results below do not accept this
new mode; numerical, compiler, formatter, scientific and resource gates remain
unverified for it.

Status: synchronous CPU `EMPIC` collections are implemented with optional
relativistic motion. Both slices are integrated after independent root source
and caller review; their exact-revision GitHub runtime/compiler/format checks
pass at `b5ce21a` (run 41683122086). Specialized planar `EMPIC1D` S1 collections
are source implemented and passed the same run; see the S1 section below.
Specialized S2 dividing per-species schedules are source implemented with
separate history and schema-2 migration and passed the same run; resource
acceptance remains pending.
Generic D1–3 dividing kick schedules are now source implemented as an explicit
constructor opt-in, sharing the S2 history owner; see the generic schedule section
below. `test_empic_species_schedule.jl` now covers them: an independent rational
Boris oracle over due/fast events for D = 1–3, bitwise equivalence between an
all-fast opt-in and the unscheduled collection, setter/copy/alias/corruption
guards, schema-2 checkpoint/archive round-trips, schema-1 migration to all-fast
rows, scheduled cold-mode refinement and relativistic failure rollback.
Complete-operation cost acceptance is deferred; an indicative measurement
(Float64, 8-cell 1-D grid, 2×64 particles, one `step_empic!`, Julia 1.12.7)
recorded 98,448 B/1,324 allocations unscheduled, 113,328 B/1,617 with all-ones
strides at `n_sub=4`, and 119,984 B/1,877 with `species_stride=[2,3]` at
`n_sub=6`. Broader collision composition, nonperiodic boundaries and resident
MPI/GPU evolution remain open.

A first [electromagnetic–local-Coulomb composition](COLLISIONAL_EMPIC.md) is now a
source candidate for all-mobile ordinary-Vector generic collections with fixed
dt, `n_sub=1`, nonrelativistic motion, no species schedule owner and consistent
collision normalization.
It stages fields, scattering and advanced RNG before publication and has a
distinct checkpoint schema. No tests were added or run; runtime acceptance and
broader collision intersections remain open.

## Ordinary-store admission (2026-11-01)

The canonical particle owner now rejects known read-only `Base.Experimental.Const`
and atomic-only storage through standard wrapper parents before collection alias
checks. It covers every population and particle leaf. Fixed populations are also
copied at publication, and checkpoint restore writes weights and provenance.
Writable ordinary vectors, non-atomic Memory and writable views/reinterpretations
remain supported. Electromagnetic collection particle leaves must also be
injective: a strided view with multiple logical entries sharing one backing
cell rejects before staging or restoration. Host storage must support ordinary writes and expose aliases;
OS protection changes and arbitrary custom setter side effects are outside the
publication guarantee. The shared check also serves ES collections and closed
Yee models. Essential constructor/init/restore controls are present but unrun;
exact-revision GitHub runtime/compiler/format and resource acceptance remain open.

## Initial nonrelativistic pre-edit brief (2026-11-01, historical)

- Mode: extend the existing periodic electromagnetic PIC owner.
- Requirement: ordered Float32/Float64 CPU populations in 1D–3D, arbitrary finite
  charge and positive mass, mobile/fixed populations, one summed charge/current
  and Maxwell solve, synchronous substeps, diagnostics and deterministic restart.
- Baseline: `full_pic.jl` restricts electrons to q=-1,m=1 and optional ions to
  q=+1,m=mi. `_charge_density_sum!` already sums arbitrary charges before range
  replay; `_correct_current_continuity!` owns the longitudinal correction and
  pure-Nyquist exclusion. ES collections already own explicit neutrality policy.
- Gap: EM collection ownership, current accumulation, phase, staged publication
  and schema are missing. This extends the documented model scope.
- API: `EMPIC(g, species; mobile=trues(length(species)), background_charge=0,
  neutrality=:require, charge_rtol=64eps(T), c=5, shape=CIC(), n_sub=1)`;
  `init_empic!(state, species)`, `step_empic!(state, species, dt)` and
  `electromagnetic_diagnostics(state, species)`. Each mobile population receives
  one kick and drift at every Maxwell substep. Fixed populations retain positions
  and velocities and deposit charge only; their stored velocities must be zero.
  `relativistic=true` is rejected in this slice.
- Invariants: finite independently stored populations; original q,m reach the
  canonical factor owner even if q/m is not representable; exact physical charge
  checks use the ES tolerance/background policy. Initialization binds ordered
  population identities. All mobile velocities and B share the fast half level;
  variable full dt/n_sub uses the preceding fast interval. Fixed particles have
  no velocity phase. Reordering/replacement requires initialization with physical
  full-time velocities; growth requires velocities at the existing half level.
- Canonical owners: the existing charge deposition, gather, Boris, priming,
  midpoint drift, continuity correction, Gauss projection, Maxwell update,
  field-energy and particle diagnostics. Extract the Maxwell centering/update
  blocks without changing their expressions; both legacy and collection paths
  call those owners. A collection workspace stages particles and all arrays in
  `_empic_transaction_fields`, plus time/step/full/fast timestep history. Stage
  fields use the same grid/plans. Public physical fields, particles and phase
  publish only after the whole operation succeeds; FFT/gather scratch may change.
- Compatibility: integer-count electron/ion constructors and payloads retain
  their existing behavior. The concrete `EMPIC` type gains a workspace parameter.
  A collection is a distinct schema-1 payload with ordered particles, mobility,
  neutrality/tolerance/background, c/shape/n_sub, initialized flag, physical
  fields/current/charge snapshots and full/fast phase. Scalar/collection payloads
  cannot be interchanged. Destination policy/mobility must match; saved species
  values and EM configuration restore after all validation.
- Strategy: add a collection workspace to `EMPIC`, stage existing canonical
  kernels, sum all mobile midpoint currents before one continuity correction,
  and replay all original contributions if their sum loses range. Count shared
  field energy once. Extend checkpoint/archive through their existing owners.
- Verification: direct source/caller review and exact extraction comparison;
  focused compatibility, arbitrary-charge current/continuity, fixed/neutral/empty,
  variable-step, ownership/error and checkpoint/archive regressions. Run Julia,
  package, Aqua/JET and formatter checks only on GitHub. No local computations,
  benchmark or measured optimization claim.

Plausible wrong patches: an electron-only current fails the arbitrary-charge
reference; updating fields between populations fails the shared-field and legacy
compatibility checks; updating fixed populations fails the fixed-state snapshot;
resetting fast phase fails variable-step checkpoint continuation; a late input
failure after earlier publication fails the complete-state rejection snapshot.

## Initial slice progress

- [x] Inspect the canonical ES/EM, deposition, phase, diagnostics and persistence owners.
- [x] Record scope, compatibility, ownership and schema before implementation.
- [x] Implement canonical collection construction, initialization and advance.
- [x] Implement diagnostics and checkpoint/archive.
- [x] Add focused essential regressions and update public documentation.
- [x] Review the exact source diff and prove unchanged legacy Maxwell extraction.
- [x] Complete independent root review and integration.
- [ ] Pass exact-revision GitHub checks and complete compiler review.

## Source review and representation boundary

The collection stages all 24 arrays in `_empic_transaction_fields` and all
particle components before advancing. Only positions and velocities change;
weights, q, m, ids and tags remain unchanged. The exact conserved charge target
is therefore computed once per public call and reused at both endpoints of each
substep. This removes repeated charge-total scans by construction; complete-call
cost and memory acceptance remain unmeasured.

The physical inventory is E, B, J, rho_n and rho_np1. The midpoint source Jraw
and gathered magnetic input Bc must also be finite after their construction.
`ne` is per-population density/current scratch, not a collection electron
density. A number-density scratch cell may overflow when its complete
q-weighted charge/current is finite (including q=0). `_charge_density_sum!` and
`_add_species_current!` overwrite it before every use; exact replay reads
original particle operands. `curlE` and `curlB` can likewise overflow when the
canonical cold derivative representation produces a finite complete Maxwell
update. Each is overwritten by `_em_curl_parts!` before its next consumer;
centering/updates use the returned complete representation. No catch suppresses
operator errors or substitutes a finite field. General FFT conditioning and
unrepresentable required mesh outputs remain outside this acceptance.

Diagnostics and checkpoint/archive consume the physical inventory, excluding
these scratches. Restart restores required E/B/J/rho_n/rho_np1 and full/fast
phase; the next substep reconstructs derivatives, Bc, midpoint current and charge
before consuming them. All public/staged fields, particle/gather storage and
Fourier buffers plus midx/kfull/ik/kvec metadata must be independent. The staged
grid remains identity-bound to the public grid; replacing it requires a new
collection state.

Source comparison against `d6d171a8ae353e6f72b75004d92b76db27eff92f` checks the
old magnetic centering and Maxwell update blocks byte-for-byte. Expanding their
new helper calls reconstructs the entire legacy substep byte-for-byte. Initial
magnetic priming retains its original expressions and order in one shared helper.
The legacy allocation body differs only in its workspace type parameter and
extra field. Root independently reproduced both whole-function expansion comparisons and
verified all 17 frozen changed-file hashes before integration onto `7189d75`.
These are source proofs, not runtime compatibility results.

Restore review confirmed that incoming fields can be permuted live arrays when
using the internal in-memory state API. Copying them directly would overwrite
later source leaves. Collection restore now validates and detaches incoming
state before any publication, following the collisional checkpoint owner;
ordinary Vector capacity and owner binding capacity are reserved beforehand.
Fixed-size writable strided views retain their storage and require matching
counts. Expected input and arithmetic failures preserve public physical state;
resource-exhaustion recovery and custom array side effects are not established.

The candidate that saving might bypass finite physical-field validation was
discarded: `_save_checkpoint` validates the payload before atomic serialization,
and `_archive_run` applies the same validation before archive publication. No
second writer-side finite-field validation pass or compiler suppression was added.

## First remote checkpoint-routing failure

The [first essentials run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41278909886)
at `213ed8c` reaches collection stepping and policy checks on both actual Julia
1.12.7 and 1.13.3, then fails collection save: the shared checkpoint entry selects its legacy required
fields (`mi`, `electrons`, `ions`). No persistence acceptance follows from the
earlier passing groups. Both failed raw logs are retained; corrected-revision
checks remain pending.

The shared EM entry now routes by payload layout and destination collection
workspace, following the existing ES checkpoint policy. A named private owner
retains every collection schema/storage/phase check; the scalar validation body
is unchanged. Scalar/collection layout mismatches reject before publication.
The existing checkpoint/archive, malformed-payload and strided-view cases must
pass remotely before this repair is accepted. Two public cross-layout load
controls also verify that both destination states remain unchanged.

## Remaining acceptance

The first relativistic essentials run at `7db86a9` passes the preceding
nonrelativistic checkpoint/archive and strided-destination cases, then fails
11 relative continuity assertions in the uniform fast-phase fixture on each
runtime. Physical fast/fine state comparisons pass. The subsequent
[scale observation](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41283466552)
records identical measurements on actual Julia 1.12.7 and 1.13.3: Float32
projected charge increments are about 6e-8–7e-7 with absolute residuals at most
5.75e-10; Float64 nonzero increments are about 1e-16–1.1e-15 with residuals at
most 1.11e-18. One exactly zero increment has exactly zero residual.

The uniform fixture has zero physical density change. Its absolute residual is
now checked against one epsilon of the noncancelling population source-density
scale. The unchanged `1024eps(T)` relative bound applies to a separate mobile
nonuniform-density fixture, with an explicit resolved-increment guard and the
same variable fast-step schedule. The correction and public diagnostic formulas
are unchanged. At `82a7499`, the [focused GitHub run](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41285023160)
passes all 194 assertions on each actual Julia 1.12.7 and 1.13.3: physical diagnostics, fast/legacy
phases, resolved continuity, C1 rollback and checkpoint/archive continuation.
The 12 uniform measurements are identical to the earlier failed observation;
all 24 observation records also match across runtimes, including the resolved
signal guards and unchanged relative bounds. Full-package/compiler/format and
later-revision acceptance remain pending; the
earlier failures remain retained.

The initial essential suite covers both precisions, dimensions and built-in shapes;
summed original currents, fixed/neutral/empty populations, original-factor
priming, shared old fields, synchronous fast-step equivalence, legacy two-species
compatibility, late arithmetic/preflight failures, grid/metadata aliases,
variable-step checkpoint/archive continuation, detached field permutation and
one strided destination roundtrip/count failure. The first GitHub run at
`213ed8c` passed its earlier charge/current/substep/policy assertions, then failed
at collection checkpoint dispatch: the shared validator entered the legacy
electron/ion payload path. The canonical routing correction and relativistic
follow-on are now composed; that failed run does not accept either change.
No tests ran locally.
The corrected exact revision must pass GitHub package, formatter, Aqua/JET and
neighboring checks, with complete compiler findings reviewed.

CAP-008 remains incomplete. Pair/unequal-species dispersion and convergence,
permutation/conservation campaigns, population growth and complete resource
acceptance, exact-revision acceptance of relativistic collections and specialized
EMPIC1D S1/S2 collections, broader species multirate evolution, collisions,
boundaries and resident backends remain open. S1's synchronous source path
does not by itself complete S2's dividing schedule; the S2 source candidate below
still requires integrated acceptance.

## Relativistic follow-on pre-edit brief (2026-11-01, historical)

- Mode: extend the synchronous collection owner and canonical diagnostics.
- Requirement: optional relativistic mobile populations with the unchanged C1
  supported velocity domain, common fast phases and staged publication; physical
  relativistic particle momentum/energy; schema-1 checkpoint/archive continuation.
- Baseline: `empic_collection.jl` rejects relativistic construction/validation,
  hard-codes nonrelativistic priming and Newtonian diagnostics. Its staged
  `_push_velocities!` already dispatches on the ordinary EMPIC relativistic flag.
  Collection schema 1 has a Boolean `relativistic` leaf but accepts only false.
  Scalar energy's existing cold owner sums rational radicals via integer
  square-root bounds and can extend across populations.
- Gap: enabling only the flag leaves Newtonian diagnostics and omits the
  relativistic configuration/domain on restart. Native scalar energy forms a
  cancelled `1-v²/c²` until an eight-epsilon domain guard. The existing priming
  residual-ball owner already encloses gamma from original stored inputs and
  includes division, square-root and arithmetic residuals.
- Invariants: charge/current, mobile/fixed, exact neutrality and shared Maxwell
  contracts remain. Mobile populations use canonical relativistic prime/kick and
  C1 stored-output checks, with original q,m. Rejection cannot publish earlier
  populations. Stored half-time diagnostics use sum(w*m*gamma*v) momentum and
  sum(w*m*(gamma-1)*c²) kinetic energy, counting shared field energy once.
- Canonical owners: existing prime/kick/domain/field staging; extend the exact
  radical energy owner to whole-collection energy and momentum, and reuse the
  bounded residual-ball owner for physical gamma and complete moments. Complete
  original factors replay on range or unresolved conditioning, including
  cross-population cancellation and subnormal contributions. No new MPFR hot loop.
- Compatibility: default/nonrelativistic behavior stays unchanged. Scalar
  relativistic low-speed energy order remains; high-speed contributions require
  a physical error enclosure. Schema 1 retains its layout, validates the saved Boolean
  and c against saved particles, and restores them; old false payloads remain
  readable. Saved mode replaces the destination mode.
- Dependency: P1 concurrently repairs normal relativistic kick accuracy in the
  canonical mover. This patch does not edit `_boris_rel*`; final accuracy
  acceptance requires the combined source and GitHub results.
- Strategy: thread flags/domain into existing collection prime/kick; extend
  shared diagnostic/persistence owners; retain open multirate, EMPIC1D, collision,
  boundary and backend intersections.
- Verification: essential both-type rational-gamma/near-light diagnostics,
  original-factor range/cancellation, fast-phase and legacy equivalence, C1
  rollback and changed-mode/c checkpoint/archive continuation. Direct source,
  caller and final-diff review; all Julia/formatter/package/compiler work on
  GitHub, with no local runtime execution.

Wrong patches rejected: gamma-free diagnostics fail rational-gamma references;
per-population rounded range sums fail whole-collection cold completion; enabling
only the public flag fails phase/restart controls; unsupported rounded output
acceptance fails C1 atomic rejection.

- [x] Inspect source and record the relativistic contract before production edits.
- [x] Extend canonical diagnostics and collection evolution/configuration.
- [x] Add essential phase/domain/range/restart regressions and update docs.
- [x] Complete independent root source/caller review and integration.
- [ ] Pass exact-revision GitHub essentials, full CI and complete compiler review.

## Relativistic source boundary and acceptance

The optional mode reaches the same canonical `_prime_leapfrog!` and
`_push_velocities!` owners as legacy relativistic EMPIC. No prime, Boris, Maxwell,
current or continuity numerical kernel is copied. The existing C1 domain is
checked at collection input and after each staged mobile kick/drift. Every
population shares the existing fast interval and variable-step centering.
Nonrelativistic constructor defaults and diagnostic formulas remain unchanged.

`electromagnetic_diagnostics` chooses physical relativistic moments when the
mode is enabled. The new collection `kinetic_energy_relativistic` and scalar or
collection `total_momentum_relativistic` require finite CPU Float32/Float64
mass, weight and velocity sources; signed finite factors have their algebraic
meaning. These read-only helpers accept the physical stored-input `|v|<c`
domain, which is wider than C1's mover-output domain. They preserve the stored
velocity time level and do not synchronize fields.

Whole-collection energy and momentum use the existing bounded residual-ball
arithmetic, including original m, w, v and c. The energy identity
`(gamma-1)c² = |v|² * gamma * (gamma/(gamma+1))` avoids a small-speed subtraction
without squaring a large gamma. The final result is accepted only when the
whole-sum enclosure lies inside one rounding cell. Otherwise the existing exact
rational-radical completion rounds the complete physical expression, allowing
cross-population cancellation, original-factor range rescue and subnormal sums.
Structural zero components remain exact zero. This is source-supported
arithmetic ownership; resource and complete-call cost acceptance remain open.
The canonical nearest/gradual context gate precedes all residual-ball work;
other arithmetic contexts use the same exact original-source completion.

Scalar `kinetic_energy_relativistic(ps,c)` retains its existing low-speed
operation and summation order. At high speeds, a finite historical weighted
dimensionless contribution is retained only when an enclosure proves its
distance from `w*(gamma-1)` is at most four times the smaller adjacent spacing
of that contribution. The bound includes the reference gamma's own residuals,
subtraction and weight multiplication. An unaccepted historical value uses a
certified rounding cell or complete exact replay. This is a point-contribution
bound; it does not claim a new whole-scalar-sum accuracy bound. Nonfinite mass
or weights keep their historical arithmetic and propagation. The unchanged
existing scalar energy corpus is registered in remote essentials alongside the
new focused collection tests.

Schema 1 already carried a relativistic flag. Its validator now requires Bool
and checks saved particle velocities using the saved mode and c. Restore
publishes mode and c with the detached species/physical fields and saved phase;
the destination's prior mode/c can differ. The legacy global `mobile` flag must
remain false for collection states, including checkpoint destinations; mobility
is always per population. Existing strict physical-field, alias and detached
publication controls remain. Root's checkpoint dispatcher correction and P1's
canonical normal relativistic accuracy repair are composed in this implementation.
Their remote acceptance remains required.

The new essential suite checks rational-gamma moments, signed cancellation,
mass/weight range rescue, near-light conditioning, common fast phases,
FTZ replay preserving original subnormal factors,
normalized legacy compatibility, a nonvacuous late C1 rollback and changed-mode/c
checkpoint/archive continuation. No local Julia, formatter, build or benchmark
was run. Neither full CAP-008 completion nor performance improvement is claimed.

## Specialized planar EMPIC1D S1 synchronous baseline (historical)

`EMPIC1D(g, species; mobile, background_charge=0, neutrality=:require,
charge_rtol=64eps(T), c=5, shape=CIC(), n_sub=1, relativistic=false)` now has
an ordered CPU Float32/64 collection path using the actual specialized
`Ex/Ey/Bz` solver. See the [dedicated S1/S2 plan](EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md).
All original stored `v_z` values must be zero, including under FTZ; fixed
populations require zero stored velocity. Generic `EMPIC{1}` already supplies
both transverse polarizations. This specialized slice does not add `Ez/By`.

Each mobile population gathers the same old fields, kicks and drifts at every
fast interval `dt/n_sub`. Total charge and the physical DC longitudinal current
include every applicable population; midpoint `Jy` is summed before the existing
Nyquist treatment. First-step priming, preceding fast-step centering and zero-step
phase retention use the canonical specialized owners. Changing `n_sub` preserves
the previous resolved fast interval. S1 is synchronous: the unchanged legacy
count API retains its distinct coarse-ion schedule. Compatibility checks use
the common planar `n_sub=1` subspace.

Initialization binds ordered identities and requires physical full-time
velocities. In-place population growth or edits do not provide a source model:
new/modified velocities must have the declared existing half phase after
stepping, and external charge changes require physically consistent fields and
source handling. S1 does not implement charge injection or supply its current.
Do not extend continuity/Gauss claims to arbitrary external particle/field edits.
Inputs require independent writable CPU strided leaves and the shared canonical
ordinary-store admission; known Const/AtomicMemory wrapper ancestry rejects.
Both model families validate fixed velocities once through the shared
stored-IEEE-zero owner, including constructor, initialization, zero-step and
saved-particle admission. Specialized mobile populations additionally require
stored-zero `v_z`. Signed zeros remain valid. Essential pre-stored subnormal
controls and arithmetic-zero observations await GitHub on both Julia versions.
The source comparison alone does not establish a pre-patch runtime failure.

The complete external operation stages particles, all specialized transaction
arrays and clocks before publication. Expected input/arithmetic rejection leaves
public state unchanged; private scratch/FFT buffers may change and are refreshed
on retry. Ordinary writable array/view identities remain. Resource exhaustion,
OS-protected memory and arbitrary custom-array side effects are not new rollback
guarantees. Concurrent calls require independent states, particles and grids.

`electromagnetic_diagnostics` reports exact charge policy, physical stored-velocity
particle moments (relativistic when selected), shared field energy counted once,
and the explicit fast velocity/magnetic half offsets. No exact fully discrete
energy, momentum, Gauss storage or general stability theorem is claimed. The
spectral Nyquist null and finite-precision FFT/source range limits remain.

Distinct schema-1 `:empic1d_collection` checkpoint/archive state requires
`polarization=:ex_ey_bz` and stores physical fields/current/charge snapshots,
ordered particles, mobility/background/neutrality policy and full/fast phase.
Saved mode/c govern saved-particle admission. All incoming leaves are validated,
detached and capacity-reserved before verbatim restoration. Legacy scalar and
generic collection layouts cross-reject; archive container 4 is unchanged.
Matched writable strided destinations retain identity and require saved counts.
Derived gathers, midpoints, centered B and transaction scratch are not persisted.

Source implementation and essential independent references are prepared; exact
integrated GitHub numerical/compiler/formatter gates and complete-step cost
acceptance remain pending. No local Julia or resource measurement ran. The S2 source slice below extends
this synchronous baseline; broader CAP-008 backend, collision and boundary
intersections remain open.


## Specialized EMPIC1D S2 dividing schedules: source candidate

`EMPIC1D(g, species; n_sub=M, species_stride=...)` now gives each mobile
population a positive integer stride r dividing M. Default mobile rows are 1;
fixed rows are 0. Every specialized collection has concrete schedule metadata,
including all-fast and empty states. Generic `EMPIC` remains synchronous by
default; its new explicit opt-in is described below. The legacy count API retains
its scalar numerical owners.

For h=dt/M and H_s=dt/(M÷r_s), species s kicks at zero-based fast indices
0,r_s,2r_s,... . Every mobile population drifts and deposits charge/current at
every h, including intervals without a kick. All populations use unchanged old
fields before one Maxwell update. Each first-step velocity is primed by its own
H_s/2; B uses h/2. At a later event the kick spans saved old H_s/2 plus new H_s/2;
interior events use the new interval. B retains its independent fast history.

`set_empic1d_schedule!(state; n_sub=state.n_sub, species_stride=nothing)` changes
the future schedule at completed external boundaries. `nothing` keeps the
current strides. The complete proposal is copied/validated before either
configuration store; previous species/fast histories remain. Direct `n_sub`
edits must leave the current strides valid before init/step/diagnostics/save.
Zero steps validate and preserve phase. Fixed/step-zero histories use original
stored-bit zero; advanced mobile histories satisfy last_sub_dt<=H_s<=last_dt.
Physical arrays, particles, clocks and candidate histories publish only after
the complete external operation succeeds. Metadata must have independent
ordinary storage. Rejection may change private candidates and shared FFT work.

Diagnostics retain physical moments and one shared field energy, with each
mobile nominal velocity offset -last_species_dt[s]/2 and B offset
-last_sub_dt/2. They also return copied `species_stride`/`last_species_dt`.
These declared stage offsets are not exact timestamps relative to rounded
`time`: H_s need not equal r_s*h and the clock retains repeated fast additions.
No synchronization or energy projection is supplied. The planar, source-range,
Nyquist/FFT, external source-edit and concurrency limits of S1 still apply.

Specialized saves now use schema 2 and required integer `schedule_policy=1`,
ordered strides and species history, with the same tag/polarization and archive
container 4. Fully validated tagged S1 schema 1 migrates explicitly to all-fast:
advanced mobile H inherits the saved actual last_sub_dt, fixed/step-zero H is
zero. Reserved S2 fields in schema 1 reject. Saved M/strides replace destination
configuration together, independently of old requested knobs/history values.
Incoming complete graphs detach/revalidate before publication; physical fields
restore verbatim and matched view identities remain. Scalar/generic layouts and
generic schema 1 retain their meanings.

Compact independent interior-event/current/field references, actual legacy
coarse-ion comparisons, resolved uniform pair/unequal-mass refinements,
relativistic rollback/retry and schema/migration/continuation checks are source
implemented. All exact-integrated GitHub numerical/JET/Aqua/formatter and
complete-step resource gates remain pending. No local Julia ran. Dividing
schedules do not establish a general nonlinear stability, temporal-order or
exact fully discrete energy theorem. Persistent phases have a separate source
candidate below; mid-step checkpoints and broader CAP-008 integrations remain open. See the
[canonical S1/S2 plan](EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md) for acceptance.


## Generic D1–3 dividing kick schedules: source implemented

`EMPIC(g,species; n_sub=M,species_stride=rows)` opts into the same dividing
schedule/history owner as specialized S2, while retaining the actual generic
three-vector force, midpoint current, continuity and Maxwell kernels. D1 retains
both transverse polarizations; D2/D3 retain their full fields/velocities. Omitted
or explicit `species_stride=nothing` keeps the original concrete Nothing owner,
common-fast ordering, diagnostic keys and schema1. Explicit all-fast rows still
opt in. Persistent nondividing phases use the separate concrete schedule owner below.

Mobile strides are positive integer divisors of M; fixed rows are zero. Bool,
noninteger/out-of-Int values, wrong axes/counts and invalid divisibility reject.
For h=dt/M and H_s=dt/(M÷r_s), kick at zero-based j divisible by r_s. All mobile
species drift and deposit charge/current at **every h**, including intervals without a due kick;
current is not reused from an earlier field step. Only due species gather common
old E/centered B and call the canonical original-factor Boris owner. Common
magnetic centering/update retains previous fast h, independently of coarse H_s.

First-step priming uses H_s/2 for each velocity and h/2 for B. Later due kicks
join the saved old species half interval and new half interval; interior events
use the new H_s. The next completed boundary publishes each species' last used H
with all particles/fields/clocks. Failure leaves authoritative history unchanged;
private interval candidates may change. Zero dt validates but does not prime or
consume a proposal. Empty and fixed rows preserve their declared zero/nonzero
phase rules. Growth requires the existing per-species velocity half phase.

`set_empic_schedule!(state; n_sub=state.n_sub,species_stride=nothing)` atomically
copies/validates the future proposal on an opted-in generic state. Nothing keeps
current strides. It cannot enable a default Nothing-owner state; construct a new
opted-in state instead. The public specialized setter remains unchanged and
shares this proposal owner. Future configuration changes preserve previous
species/fast histories; do not infer old H from current knobs. Nominal diagnostic
offsets are -last_species_dt/2 and magnetic -last_sub_dt/2, with copied schedule
arrays. Rounded H_s need not equal r_s*h or exact public-clock timestamps.

Scheduled generic state uses tag `:empic_collection`, schema2, Int policy1 and
ordered stride/history, with the same complete generic physical payload.
A scheduled destination explicitly migrates fully validated generic schema1
without schedule keys to all-fast rows; advanced mobile H comes from its saved
actual last_sub_dt, even if future n_sub had changed. Schema1 with reserved
schedule keys and unknown versions/policies reject. Default destinations reject
schema2 instead of discarding its phase. Scalar and specialized tags remain
separate; archive container4 stays. Incoming saved M/stride/history replace old
requested settings together after graph detachment and metadata/capacity
preparation, with verbatim field/current restore and saved-mode particle admission.

`CollisionalEMPICStepper` explicitly requires no schedule owner, even for all-fast
M1, matching its existing common-fast collision and schema1 contract. Scheduled
collision composition remains open. Source reduction checks cover default and
specialized expressions; no runtime/JET/compiler/formatter, numerical or resource
acceptance is claimed. No tests or executions were added for this implementation.

## Persistent nondividing species schedules

`EMPIC(g,species; species_dt=Hs,event_budget=K,n_sub=M)` and the specialized
`EMPIC1D` collection opt into a persistent calendar. `species_dt` and
`species_stride` are mutually exclusive. Supply one finite resolved positive
cap per mobile population and stored zero for fixed rows. H may exceed the
requested macro duration, need not divide it and need not be an integer multiple
of a field interval. Nonzero scalar inputs that disappear in grid precision
reject. `event_budget` is an explicit positive Int bounding field segments plus
recorded kick decisions in each positive call; exhaustion rejects the whole
candidate. Typed empty populations keep the same calendar rules.

The exact rational source clock `source_time` orders base edges
`tau_start+j*Q(dt)/M`, pending species deadlines and the macro endpoint. Equal
deadlines coalesce; unequal deadlines remain separate. Each positive gap narrows
once to an actual field interval h and must survive the canonical resolved-dt
gate. Every actual repeated floating clock addition must advance. Tiny distinct
events reject; they are not merged or moved onto a tick grid. `last_dt` remains
the requested macro duration and `last_sub_dt` is the last actual field h.

Each active row carries its actual H, last kick a and next deadline b with
`b-a=Q(H)` and `a<=source_time<=b`. Its nominal velocity source timestamp is
`(a+b)/2`, so the reported offset can be positive or negative. A right-edge due
kick remains pending until the next positive call; zero dt validates and consumes
no event. The public floating clock retains repeated actual-h addition.
`public_clock_defect=Q(time)-source_time` and per-interval duration defects expose
materialization differences. These nominal source coordinates are not an exact
sample-time theorem for the rounded particle trajectory.

At step zero choose each initial H once on unprimed full-time values. Delegate
the row's H/2 velocity prime and the first actual field-h/2 magnetic prime, then
reuse the initial decision for its first kick. Later due kicks consume the
canonical represented old-H/2 plus new-H/2 duration; the receipt records that
actual duration and its exact-factor discrepancy. B independently retains the
actual previous field h. At every field segment, all due rows gather the same
old E/centered B before any field update, and **every mobile row** drifts and
deposits fresh midpoint current and charge. The existing original-factor Boris,
relativistic C1/output admission, source replay, DC/Nyquist and Maxwell owners
remain authoritative. Specialized `EMPIC1D` retains its planar restriction.

Use `set_empic_schedule!` or `set_empic1d_schedule!` at a completed boundary to
change future `n_sub`, `species_dt`, `species_adaptive`, `species_dt_min`,
`schedule_limits` or `event_budget`. `nothing` retains that value. Prepare and
validate the whole copied proposal before publication. A setter preserves
pending H/a/b and preceding magnetic h; the next due kick uses its new proposal.
It cannot switch between persistent, dividing and Nothing owners.

Optional `species_adaptive` is a per-row Boolean vector. Adaptive rows require
explicit positive resolved `species_dt_min` and a NamedTuple
`schedule_limits=(rotation=...,electric=...,travel=...,growth=...)` of positive
finite values with growth>=1. Fixed rows remain false/zero. At an actual due
sample, E1/B1/V1 are maxima of particle L1 magnitudes of gathered fields and
stored velocity. Starting at the cap, limited by growth*previousH after
bootstrap, halve H until all exact original-factor bounds hold:

```
H*abs(q)*B1 <= rotation*m
H*abs(q)*E1 <= electric*m*c
H*V1        <= travel*minimum(dx)
```

Reject below the minimum or halving budget. These are explicitly supplied
sampled-state resolution limits. They do not estimate local truncation error,
guarantee future bounds during a pending interval, or provide coupled error
control. No callback or RNG state is introduced. Full error control remains
required work within CAP-008.

New-mode source inventory is fixed between initializations: ordered counts,
q/m, weights, provenance and ordinary component identities, plus committed
particle/physical-field/clock values. Detached shadows detect external edits.
Reinitialize using physically full-time inputs after forcing, injection, growth,
replacement or other source edits; the implementation supplies no recentering.
Existing default/dividing growth responsibilities remain as documented above.
Candidates contain all calendar/controller/report/shadow changes, so rejection
does not consume committed events. Scratch and shared FFT work may change.
Concurrent calls require independent states, populations and grids.

`last_schedule_report` retains actual partitions, due decisions, gathered
resolution operands, both deposited-charge endpoints and actual E/current
operands per segment. Generic raw and used current impulses are exact-factor
`sum(Q(h)*Q(J))`; planar reports retain actual Jx/Jy without inventing a raw
longitudinal channel. Work uses actual stored currents and exact midpoint E.
Report validation recomputes these impulses/work, checks historical calendar
and half-join algebra, and links the final source fields to the saved endpoint.
`field_defect=field_after-field_before+work_used` is measured, not forced zero.
Continuity observations retain the canonical spectral null-space and FFT
conditioning qualifications; `charge_source_scale` retains the canonical exact
absolute physical charge of the fixed source inventory, separately from the
normalized residual. No exact total-energy, temporal-order, nonlinear
stability, momentum, passivity or performance theorem follows from this mode.

Schema3/policy2 and `phase_policy=:persistent_event` persist the complete future
proposal, exact source clock, H/a/b, actual previous field h, current/source
snapshots and report. Generic Jraw joins the new-mode physical inventory.
Capture/output consumes no due event. Raw and checksummed archives share the
canonical transport and restore owner. Restore admits stale overwritten old
values, checks the original incoming graph against every retained live role
before copying, builds an independent phase/shadow, and publishes saved species,
provenance, fields and clock verbatim. It does not prime, project or silently
rebase time. Old schema1/2 meanings and schema1-to-dividing migration remain;
loading across persistent and old owner kinds rejects.

The bounded report retains O(K*mesh) actual operands, and validation/shadow
copies add complete-operation work. Allocation traffic, peak memory and elapsed
time require remote measurement. The compact
`test_empic_persistent_schedule.jl` supplies independent fixed-calendar/Boris,
all-fast compatibility, adaptive-sample, all-row current, rollback and
checkpoint/archive controls. It is unexecuted locally. Exact integrated remote
runtime, compiler and formatter checks, coupled refinement/error campaigns,
resource acceptance, source-growth/collision/boundary intersections and resident
distributed/device evolution remain open. This source candidate does not close
CAP-008.
