# CAP-008 periodic electromagnetic species collections

Status: synchronous CPU `EMPIC` collections are implemented with optional
relativistic motion. Both slices are integrated after independent root source
and caller review; exact-revision GitHub runtime/compiler/format acceptance
remains pending. Specialized
`EMPIC1D` collections, species-specific multirate schedules, collisions,
nonperiodic boundaries and resident MPI/GPU evolution remain open.

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

The initial essential suite covers both precisions, dimensions and built-in shapes;
summed original currents, fixed/neutral/empty populations, original-factor
priming, shared old fields, synchronous fast-step equivalence, legacy two-species
compatibility, late arithmetic/preflight failures, grid/metadata aliases,
variable-step checkpoint/archive continuation, detached field permutation and
one strided destination roundtrip/count failure. The first GitHub run at
`213ed8c` passed its earlier charge/current/substep/policy assertions, then failed
at collection checkpoint dispatch: the shared validator entered the legacy
electron/ion payload path. Root owns the canonical validator-routing correction;
the relativistic follow-on depends on that composition. No tests ran locally.
The corrected exact revision must pass GitHub package, formatter, Aqua/JET and
neighboring checks, with complete compiler findings reviewed.

CAP-008 remains incomplete. Pair/unequal-species dispersion and convergence,
permutation/conservation campaigns, population growth and complete resource
acceptance, exact-revision acceptance of relativistic collections, species-specific multirate evolution,
EMPIC1D collections, collisions, boundaries and resident backends remain open.

## Relativistic follow-on pre-edit brief (2026-11-01)

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
