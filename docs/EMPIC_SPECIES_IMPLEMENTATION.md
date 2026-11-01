# CAP-008 periodic electromagnetic species collections

Status: source implementation and independent root review are integrated; execution
and complete compiler acceptance are pending GitHub Actions. This is a synchronous,
nonrelativistic CPU `EMPIC` slice. The
specialized `EMPIC1D`, species-specific multirate schedules, relativistic
collections, collisions, nonperiodic boundaries and resident MPI/GPU evolution
remain open.

## Pre-edit brief (2026-11-01)

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

## Progress

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

## Remaining acceptance

The essential suite covers both precisions, dimensions and built-in shapes;
summed original currents, fixed/neutral/empty populations, original-factor
priming, shared old fields, synchronous fast-step equivalence, legacy two-species
compatibility, late arithmetic/preflight failures, grid/metadata aliases,
variable-step checkpoint/archive continuation, detached field permutation and
one strided destination roundtrip/count failure. These tests have not run
locally; the exact pushed revision must pass GitHub package, formatter, Aqua/JET
and neighboring checks, with complete compiler findings reviewed.

CAP-008 remains incomplete. Pair/unequal-species dispersion and convergence,
permutation/conservation campaigns, population growth and complete resource
acceptance, relativistic collections, species-specific multirate evolution,
EMPIC1D collections, collisions, boundaries and resident backends remain open.
