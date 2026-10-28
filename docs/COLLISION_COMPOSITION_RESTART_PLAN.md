# CAP-010: collision/field composition and stochastic restart

Status (2026-10-26): ported onto `main` after `cda63f1` on the landed cell-local
Coulomb operator; the provisional V10 collision copy is not used. Both recorded
open findings are resolved as test-expectation errors with evidence; see
[Port onto main](#port-onto-main-2026-10-26) at the end. The sections below are
the design and evidence history of the frozen candidate (design base
`a18f09df5c7aeff4c68c3f3252460197fd2fb7c5`, provisional V10 assembly `3e1584a`).
Evidence: `../HybridPlasmaPIC-audit-evidence/2026-10-25/collision-composition-v1`.
Case 35 remains frozen separately with its original Julia 1.12.7 RSS acceptance
failure; this work does not resolve that failure.

## Original design pre-edit brief

- Mode: evidence-backed design with executable, independent witnesses. No new
  sampler, collision kernel, field algorithm or production composition wrapper.
- Requirement: define a bounded first restartable collisional hybrid path, with
  explicit physical/discrete time levels, order, RNG ownership, state publication,
  failure semantics and an implementation/validation sequence.
- Baseline: `extrapolated_leapfrog.jl:3-29,739-807` documents and implements
  integer positions/B/carried E and half-time particle velocities. `last_dt`
  sets the next kick interval. The local Coulomb owner stages every velocity
  before publication but permits RNG/scratch advancement on failure. Current
  checkpoint capture stores deterministic state and timestep history, not RNG.
- Design gap: externally calling these operators in sequence supplies neither
  a synchronized physical-state map nor a whole-composed-step transaction.
  A seed cannot describe an already advanced random stream. These are source
  contracts; numerical consequences below require the named witness runs.
- Invariants: one field clock and ordered population list; no hidden particle
  growth/reordering; every accepted collision interval and RNG advancement
  belongs to one committed step. Invalid/mid-step restart states fail visibly.
- Owners: integrator orchestration, existing moment/Ohm/push owners, accepted
  local collision API, and canonical checkpoint capture/validate/restore plus
  archive publication. No checkpoint-only private reconstruction of physics.
- Compatibility: preserve existing collisionless step/checkpoint APIs and
  seeded standalone collision semantics. New composition must identify its
  own schema and RNG policy; old deterministic bytes cannot be inferred to
  contain a stochastic continuation. No global Serialization overrides.
- Strategy: first establish witnesses and contracts, then implement the accepted
  CAP-001 composition in separate reviewed slices. Do not modify held candidates.
- Verification: both Julia 1.12.7 and 1.13.2; exact restart/failure/zero-step
  state comparisons, independent time-level and split-order oracles, stochastic
  weak-error confidence budgets, complete-step resource/compiler gates and CI.

## Verified ownership map

| Concern | Current owner and source contract | Composition requirement |
| --- | --- | --- |
| State time levels | `src/integrators/extrapolated_leapfrog.jl:3-29,83-87,739-807` | Track integer field/position time and velocity half level; do not label stored v as synchronous |
| First kick / changing dt | `_prime_leapfrog!`, `_push_hybrid_population!`, `step!` at lines 562,666,766 | Priming is once; later force kick is `(last_dt+dt)/2`, not always dt |
| Field/current refresh | `_hybrid_integer_moments!`, `ohms_law!`, `_recenter_carried_E!` at lines 338,658 | A velocity-changing source can invalidate carried E; refresh with the owning half-level convention |
| Local collision | Held `src/particles/coulomb_local.jl:651-775` | Delegate accepted self/cross sampling, normalization, locality and staging; no duplicate sampler |
| Collision scratch | Held `CoulombCollisionWorkspace` at lines 13-24 | Derived bins/indices/velocities only; recreate after restart, never treat scratch as RNG state |
| Current transaction | `_save_pressure_transaction!`, `_restore_pressure_transaction!` at lines 194-239 | Present only for evolving pressure; algebraic hybrid fields lack full-step rollback |
| Checkpoint | `src/io/checkpoint.jl:3-15,660-713,961-1020,1142-1171` | Extend capture/validation/restore for composition, preserving array identities and atomic file publication |
| Archive | `src/io/restart.jl:28-59,92-112` | Seed remains provenance metadata; persist actual RNG state in validated state payload |
| Parallel/device | CAP-003/004, MPI checkpoint and cell-ownership contracts | Deferred; no serial-bitwise claim across species/rank/device/layout changes |

The exact read sources and hashes are retained under `baseline` and
`held-local-v10`. Line references are to that frozen base, not current moving main.

## Original bounded first-slice proposal

Recommend an explicitly **first-order** fixed-step CPU collision composition for
`HybridStepper`, ideal massless algebraic `IsothermalElectrons`
(eta=etaH=de2=0), built-in NGP/CIC/TSC shapes, ordered positive-ion
populations, periodic Fourier grids, and Float32/Float64 in D=1:3. Use the accepted
local Coulomb API with explicit `coulomb_log`, `epsilon0`, `weight_to_number` and
`u_floor`. Fixed species and particle counts are the first boundary. Neutral MCC,
BGK, ionization, radiation, evolving/anisotropic pressure, CAM-CL, ES/EM, MPI and
GPU remain separately planned stages.

A provisional `CollisionalHybridStepper` would own the existing HybridStepper,
collision configuration, a private explicit RNG, local collision scratch and
reusable transaction storage. Names are provisional until implementation review;
do not add a generic callback/plugin framework. The normal collisionless API
continues unchanged and pays no composition work/storage cost.

Proposed sequence for nonzero constant h:

1. Validate supported storage/configuration, independent populations, clock,
   required h, normalization, and completed-boundary ledger before mutation.
2. Save public particle/field/clock state and live RNG to reusable transaction
   storage; set a private in-flight phase. No user callbacks execute in this slice.
3. Call the existing deterministic `step!` once, retaining its documented
   staggered output and previous-step history.
4. Call the accepted local collision operator for interval h, at current
   positions and stored half-time velocities, using the privately owned RNG.
5. Recompute integer-position moments/Ohm field and invoke the existing carried-E
   recenter owner with h. Validate the complete resulting public state.
6. Commit the completed collision interval/phase ledger; return the original
   wrapper. Any failure restores the whole entry state and RNG before rethrowing.

This is a discrete Lie-style, first-order proposal with explicitly staggered
source evaluation. The end collision uses positions at integer time and velocities
at a half level; it is **not** a synchronous exact collision subflow. Its leading
coordinate/time-level discrepancy is O(h), so a symmetric-looking sequence alone
cannot establish second order. Actual weak convergence and self-consistent field
coupling remain required implementation gates. A passing manufactured drag
witness establishes only its stated reduced problem.

Reject timestep changes in this first public slice before mutation. A subsequent
variable-step implementation must derive collision time accounting together with
`(last_dt+dt)/2` kicks, abrupt changes and restart phase; it cannot merely use the
new field dt for every operation by analogy. Zero dt is a true no-op, including
RNG/phase/scratch, after validating the scalar dt contract; it must not prime.
A configured zero collision rate is a full physical negative control and must
match collisionless results. Resource tests must report its transaction overhead
separately; absent composition must preserve the existing zero-overhead path.

Second-order/Strang composition is a later slice. It needs a demonstrated
second-order *synchronized* deterministic map or a derivation native to the
staggered state, and an appropriately accurate stochastic transition kernel.
Conservation and two half collision calls do not prove weak order two. Resetting
`init!` every step and treating an Euler half-kick as an exact synchronization
are specifically forbidden until independent order witnesses support them.

## RNG, checkpoint and failure contract

- Accept an explicit `Random.Xoshiro` or `Random.MersenneTwister` initially, deeply copy
  it into owned state at construction, and retain the algorithm type plus full
  advanced state, including cached draws. Do not take ownership of the caller's
  RNG or mutate a task-global RNG. The installed stdlib `copy(::MersenneTwister)`
  shares its mutable seed field (`Random/src/RNGs.jl:115-132`); construction must
  detach that metadata too, while reusable private rollback uses complete RNG
  state copying. The witness includes an external-seed mutation control. Custom RNG/native callback serialization is
  outside the first contract and must reject explicitly.
- Reproducibility is exact for identical stored state, species order, algorithm,
  Julia version and environment. Test each supported version independently;
  do not promise bitwise continuation across Julia versions or platforms.
- Save only completed boundaries in the first schema. Persist inner deterministic
  checkpoint, collision configuration/normalization, RNG state/algorithm, fixed h,
  completed source interval count and phase `ready`. Reject an in-flight phase;
  do not invent an idle state from incomplete bytes. Future resumable partial
  phases require their own schema and persistent stage state.
- Validate every incoming component, source count/clock consistency, RNG type,
  configuration graph/ownership, species order and destination capacity before
  restoring any destination. Do not merge a stochastic payload into the old
  deterministic integrator tag. Recreate only derived collision scratch.
- Reuse canonical `_atomic_serialize`/archive checksum publication. Serializer,
  write, close and rename failures preserve prior checkpoints and propagate the
  original error. Do not serialize FFT handles or transient transaction buffers.
- Whole-step rollback covers particle x/v, all public field arrays/counters,
  time/step/last_dt, source ledger and RNG. Preserve existing array identities.
  Private scratch may change. Standard owned CPU buffers keep rollback bounded;
  allocation failure or arbitrary external callback side effects are not hidden
  behind a false transaction guarantee. No automatic retry changes random draws.
- First tests inject failure after field advance, after RNG draws, after one
  internal species pair, after collision publication, during carried-E refresh
  and before final commit. Check exact state/RNG equivalence and identical retry.
  The standalone collision API's weaker RNG rollback semantics remain unchanged.
- No public hooks/checkpoint calls while in flight; reject reentrancy/concurrent
  use. Independent simulations use independent wrappers/RNGs/workspaces.

## Premeasurement witness contract

These are diagnostics, not production acceptance or a competing collision model.
Retain exact source/env hashes, commands, logs, raw errors and every failed run.

1. Actual HybridStepper constant-force variable-step witness: stored v matches
   analytic `v(t-last_dt/2)` and x matches integer time to roundoff; zero dt
   leaves all observed state unchanged. A physical-time-v comparison must fail.
2. Actual uniform-B reinitialization witness: ordinary one-time priming plus a
   final predictor has order 2; repeated init/step/Euler-predict has order 1.
   For h=1/N, N=(16,32,64,128), require consecutive orders 1.8–2.2 and 0.8–1.2.
   The independent matrix identity is `(I+hA/2)R_h(I-hA/2)=(1+h²/4)R_h` for
   A²=-I; repeated application accumulates the amplitude defect.
3. Independent noncommuting linear subflows: exact rotation F and anisotropic
   damping C; compare to an independently closed-form exponential of F+C.
   Constant and alternating steps, N=(16,32,64,128): Lie order 0.9–1.1 and
   symmetric exact-subflow order 1.9–2.1. Double-duration-C mutation must fail.
   This witnesses order logic, not the Coulomb kernel's stochastic weak order.
4. Actual forced field step plus manufactured exact drag: compare x and stored
   half-level v with analytic acceleration/drag at integer/half times. Constant
   and alternating steps, N=(16,32,64,128), expected first-order errors 0.8–1.2.
   This probes the proposed scheduling boundary; the manufactured drag is never
   production code or a replacement physical collision operator.
5. Prescribed equal-weight, unequal-mass conservative pair rotation: verify total
   kinetic energy/momentum, demonstrate changed charge-weighted bulk velocity,
   unchanged stale carried E and corrected canonical refresh. This fixture has
   no random sampler and supports no scattering-rate claim.
6. Actual checkpoint and standalone collision continuation with both supported
   RNGs: advanced-state serialization must match complete particles and future
   mixed scalar/bulk/shuffle draws exactly; seed-only restart must disagree.
   Confirm existing deterministic payload has no RNG field.
7. Actual rejected collision sample: all particle velocities unchanged while
   RNG advances; retain exception and draw counts. This establishes the boundary
   the future composed transaction must strengthen, not a standalone API bug.

## Implementation stages and acceptance matrix

- [x] Read completeness, local collision, hybrid time-level and checkpoint owners.
- [x] Freeze baseline sources, this premeasurement witness matrix and the bounded
  first-slice proposal before diagnostic runs.
- [x] Run/independently replay the seven witness families on both Julia versions;
  retain failures and revise only unsupported design assumptions with evidence.
- [x] Parent review of first-order scope, time levels, fixed-step decision and
  source-ledger/RNG ownership. This is technical integration review, not a new
  user permission requirement.
- [x] Implement the bounded composition owner on the separate provisional V10
  dependency, with full rollback and exact disabled-source/component controls.
  This describes the historical V10 assembly; final sampler/phase recomposition
  and the landed port are recorded below.
- [x] Extend canonical checkpoint/archive schema and validation; exact same-env
  uninterrupted vs restart at every permitted boundary, corrupted/missing/type-
  mismatched RNG/phase/configuration, rejected incoming state and failed publish.
- [ ] Freeze and run actual composed physical gates: D=1:3, Float32/64, 1/2/4
  species, empty/singleton cells, unequal counts/weights/masses/charges, periodic
  crossings, weak/strong collision regimes and nonzero Hall field feedback.
  Require analytic weak-generator budgets and ensemble standard errors chosen
  before measurement; fixed seeds/common-path replay alone do not prove order.
- [x] Derive total mass/charge and equal-weight collision momentum/energy budgets
  separately from the existing deterministic integrator truncation errors.
  Use expectation budgets for unequal macro-weights, never exact conservation.
- [x] Complete failure injection at every owned phase, zero/invalid steps, limits,
  aliases, reentrancy, shared config and source-ledger tampering before promotion.
- [x] Freeze and run the complete 180-profile allocation and 56-profile retained-
  memory inventory. Original zero additional byte, 1 MiB Julia live and 16 MiB
  RSS limits remain unchanged; include all transaction copies and full RNG state.
- [x] Run targeted registered tests, actual Pkg test sandboxes, compiler/Aqua and
  formatting checks; independently replay raw physical, rollback and resource data.
- [x] Complete neighboring regression receipts and final documentation freeze.
- [x] Measure the complete 180-profile checkpoint/output allocation and size
  inventory; preserve its byte-stability failures and independent owner diagnosis.
- [x] Resolve the recorded output byte-stability limitation: the byte-identity
  expectation was wrong; bounded retention is the enforced invariant (port record).
- [ ] Acceptance latency requires a frozen exact-input campaign and continuously
  monitored exclusive window.
- [x] Repeat scientific/restart/resource/neighbor gates on `main`'s landed sampler
  (port record below).
- [x] Root integration and individual push as `2b2420b`; [all six exact-SHA jobs](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40326554556) and complete compiler review pass for the bounded fixed-step CPU isothermal slice.

Dependent expansion remains explicit: variable-step collision scheduling;
second-order stochastic splitting; evolving-pressure heat exchange; CAM-CL
persistent magnetic phase; ES/EM families; BGK/neutral reservoirs; ionization
particle growth, unique IDs and counters; radiation/antenna energy destinations;
MPI cell ownership/stream layout and genuine device RNG/storage. Each needs an
independent physical oracle, commit/restart semantics and complete-run cost gate.
No CAP-010 completion claim follows from this design or its witnesses.

## Witness iteration ledger

`witnesses-v1` retains 64 passes and four failed snapshot comparisons per runtime.
The snapshot deep-copied `fields.floor_count::Ref` and compared Ref identity;
printed counters and arrays match, and the source's zero-step return is before
mutation. V2 adds a minimal Ref identity control, compares the counter value,
checks preserved field identities and requires counter/velocity mutation detection.
No production source or physical/order gate changes. V1 also retains a harness
soft-scope warning: its receipt stays `running` after the failed testset. Exit 1
and failure logs are authoritative. V2 places status/receipt handling in a local
function and will verify deliberate failure status as a separate control.

The installed RNG source shows shallow MersenneTwister copies share the mutable
seed metadata. The proposed constructor now explicitly requires a deep copy;
an independent mutation control is added before accepting private RNG ownership.
This is a design refinement, not a change to any existing RNG implementation.

## Completed design witnesses

`witnesses-v2/commands.json` records exit 0 and **88 assertions per runtime** on
Julia 1.12.7 and 1.13.2. `replay-witnesses-v1.py` independently reconstructs the
numerical errors, proves the noncommuting Lie/Strang Taylor coefficients and
gyro gain identity with exact rational matrices, verifies archived source/env/log
hashes, and rejects 11 semantic mutations per runtime. Results are in
`witnesses-v2/independent-replay-v1.json`. The deliberately injected receipt
failure exits 1 with `status=fail` and its original cause on both versions
(`receipt-mutation-v1`); no failure is converted to a passing result.

Measured consecutive orders are 1.999989–1.999999 for ordinary gyro priming and
1.001391–1.005348 for repeated initialization/Euler synchronization. Exact Lie
subflows give 0.999148–1.001219; exact symmetric subflows give 2.000006–2.003979.
The actual forced step with manufactured end drag gives 0.991929–0.998120.
These reduced witnesses support the time-level/order design decisions only.
Actual local Coulomb weak convergence and self-consistent field coupling remain
unimplemented acceptance gates. Both supported RNG types reproduce continued
collisions and future mixed draws with saved advanced state, while seed-only
restart disagrees. The actual rejected sample preserves velocity arrays and
advances the RNG, confirming the stronger composed transaction is necessary.

Parent authorizes the bounded implementation after this design freeze. The ES
phase owner reports a provisional repair in
`HybridPlasmaPIC-nonrel-phase-20261025`: private push helpers now receive
`preceding_dt`, and `_predict_half_kick!` receives full dt. Its public `step!`
and `_recenter_carried_E!(...,full_dt)` contracts remain unchanged. Composition
will use these stable owners, not those changing private interval arguments.
The witnesses above are explicitly on a18f09d; their private predictor calls are
not claimed verified on the future phase dependency. Reassembly must reverify.

## Bounded implementation contract — frozen before production edits

Parent accepts the bounded interface and directs the frozen V10 collision
source as a provisional dependency. This is now the implementation contract:

```julia
composed = CollisionalHybridStepper(inner::HybridStepper, species::AbstractVector;
    dt, rng, coulomb_log, epsilon0, weight_to_number, u_floor=1e-3)
init!(composed, species)
step!(composed, species, dt; NB=1)
save_checkpoint(path, composed, species)
load_checkpoint!(composed, species, path)
archive_run(path, composed, species; rng_seed, kwargs...)
```

- The constructor requires an unadvanced inner stepper with zero time/step/
  last_dt. It validates the ideal isothermal model and built-in NGP/CIC/TSC
  shape, D=1:3 and Float32/64, standard CPU Vector particle storage, independent
  populations, fixed counts, positive ion charges/masses and explicit collision
  normalization. It stores the ordered population object identities and their
  counts/charges/masses. Replacing/reordering species or changing those physical
  identities must reject before mutation. Copy the species collection and
  collision matrix; deeply copy the supplied supported RNG. Constructor work
  must not mutate the inner state, particle values or caller RNG/configuration.
- Initialization is explicit. Active `step!`/save/archive reject an uninitialized
  wrapper. `init!` delegates to canonical initialization under the wrapper's
  full transaction, sets phase ready, and resets the composed clock/count.
  It **preserves the current RNG state**; reinitialization never silently reseeds.
  The first API does not wrap an already advanced collisionless run, since that
  needs a separate physical initialization/phase contract.
- A nonzero step must equal the fixed dt exactly at grid precision; reject a
  nonzero real dt that converts to zero. Existing scalar dt/NB validation owns
  numeric conversion. Zero dt is a true no-op and consumes no RNG or priming.
  Enforce species identities/order/counts, closure/backend/configuration and
  ready phase before active work. Track source time/count separately and require
  exact agreement with the inner completed clock before stepping or saving.
- Use canonical `step!`, public local `collide_coulomb!`, then existing integer
  moments, Ohm's law and `_recenter_carried_E!(inner,species,full_dt)`. No private
  predictor or kick-interval arithmetic is copied. If every channel is disabled,
  skip collision and post-collision refresh; results/RNG match the collisionless
  call exactly. The wrapper still owns its declared transaction buffers; report
  configured-inactive overhead separately from the absent-wrapper baseline.
- Snapshot x/v, every public HybridFields value/counter and inner clock into
  reusable storage; snapshot full owned RNG state. Keep shape/model/metadata
  immutable across the active operation. Restore entry state and RNG in cleanup
  on any exception, including interrupts, then preserve the original exception.
  Scratch may change. No callback hook, silent retry, random reseed, clipping,
  test-only branch or custom sampler enters production. Zero/failed steps leave
  the source ledger and phase unchanged. Full-step finite-result validation
  catches unrepresentable proposals before successful return.
- Add composition checkpoint dispatch beside canonical Hybrid owners, using
  existing deterministic capture/validation/restore, configuration graph and
  atomic publication helpers. A distinct versioned tag stores inner state,
  fixed dt, phase/initialized flag, source time/count, collision parameters and
  complete RNG state plus algorithm/runtime identity. Only ready initialized
  boundaries are serializable. Old deterministic payloads and mismatched RNG,
  runtime, phase, configuration/layout or fixed population counts reject before
  destination mutation. Preserve current public array identities; recreate only
  derived workspace. Validate RNG storage/cache dimensions/indices through the
  installed stdlib contracts before random use/copy publication. Constructor
  and incoming-state ownership must not retain live links into caller payloads.
- Save/restore works for fresh destination wrappers of the same dimension,
  precision, population identities/counts/charge/mass layout and RNG algorithm.
  It does not resize populations in this first slice. Archive delegates to the
  same state owner and existing metadata/checksum publication; `rng_seed` remains
  provenance, not a replacement for advanced RNG state. Generic unrelated I/O
  APIs and existing deterministic checkpoint formats are unchanged.

Before promotion, freeze/run 144 ordinary resource profiles: D=1:3,
Float32/64, 1/2/4 species, 0/1/256/4096 particles per species, and both supported
RNGs on 4^D grids. Add 36 symmetric-matrix-log profiles at 256 particles with
those dimensions/precisions/species/RNGs. Use eight warmups, full GC setup, three
stable identical-input allocation samples and zero measured compilation. New
orchestration must add **0 warmed allocated bytes** beyond the exact component
sequence; transaction copies reuse storage. Record the full component and
composed bytes, not only a helper. For largest fixed-state profiles, require
Julia live growth <=1 MiB and RSS growth <=16 MiB after 20 further complete
calls and full GC; keep raw growth and explicit retained-buffer sensitivity
controls. These budgets precede implementation measurements and cannot be
expanded after failures. No latency claim is accepted without a separately
frozen exact-input campaign and monitored exclusive window.

Implementation tests must reject three plausible wrong patches: a seed-only RNG
snapshot (restart draws diverge), collision-only rollback (field/clock remain
advanced), and stale-E reuse (next coupled force differs). They must also reject
reordered populations, a changed nonzero dt, hidden phase loss, duplicate source
intervals and silently accepted malformed incoming RNG state. The independent
scientific weak-order/physical gates above remain prerequisites, not satisfied
by this interface or the reduced design witnesses.

Dependency assembly is a separate local commit/patch containing exact V10 source
hashes from `composition-v10-permanent-source.json`. Root owns its eventual
sampler assembly with the accepted null rotation plus pending angle/weight-mean
repairs. The V10 legacy sampler predates those repairs; this provisional assembly
must not be mistaken for main-compatible numerical acceptance. The composition
patch is measured relative to that dependency commit, and all scientific,
restart, resource and neighboring gates must repeat on the final root assembly
with the separately pending ES phase repair.

Implementation boundary detail (before code): composition checkpoints will store
all public HybridFields arrays and floor_count, in addition to the canonical
inner checkpoint. This preserves immediate post-load diagnostics exactly. The
existing Hall field-name owner will delegate to a shared HybridFields descriptor;
no deterministic payload changes. B/E duplicated by the inner record must agree.
Incoming configuration must match the constructed wrapper, and full RNG state is
deep-copied before publication so no incoming metadata links survive restore.
The installed Random sources on both supported versions retain the same
MersenneTwister storage/cache contract, although 1.13 moved it to
`Random/src/MersenneTwister.jl`. Validation will use its native constructor/cache
bounds before any random call. Only live Float64 cache entries must satisfy their
native [1,2) interval; exhausted cache bytes are retained verbatim. Xoshiro's
all-zero four-word recurrence is rejected (it stays zero under the actual shift/
xor update and is not a random stream). The fifth split state remains unrestricted.
The scalar timestep owner currently permits conversion underflow, so this wrapper
must explicitly reject nonzero dt converted to zero as required by its contract.

Implementation V1 is retained under `functional-v1`: each runtime reports
356 passes and 40 errors, all errors arising in the sparse diagonal-load
multidimensional field histories (including direct canonical component runs).
`preflight-v1` records zero-density cells, large public B/E growth and exact
exceptions for three dt values; it does not establish a new solver defect or a
physical stability limit. The V2 ordinary fixture covers every grid cell with
Cartesian mixed-radix positions. The original diagonal load is retained as an
explicit whole-state rollback/retry witness; no solver tolerance is changed.

The same preflight run positively reproduced two missing wrapper invariants on
both runtimes: replacing the inner grid and aliasing a public B array with a
particle position array were accepted. V2 binds the original grid and rejects
public field/field or field/particle overlap before any transaction mutation.
These are repairs to this unaccepted composition implementation, not changes to
the canonical field or collision operators. Concurrency/reentrancy uses an atomic
nonblocking acquisition; ordinary rejection preserves the original error and
releases ownership in `finally`. Zero dt remains an unconditional validated no-op.

`functional-v2` now passes 628 assertions on each runtime (1.12.7 and 1.13.2).
Coverage includes all D/T/RNG ordinary component-history combinations, exact
inactive-source/RNG equivalence, initial and advanced restart boundaries,
public array identity, grid/storage rejection and the retained failed-field
rollback/retry. This is functional evidence only: full malformed-payload,
phase-injection, physical weak-order, resource and final assembly gates remain.

### Allocation repair pre-edit brief

- Verified baseline: `resource-preflight-v1` fails all 24 samples per runtime,
  each adding 17,792 B with zero compile/recompile. Full states still agree.
  `allocation-owners-v1` measures 8,128 B in complete field validation and
  768 B in population validation, each called twice. RNG, model, boundary and
  transaction save/restore individually allocate zero. Full allocation traces
  locate heterogeneous tuple zip/flatten iteration and the repeated canonical
  `_validate_hybrid_species` call; the raw profile and inference are retained.
- Contract stays zero added warmed bytes and identical complete outputs/errors.
  Replace heterogeneous runtime tuple iteration with statically flattened
  homogeneous public-array tuples, preserving every dimension/finite/alias check.
- The constructor still invokes the canonical hybrid validation for positive
  charge, q/m representability and total-count overflow. Those values are bound
  and checked unchanged on every call. The existing local collision validator
  owns changing storage, axes, independence, positions, velocities and weights;
  use it directly in repeated wrapper preflight/postflight without redundantly
  allocating the hybrid count/identity validation. Canonical `step!` still calls
  its own full hybrid validation. No sampler/field-owner change or weakened guard.
- Reverify the complete byte comparison on identical states, the permanent
  functional/malformed/rollback cases, and invalid metadata/storage mutations.
  Preserve this baseline and do not subtract helper costs from acceptance.

`checkpoint-v1` retains 168 passes and five failures per runtime. The independent
`checkpoint-transform-v1` checks 11 assertions per runtime and confirms the two
unsupported test assumptions: Serialization.jl:193 emits a null pointer, and
serialized external arrays become detached storage. The canonical graph guard
rejects both original live payloads before serialization; the transformed
values are no longer forbidden links. Test those at their actual incoming-graph
boundary, with explicit transformation controls; production I/O policy is unchanged.

`sampled-failure-probe-v1` uses only finite inputs and stock supported RNGs. Of 32
fixed cases per runtime, 14 reject unrepresentable sampled velocities after the
canonical field step and RNG advance. All 14 composed operations restore the
full public state and RNG exactly; direct components expose both advances.
The 18 accepted proposals remain finite. These are actual exception-path
witnesses (46 assertions/runtime), not a production injection hook.

`resource-preflight-v2` passes 104 assertions/runtime: all eight dimension-1
precision/RNG/log-form profiles have three identical-input samples with exactly
zero extra complete-step allocated bytes, zero compile/recompile and identical
full outputs. `checkpoint-v2` passes 178 assertions/runtime. These resolve their
bounded preflight failures; the complete 180-profile resource inventory and
physical, quality, final dependency and integration gates remain open.


## Historical provisional assembly

The results in this section describe the isolated V10-based candidate before
its port to the landed sampler. The port and its verification appear below.
The bounded fixed-step CPU isothermal API landed as `2b2420b`, with six green
exact-SHA jobs and complete compiler review. The current contract and restart
example are in [COLLISIONAL_HYBRID.md](COLLISIONAL_HYBRID.md). Expanded physical
regimes, higher-order/variable-step composition, other closures and exclusive
latency acceptance remain separate open gates.

- `functional-v3`:817 assertions/runtime, including finite-input post-draw
  failure, exact retry and dynamic storage/value rejection. `checkpoint-v2`:178
  assertions/runtime. The final `sandbox-v1` actual Pkg test target passes 1,945
  assertions/runtime: 866 functional, 178 checkpoint, 896 analytic production
  invariants and 5 guide-example checks. This is not all-repository Pkg.test.
- `resource-full-v1`:2701 assertions/runtime; all180 profiles ×3 identical-input
  samples have exactly zero additional complete-step bytes and zero measured
  compile/recompile. `resource-memory-v2`:398 assertions/runtime,56 profiles,
  maximum raw Julia live4496B and RSS81920B, under the original1MiB/16MiB gates.
  Both64MiB retained-buffer controls detect the intentional retention. V1's
  memory harness parse error remains retained; no production memory workaround.
- `physical-angular-v1`:4,194,512 assertions/runtime; fixed8192 paths ×16 profiles
  compare actual public histories with independent55/90-digit angular quadrature.
  `replay-angular-v2.py` replays32 profiles, all retained first16 histories/profile,
  moment statistics and wrong-rate/phase/gyro/initialization controls.
- `physical-weights-v1`:655457 assertions/runtime; equal-weight conservation and
  unequal-weight expectation at fixed16384 paths/profile. Independent V2 replay
  reconstructs budgets from32 profile CSVs and rejects accept-both unequal-weight
  mutations. V1 failed to round parsed Float32 decimal text back to its declared
  precision; the independent witness measures3.404e-8 before that conversion and
  2.776e-17 after it. No production values or physical tolerance changed.
- `physical-feedback-v1` retains833pass/32fail per runtime: the h=.01 Float32
  next-force stale-E mutation signal falls below the fixed128eps threshold. The
  independent dense-DFT/NGP/predictor/Ohm oracle passes all288 actual states and
  rejects wrong curl/Hall/predictor mutations. A separately frozen h=.1 witness
  (`physical-feedback-v2`) passes865 assertions/runtime with unchanged limits;
  this stronger sensitivity witness does not erase the original failed gate.
  Its complete independent replay also passes all 288 actual states and the
  wrong-curl/Hall/predictor mutation controls.
- `quality-v1`: registered Aqua/JET/inference gates pass20 assertions/runtime.
  All20 concrete init/step/capture/validate/restore JET reports/runtime contain
  no errors. Complete package report trees in `package-jet-full-v2` match the
  reviewed V10 baseline (29 minimum / 27 latest abstract findings), after verified
  source-line/path normalization. Existing reports remain visible.

The reduced angular campaign supports first-order behavior in its uniform
production problem. It does not establish general coupled plasma weak order.
No acceptance latency campaign has run. The output inventory below is measured
but its original byte-stability gate remains failed. Neighbor regression passes 39,873 / 39,873 assertions (Julia 1.12.7 / 1.13.2). Final sampler/ES
phase reassembly and repeat gates, parent promotion and CI remain required.
Historical failed auditor versions (NumPy boolean JSON output) and raw data remain
alongside their corrected V2 replayers.


Final local verification details:

- `phase-injection-v2`: 544 assertions/runtime at eight owned phases with both
  ordinary exceptions and explicit InterruptException. Independent replay of all
  64 serialized entry/event/rollback records/runtime passes 586 assertions/runtime
  and rejects one-ULP field/velocity/RNG mutations. Original exceptions propagate;
  no process signal or production hook is used.
- Formatting preserves parsed AST and is idempotent: 20 assertions/runtime.
  The final formatted-source coupled allocation preflight passes 104/runtime,
  still with zero additional bytes. `git diff --check` passes.
- `final-evidence-review-v1.json` verifies 3,717 archived source/input pins across
  12 campaigns and maps every production file to final source, with only verified
  formatting differences. The final source audit found no additional confirmed
  composition-owned defect within this bounded implementation.
- The accepted main population-mean repair is not part of this frozen provisional
  V10 dependency. Parent owns final null/angle/mean/original-factor integration;
  all composed scientific, restart, resource and neighbor gates must repeat on
  that accepted sampler together with the pending ES phase repair.


## Expansion sequence after the bounded slice

These are planned implementations, not capabilities delivered by this patch.
Each starts with its own equation/time-level and failure contract; acceptance of
one row does not imply acceptance of the following rows.

| Stage | Canonical integration boundary | Required state and verification before implementation acceptance |
| --- | --- | --- |
| Final CPU dependency assembly | Accepted local `collide_coulomb!`, shared sampler and public Hybrid `step!`/carried-E owner | Repeat exact trajectories/RNG/restart, independent physical histories, all resource gates and neighbors; retain original rate-boundary failures |
| Variable-step CPU collision composition | Existing `last_dt` field/particle contract and this source ledger | Derive preceding/following collision intervals and velocity centering; abrupt changes, alternating steps, zero steps, restart at each allowed boundary and independent variable-step weak-order evidence |
| Higher-order stochastic composition | Native staggered deterministic map plus accepted stochastic transition kernel | Establish synchronized or explicitly staggered order conditions, positive substep intervals and weak-kernel accuracy; noncommuting coupled production oracles must distinguish first/second order before a Strang claim |
| Evolving/anisotropic electrons | Electron closure, moment/Ohm owner and canonical pressure transaction | Explicit receiving energy reservoir and heat-exchange equations; full pressure/source rollback and restart; analytic relaxation, entropy/energy and positivity gates |
| CAM-CL and ES/EM families | `camcl.jl`, `electrostatic.jl`, `full_pic.jl`, their canonical checkpoint owners | Respect persistent cyclic magnetic phase, each species' half-time history and fast subcycles; do not substitute Hybrid recentering. Re-derive current/charge continuity and source-time rules per family |
| BGK and neutral collisions | Existing `collide_bgk!` / `collide_neutral_mcc!` | Explicit bath state and exchanged momentum/energy; full RNG/reservoir continuation and rollback; independent equilibrium/rate tests and inactive-source cost controls |
| Ionization and particle growth | Existing `ionize_mcc!`, particle append/provenance and checkpoint owners | Stage all species, charge changes, unique-ID allocation and counters atomically; preserve growth/source clocks and RNG in checkpoints; threshold/rate, conservation, capacity failure and same-environment continuation gates |
| Radiation, antenna and Raycon deposition | Existing radiation/antenna/ray owners; receiving particle or electron-energy model | Freeze coupling time and energy destination, account for removed/injected power and source state; no relabeling of reduced radiation as a complete radiation model |
| MPI and devices | Accepted cell-complete ownership/migration, distributed checkpoint stages and native device streams | Specify same-layout exact versus changed-layout statistical continuation; collective error agreement and rollback, stream/ID assignment, real device execution and transfer/scaling/resource measurements |

External side effects and arbitrary callbacks remain outside the bounded
transaction. A future callback interface requires a concrete publication and
failure contract rather than claiming rollback can undo external effects.


## Output-cost inventory and retained stability failure

The separately frozen `output-cost-contract-v1.md` measures actual public save,
load and archive allocations and serialized sizes over the same 180-profile
matrix, eight warmups and three repeated states. Julia 1.12.7 reports **6,065 pass /
236 fail**; Julia 1.13.2 reports **6,073 pass / 228 fail**. Every failure is the
unchanged identical-byte stability check: save/load/archive fail in 50/120/66
minimum and 44/119/65 latest profiles. Full source/destination/RNG equality,
array identity, archive checksum/state, file sizes and zero compilation checks
pass. The inventory supplies raw costs; it supplies no allocation ceiling or
latency acceptance. No production code or original gate was changed.

`output-owners-v1` passes 642 diagnostic assertions/runtime for D1 Float32, one
empty population, cold zero fields. For each RNG, 64 semantically equal fresh
checkpoint graphs give varying canonical configuration-
validation bytes; all three measurements of each fixed graph are identical.
Xoshiro ranges are 39,488–44,704 B minimum and 40,576–45,792 B latest; MersenneTwister
is 43,280–48,496 B on both. The full validator minus this graph owner is constant
32,704 B minimum / 30,592 B latest across every row. Complete allocation stacks
include the existing IdDict path; all first-use profile records remain retained.
These observations do not establish a violated numerical/public-state contract or
justify rewriting the existing configuration owner.

`output-cost-v1/independent-failure-analysis-v1.json` independently reconstructs
the complete matrix and exact failure counts, retains the gate as **failed**, and
rejects state/compiler/missing-profile mutations. The original strict V1 replayer
also exits 1 as expected; it was not changed to report success. This limitation is
separate from the passing full-step zero-additional-byte and retained-memory gates.

Final frozen source, own patch and complete evidence inventory are retained in
`final-handoff-v1` under the evidence directory above. This is a review handoff,
not final dependency acceptance, main promotion or CI.


## Port onto main (2026-10-26)

Base `cda63f1`. The composition calls the landed `collide_coulomb!` and
`CoulombCollisionWorkspace`; none of the provisional V10 collision code is carried.
Probe scripts and logs are under `/tmp/hpic-sub-cap010/{probes,logs}`.

- `HybridStepper` `step!` now delegates its body to
  `_step_hybrid!(st, ps, dt, NB, source)`. Plain stepping passes `nothing`; the
  composition passes itself as the only source, so the collision runs after the
  push and magnetic advance and before the integer moments, Ohm's law and the
  carried-E recentering. This equals the candidate's step, collision and refresh
  bitwise (permanent canonical-sequence test) without forming E twice. Plain
  stepping is bitwise unchanged against `cda63f1` (scalar and collection,
  isothermal and evolving pressure, NGP/CIC/TSC, D=1–3, Float32/64, variable steps
  with zero steps; identical SHA-1 per scenario on both runtimes).
- Internal calls use the bound, concretely typed species vector. Test files run
  standalone through `test/collisional_fixtures.jl`. Added tests: zero logarithms
  against a separate unwrapped `HybridStepper` with interleaved zero steps; the
  fixed-step contract after stepping; file restart at step 4 of 10 with a
  seed-only negative control; next-kick carried-E feedback; bounded output
  retention.

**Finding 1, `physical-feedback-v1` (833 pass / 32 fail).** Reproduced on the port:
the same 48-profile gate at h=0.01 fails the same 32 Float32 next-step detection
assertions on Julia 1.12.7 and 1.13.2. After one step the composition and the
stale-E variant share x, v, B and RNG; the Boris kick is affine in E, so the next
kick differs by `(h/2)(q/m)(I+R)dE` with `|R-I| <= h*Omega`. The measured
pre-collision difference equals `h(q/m)dE` within that bound in all 48 profiles
(ratio of maxima 0.9993–1.0040 on both runtimes). The signal is 2.1e-6–2.3e-5 in
both precisions: a physical quantity set by the collisional change of the
charge-weighted bulk velocity, not by rounding. A fixed
`128 eps(Float32) max(1,|v|) = 1.53e-5` detection threshold at h=0.01 is therefore
not implied, and every failing row has its measured signal below it. The
expectation was wrong; no implementation defect. The permanent replacement checks
bitwise equality with the refreshed component sequence (a stale-E implementation
fails it at any precision) and the derived kick identity.

**Finding 2, `output-cost-v1` byte stability.** Reproduced: save/load/archive
bytes differ between identical-state samples in 52/138/71 of 180 profiles on
1.12.7 and 46/131/68 on 1.13.2 (largest spread 8.8% / 9.0%); all state, RNG,
array-identity, archive and zero-compilation checks pass. `Base.IdDict` alone shows
the mechanism: inserting 64 fresh sets of 600 mutable keys allocates either
24,576 or 90,144 bytes, while one fixed key set always allocates the same amount
(both runtimes). Its table grows when a probe sequence of address-based identity
hashes overflows; `deepcopy`, Serialization and the canonical checkpoint graph
validation all use `IdDict` on freshly allocated payload objects. The plain
`HybridStepper` species checkpoint on `cda63f1`, with no composition code, also
varies (save/archive in 3/15 of 72 profiles on 1.12.7 and 2/19 on 1.13.2). Byte
identity is not a valid expectation for these operations. The output invariant
that holds is bounded retention: at most 1,856 B of live-heap growth after 64
further save/load/archive cycles (1 MiB budget) on both runtimes; a permanent
test now enforces it.

**Verification on the port (Julia 1.12.7 and 1.13.2, identical counts).**
Collisional suites 1,220 + 181 + 956 assertions; neighbors `test_hybrid` 122,
`test_multispecies_hybrid` 1,014, `test_variable_steps` 277, `test_step_clocks`
676, `test_evolving_pressure` 316, `test_local_coulomb` 1,260, `test_collisions`
518, `test_io` 36,471, Aqua 9, JET 3, all passing. JET package reports equal the
`cda63f1` baseline in the same environments (29 / 27, identical after removing
source locations). The original angular (4,194,512) and weight (655,457)
campaigns pass unchanged on the landed sampler. Evidence-only overlays inject an
ordinary exception and an `InterruptException` after the field advance, the first
pair, the collision, the integer moments, Ohm's law, the recentering, the final
validation and during restore: 616 assertions per runtime confirm the same error
object, complete in-place rollback including the RNG, and identical continuation.
A file restart in a new process reproduces the uninterrupted 12-step digest in 12
of 12 profiles. The composed step allocates exactly as much as the plain step plus
a direct `collide_coulomb!` once the wrapper stores the concrete grid type (an
abstract field had cost a 64–80 B dynamic keyword call); its process CPU time is
1.006–1.082 times that sum (minimum of 21 interleaved samples, heavily loaded
machine; wall-time ratios were not reliable under that load).

Remaining CAP-010 work is unchanged: variable-step collision scheduling,
higher-order splitting, other closures/models and sources, distributed cell
ownership and device RNG.
