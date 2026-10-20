# Particle-growth audit

Status: append defects confirmed on both supported runtimes; repair planned.
Ionization mutation also reproduced; the isolated candidate is under verification.

Current `append_particles!` checks species charge and mass, then appends
position, velocity, weight, id and tag arrays sequentially
(`src/particles/migration.jl:36`). Its implementation comment requires resizable
fields, but the public `ParticleSet` constructor accepts matching one-based
`AbstractVector` storage, including views. `ionize_mcc!` cools selected primaries
before appending the electron and ion batches (`collisions.jl:1177`).

## Candidates

- A nonresizable late destination field may throw after earlier fields grew,
  leaving inconsistent particle-array lengths.
- Internally aliased destination fields may be appended more than once and
  corrupt the structure-of-arrays length invariant.
- Distinct source and destination populations sharing components may change
  a later source field's length before it is appended. Check rejection before
  either population changes, while retaining self-append separately.
- An ion destination that cannot grow may fail after primary cooling and
  electron creation. Determine whether both species and id counters remain
  unchanged on that input error.

**VERIFIED — GROWTH-001:** the public probe at
`/tmp/hybrid-particle-append-review-20261020` reproduces all three append defects
on Julia 1.12.7 and 1.13.2 against source byte-identical to current main.
For a 3D two-particle destination and one appended particle, a nonresizable tag
leaves x/v/weight/id lengths 3 and tag length 2. A destination sharing x₁ and
weight ends with lengths 4 versus 3, without throwing. Distinct populations
sharing destination x₁ with source weight end with destination lengths 4 versus
weight length 6, also without throwing; source weight changes too. The same
invariants fail for both precisions in 1D–3D (36 failed assertions per runtime).

The initial ionization fixture instead raises an ambiguous `randn` dispatch
before mutation; its unchanged-state assertions therefore establish no
ionization safety result. V2 preserves the complete original assertions,
resolves only that RNG method ambiguity and adds direct RNG plus successful
birth controls. Its frozen source and baseline logs are retained separately at
`/tmp/hybrid-particle-append-review-v2-20261020`.

**VERIFIED — GROWTH-001 ionization path:** V2 completes on both Julia versions
with 172 passes and 48 failed invariants each. All RNG/success controls pass.
The ionization call now reaches `resize!(tag_view, …)` and fails after creating
the electron and partially creating the ion: electron arrays have length 2,
ion x/v/weight/id have length 1, and ion tag remains length 0. Both species
change; id counters remain unchanged. All six precision/dimension fixtures
reproduce this. `summary.json` retains commands and complete log hashes.
The source mutation order at `collisions.jl:1167–1178` accounts for the observed
partial update; no RNG error is used as evidence for this finding.

## Inspected callers

`inject_face_1d!` publishes its accumulator/id counters after append, so the
canonical append preflight can protect its particle state before growth.
Ionization needs both destinations checked before cooling or the first append.
Serial `migrate_particles!` wraps positions and calls `_replace_particles!`
before appending arrivals. Real `mpi_migrate_particles!` also replaces local
storage before collective transport. Their replacement contracts require a
separate check; an append-only repair cannot establish their rejection safety.

The checkpoint path already checks component aliases and `applicable(resize!, …)`
before its restore writes (`src/io/checkpoint.jl:894`). Inspect whether that
storage preflight can have a shared particle-storage owner, preserving its
existing exception text and supported independent views. Do not introduce a
particles-to-I/O dependency or duplicate an existing validator.

## Pre-edit repair contract

The canonical owner is particle storage, with append remaining in migration.
Move the existing alias predicates from collisions to `particle_set.jl` without
changing default behavior. Growth additionally checks backing storage and
mutable view indices, because resizing can affect disjoint or empty wrappers;
checkpoint restoration must keep its existing allowance for independent
matching-size views and its existing exception text. Share its resize-support
loop in particle storage rather than duplicating it or depending on I/O.

For nonempty append, validate consistent one-based component axes, destination
capacity without integer wrap, all destination resize support, and independent
destination components before any append. Reject shared backing between
distinct populations before mutation; preserve self-append of independent
components, internally aliased read-only sources, and empty-source no-ops.
Retain the original amortized in-place append loop: no full destination copy,
no per-particle validation or temporary payload for ordinary append.

If the corrected ionization probe confirms mutation, preflight both destination
populations after finding a nonempty selected batch and before cooling any
primary or appending either species. Retain valid no-birth calls with view
storage. Counters publish only after both appends as before. The contract covers
these detectable input errors, not rollback from arbitrary user-defined
mutator failures or allocation failure.

Verification before promotion: reproduce the original failures on the repaired
public APIs; cover every destination component, empty/nonempty aliases,
backing/index aliases, malformed lengths, custom resizable storage, independent
and self-append, read-only sources, ionization/no-event controls and injection
rollback. Run neighboring particle/collision/migration/injection/checkpoint
regressions on Julia 1.12.7 and 1.13.2, full compiler comparisons and concrete
calls, Aqua and formatting. Require zero allocation for warmed ordinary
preflight and unchanged ordinary append allocation with retained capacity;
measure append cost without claiming a timing ceiling from noisy ratios.
Push the verified repair to main with its tracker evidence, then require all
six exact-SHA GitHub jobs. Migration replacement remains a separate open audit.

## Candidate V1

`/tmp/hybrid-particle-growth-v1-20261020/source` is a clean archive of `0da4e09`
with the planned shared storage preflight, append validation and both-species
ionization check. The held radiation work is excluded. Default alias behavior
and checkpoint error messages are retained; only growth requests the backing-
storage predicate. On Julia 1.12.7, the original probe passes all 220 assertions.
The initial
allocation check fails all six configurations at 384/560/640 bytes for 1D/2D/3D.
A concrete helper-by-helper measurement places all those bytes in resize
support; axes, component-alias and cross-population-alias checks each allocate
zero bytes. `allocation-diagnosis112.log` retains these measurements.

The expanded V1 suite records 553 passes and two setup errors: mixed vector/view
position or velocity tuples are rejected by the convenience constructor before
the intended operation. Those errors leave most of the new append matrix
unexecuted; they do not establish a production defect or successful coverage.

V2 separates homogeneous position/velocity resize checks from weight/id/tag
checks, preserving predicates and error strings. The isolated prototype first
passes six zero-allocation assertions on minimum Julia. V2 passes all 220
original checks and all six zero-allocation checks on each
runtime. Its 36-profile resource matrix passes 3,060 assertions per runtime:
ordinary preflight and retained-capacity appends allocate zero bytes; complete
outputs and source storage match the baseline. Twelve alternating measurements
per profile retain raw times. Across these small append cases, observed median
ratios are 3.75–8.12 on minimum Julia and 3.00–8.08 on latest Julia; largest
absolute medians are 6.44/3.56 microseconds. These are cost observations, not an
optimization claim or a new elapsed-time acceptance threshold.
V2 then retains the same two known fixture errors (553 passes); later gates do
not run on that candidate. V3 keeps V2 production code unchanged and fixes
only the mixed-storage test fixture by using the existing explicit parameterized
`ParticleSet` constructor; it adds a shape assertion and retains every original
assertion. V3 passes 1,784 focused checks on each runtime, 347 formatting checks and
108 source-migration checks. The latter compare all moved helper bodies with
the original default predicates, preserve the checkpoint resize error and
preflight behavior, and compare every untouched production syntax tree.
Both runtimes also pass 226 smoke checks and 3,060 retained-capacity append
checks. Neighboring integration and compiler gates remain running or queued. V3 production is byte-identical to V2 across all 88 source files.
Exactly six source/test paths differ from the pinned main baseline. Nothing
is promoted or accepted yet.


A further checkpoint resource probe initially assumes that the entire old
checkpoint validator allocates zero bytes. That assumption fails in all 54
profiles on each runtime, producing 108 failed assertions per run; the logs at
`/tmp/hybrid-particle-growth-checkpoint-resources-20261020` are retained.

Direct measurements of the actual old and new package owners resolve that
question at `/tmp/hybrid-particle-growth-checkpoint-resources-v2-20261020`.
Each of four commands passes 108 behavior/ownership checks. For every one of
54 matched profiles per runtime, candidate allocation is lower: concrete
vectors go from 384–640 bytes to zero, matching-length views from 1,776–2,864
to zero, and abstract vector fields from 256–688 to 192–464. This is no
checkpoint allocation regression; the first probe's zero-baseline assumption
was incorrect. `summary.json` records all profiles and comparisons.

Serial migration's separate storage hypotheses are being tested at
`/tmp/hybrid-particle-migration-storage-review-20261020` against the pinned main
owner. The probe checks both early and late view failures, successful transport
controls, and preservation of weights when positions share backing storage.
Both runtimes now complete with 36 control passes and 24 integrity failures.
No invalid-input baseline MPI collective is launched. The confirmed GROWTH-002
finding has its own [repair plan](MIGRATION_STORAGE_PLAN.md) and remains outside
the append repair.
