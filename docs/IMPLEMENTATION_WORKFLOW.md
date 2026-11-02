# Parallel implementation workflow

Updated 2026-11-02. The full simulator goal remains active and incomplete.
The [completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) owns the full
CAP-001–011 backlog; the [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress)
owns checked commits, failures and acceptance. This page owns current scheduling,
file ownership and handoffs.

## Working rules

1. Keep three independent implementation tracks active. Root is the fourth
   worker and owns integration, shared registration files and remote CI.
2. Give each repository path one writer, including across isolated worktrees.
   A worktree prevents accidental interference; it does not resolve conflicting
   implementations of the same owner. Request an ownership transfer before
   editing another track's file.
3. Build complete bounded solver families: real evolution, diagnostics, failure
   handling, physical phase and restart. Do not count a type, exported function,
   prescribed-field surrogate or placeholder as a completed solver.
4. Freeze the small actual interface needed by a dependent track early. Continue
   independent work while that interface is prepared; do not duplicate its owner.
5. Keep source implementation moving while GitHub jobs queue. Preserve failures
   and investigate them in parallel. No local Julia/test/compiler/formatter,
   numerical-proxy, device, build or profiling workloads.
6. Submit reviewed candidates through GitHub, retaining every existing numerical,
   compiler, allocation and cost gate. Root pushes each completed fix or
   implementation to `main` and checks its exact-revision CI. A failed candidate
   stays unaccepted; a billing exception requires actual GitHub evidence.

## Current assignments

| Track | Worker | Complete deliverable | Current stage | Dependency |
| --- | --- | --- | --- | --- |
| Hall A | `cpu_model_factor_admission` | `OpenHallHybrid` core, mixed CIC, two-face sources/IDs/bounded RNG, projected field controller, exchange/energy diagnostics and guide | Complete source frozen and read; root integration/remote essentials prepared | Actual v5 source/report/commit interface and mixed operators |
| Hall B | `yee_longitudinal_coupled_implementation` | Complete Hall physical/source/RNG checkpoint, archive, restore and guide | Revised complete IO frozen and read; root integrated against final A v5 | Actual A state/configuration/source/report/commit hooks |
| Gyrokinetics | Root integration of `self_consistent_gyrokinetic_design` handoff | `GyrocentreSet`, self-consistent full-f slab field/distribution/RK4, diagnostics, complete checkpoint/archive and guide | Complete bounded source/callers reviewed; integrated essentials prepared for GitHub | Registered SpectralOperators 0.1.8; runtime/scientific/resource acceptance remains open |
| Curvilinear next slice | `self_consistent_gyrokinetic_design` | A complete geometry/coupling/field/phase/persistence contract before production edits | Mapping actual owners and choosing bounded model | Coordinate/component/volume conventions and compatible operators; toroidal axis/walls remain explicit |
| Implicit EM next slice | `yee_longitudinal_coupled_implementation` | A complete joint collection/field residual, nonlinear solve and physical continuation contract | Reading canonical current correction, Boris, Maxwell and nonlinear owners | EM continuity/work residual review; Hall persistence amendments take priority |
| Static AMR next slice | `cpu_model_factor_admission` | A real static composite ES hierarchy/field/coupling/phase/persistence contract | Reading existing hierarchy, cloud and field owners; no production edits | Global field/interface work and particle crossing/self-force contract before dynamic refinement |
| Integration and CI | Root | Compose shared registrations/docs, source review, dependency publication, commits/pushes and exact-revision CI closure | Active | Complete candidate handoffs and actual remote results |

Hall and gyrokinetic implementations start from published `531f5ba` in separate
worktrees. Current repair candidates and historical primary work remain isolated.
No unfinished family is represented as runtime accepted.
Hall A/B are composed above the slab candidate `08e1bf2`, retaining exact frozen
author receipts and legacy inverses. The slab integration candidate builds on
the separate `bb7da93` correctness
candidate; neither candidate has replaced `main`. The next three contracts do not
authorize edits to current-family files or shared registrations.

## Exclusive file ownership

**Hall A** owns these existing canonical seams:

- `src/coupling/shapes.jl`
- `src/coupling/moments.jl`
- `src/coupling/gather.jl`
- `src/fields/mixed_shock_operators.jl`
- `src/electrons/isothermal_energy.jl`
- `src/boundaries/particle_reservoir.jl`

Hall A also owns new `src/models/open_hall_hybrid.jl`,
`src/coupling/mixed_cic.jl`, `src/boundaries/open_hall_reservoir.jl`,
`src/diagnostics/open_hall_energy.jl` and `docs/OPEN_HALL_HYBRID.md`.

**Hall B** owns only new `src/io/open_hall_checkpoint.jl` and
`docs/OPEN_HALL_HYBRID_CHECKPOINT.md`.

**Gyrokinetics** owns only new `src/coupling/tsc_derivatives.jl`,
`src/models/slab_gyrokinetic_pic.jl`, `src/io/slab_gyrokinetic_checkpoint.jl`
and `docs/SLAB_GYROKINETIC.md`. Its TSC derivative reuses the existing stencil
base/values; it does not edit Hall A's shape file or duplicate TSC weights.

**Root** owns `src/HybridPlasmaPIC.jl`, shared checkpoint/archive registrations,
`Project.toml`, workflows, README, CHANGELOG, global trackers/documentation map,
and current CI repair/observation files. Workers supply requested registration
hunks rather than editing these files. SpectralOperators dependency fixes and
release/registration also remain root-owned.

## Interface and handoff

Hall A freezes actual configuration/binding/storage, source/phase/report,
replacement/reserve and final commit hooks for B. B saves the full physical
state and complete owned RNG, validates the original incoming graph before
copying, prepares replacements privately, and restores stored state verbatim.
A seed or reinitialization is not a replacement for continuation state.

Each worker hands root the complete patch, owned path list, model guide,
source/caller review and any unresolved issues. Existing-owner extractions must
show that legacy arithmetic, checks and publication survive the extraction.
Keep evidence proportional to the change; do not create a new tracking framework.
Root reviews canonical seams once, composes registrations, checks local links,
and submits the exact candidate to GitHub.

## Implementation queue after current handoffs

Assign the next ready item when a worker completes its current deliverable.
Keep ownership disjoint and use the existing completeness plan for detailed
model contracts and all remaining combinations.

| Wave | Parallel work | Required boundary |
| --- | --- | --- |
| Current | Complete Hall A/B; self-consistent slab gyrokinetics; known CI/dependency repairs | Complete physical lifecycle and persistence; no shared-file writers |
| Next | First resident CUDA hybrid; curvilinear coupled solver; implicit EM model and residual | CUDA needs reviewed shared scalar seams and real device acceptance; geometry and EM need their own field/coupling/phase contracts |
| Following | Static composite-grid ES hierarchy then dynamic AMR; broader electron transport/anisotropy/inertia; broader open EM particle/source boundaries | Composite field/particle ownership, term-specific energy balances and boundary charge/current accounting |
| Composition | Extend accepted families across species, collisions/sources, schedules, MPI/device storage and restart | Explicit supported-combination matrix; no automatic claim from component availability |
| Continuous | Independent physical validation, latest Julia, failure/restart checks and measured optimization | GitHub or suitable remote hardware; retain failed evidence and fixed acceptance limits |

Large new families may use separate core and persistence workers after their
interface freezes, as Hall does now. Do not start more active tracks than the
four available workers can own and review.

## Progress stages

Record these separately for every bounded slice:

| Stage | Meaning |
| --- | --- |
| Planned | Equations, public contract, dependencies and file owner are specified |
| Implementing | Real source is being written; missing behavior remains explicit |
| Implemented | Complete admitted evolution, diagnostics, failure handling and persistence are present |
| Reviewed | Canonical owners/callers, edge cases and full patch have been inspected |
| CI accepted | Exact-revision remote checks pass; record run, Julia versions and scope |
| Scientifically/resource accepted | Required physical comparisons, refinement, real ranks/devices and fixed cost gates pass |

The existing **14/79 (17.7%)** number counts checked roadmap milestones; it is
not a source-code completion percentage. Do not increase it for drafted code,
line counts or queued CI. Update stage/owner/failures here and in the audit
tracker at each handoff, integration and acceptance change.
