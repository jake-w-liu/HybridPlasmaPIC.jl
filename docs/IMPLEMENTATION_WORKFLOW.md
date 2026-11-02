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

## Current source candidates

| Family | Candidate | Source stage | Remote acceptance |
| --- | --- | --- | --- |
| Self-consistent full-f slab gyrokinetics | `08e1bf2`, branch `implement/slab-gyrokinetic-20261102` | Complete core, coupling, diagnostics and physical checkpoint/archive reviewed and integrated | [Full CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41446645551) and [essentials](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41446663047) submitted; no accepted result yet |
| Synchronous open Hall numerical reservoirs | `cd78339`, branch `implement/open-hall-hybrid-20261102` | Complete Hall A/B source composed, including last/cumulative source records, full RNG, diagnostics and physical restore | [Full CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41448552638) and [essentials](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41448552969) submitted; no accepted result yet |

Hall `cd78339` contains slab `08e1bf2` and correctness candidate `bb7da93`.
Main remains `531f5ba`; none of those candidates is runtime accepted. Frozen
author worktrees remain intact. Root owns their integration and any subsequent
repairs; new workers do not edit those family files.

## Active assignments

| Track | Worker | Complete deliverable | Current stage |
| --- | --- | --- | --- |
| Periodic curvilinear hybrid | Root integration of frozen `self_consistent_gyrokinetic_design` source | `PeriodicShearGrid` and `ShearedPeriodicHybridPIC`, physical field/component transforms, collective kinetic-ion/field RK4, diagnostics, quiet loader and full physical persistence | Complete seven-file v2 source reviewed/integrated; receipt invariant repaired; remote acceptance pending |
| Implicit electromagnetic PIC | `yee_longitudinal_coupled_implementation` | `ImplicitEMPICStepper`, joint velocity/E/B residual, canonical continuity current, shared bounded nonlinear solver, signed work diagnostics and full physical persistence | Complete contract approved; core/driver/IO/guide drafted from `cd78339`, final source/inverse review in progress |
| Static composite-grid electrostatic PIC | `cpu_model_factor_admission` | A real static hierarchy, one composite field solve, compatible finite-cloud coupling, particle crossings, diagnostics and full physical persistence | Complete 373-line contract approved and 33 source pins checked; implementation active from `cd78339` |
| First resident CUDA hybrid | `self_consistent_gyrokinetic_design` | Complete device field/particle lifecycle and physical persistence with actual CUDA/CUFFT APIs | Contract/source/seam preparation; no production edits |
| Integration and CI | Root | Shared registrations/docs, canonical source review, dependency publication, commits/pushes and exact-revision CI closure | Active alongside all three source tracks |

The curvilinear contract keeps logical positions and Cartesian velocities,
primary contravariant magnetic flux and a compatible covariant electric curl.
Its canonical Ohm source-risk flag must be false before differentiation; the
continuous Lorentz rate retains original q/m and complete field factors. This
periodic nonorthogonal map does not supply toroidal axis, wall, equilibrium or
variable-J support. Those remain explicit CAP-006A work.

The implicit EM contract solves velocities and both fields together and
recomputes charge/current inside every residual and Jv probe. It exposes actual
Gauss/divB/continuity defects and four signed work terms, including the spatial
work of current correction. It does not inherit an exact energy theorem from a
different current/gather scheme. The first scope is nonrelativistic all-mobile
CPU CIC, exact-neutral collections and odd periodic grids; broader conservative,
relativistic, source and backend families remain open.

## Exclusive file ownership

**Curvilinear** owns only these seven new paths:

- `src/meshes/periodic_shear_mesh.jl`
- `src/fields/sheared_hybrid_operators.jl`
- `src/coupling/sheared_hybrid_coupling.jl`
- `src/models/sheared_periodic_hybrid_pic.jl`
- `src/initial_conditions/sheared_hybrid.jl`
- `src/io/sheared_hybrid_checkpoint.jl`
- `docs/SHEARED_PERIODIC_HYBRID.md`

It reuses the existing complete periodic gather, exact contribution iterator
and node-index policy from the Hall candidate without editing those owners.

**Implicit EM** owns three existing canonical seams:
`src/models/full_pic.jl`, `src/models/empic_collection.jl` and
`src/integrators/implicit_hybrid.jl`. It also owns only these five new paths:

- `src/integrators/implicit_em.jl`
- `src/integrators/implicit_pic_nonlinear.jl`
- `src/diagnostics/implicit_em.jl`
- `src/io/implicit_em_checkpoint.jl`
- `docs/IMPLICIT_EM.md`

Existing-owner extraction must recover each complete old file by literal
expansion, preserving arithmetic, defaults, validation and try/finally order.
The nonlinear driver serves these two actual implicit families; it is not a
new public solver framework.

**Static AMR** owns only these five new paths:

- `src/meshes/static_amr_hierarchy.jl`
- `src/coupling/composite_cloud.jl`
- `src/models/static_composite_electrostatic.jl`
- `src/io/static_composite_electrostatic_checkpoint.jl`
- `docs/STATIC_COMPOSITE_ELECTROSTATIC.md`

Its only existing-owner change is the byte-identical `_yl_source_scalar` move
from `models/yee_longitudinal_pic.jl` to `utils/validation.jl`; that fragment is
already frozen and released to root. Force/plate/cloud owners remain in place.
The CUDA worker currently owns design evidence only; root reviews any proposed
shared scalar fragment before implementation ownership is assigned.

**Root** owns `src/HybridPlasmaPIC.jl`, shared checkpoint/archive registrations,
`Project.toml`, workflows, README, CHANGELOG, global trackers/documentation map,
and current CI repair/observation files. Workers supply registration hunks.
SpectralOperators fixes and release/registration remain root-owned. Every path
has one active writer, including across worktrees.

## Interface and handoff

The completed Hall split demonstrates the handoff boundary: its actual
configuration/binding/storage, source/phase/report, replacement/reserve and final
commit hooks were frozen for the independent IO owner. Every new family must
save its full physical phase, validate incoming original roles before detachment,
prepare replacements privately and restore stored state verbatim. A seed or
reinitialization is not a replacement for stochastic continuation state.

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
| Current | Curvilinear integration; implicit EM and static composite ES implementation; resident CUDA contract; known CI/dependency repairs | Complete physical lifecycle and persistence; no shared-file writers |
| Next | First resident CUDA hybrid after contract review; broader geometry and implicit-method slices | CUDA needs reviewed shared scalar seams and real device acceptance; every new method keeps its own field/coupling/phase contract |
| Following | Dynamic AMR and multidimensional hierarchy; broader electron transport/anisotropy/inertia; broader open EM particle/source boundaries | Composite field/particle ownership, term-specific energy balances and boundary charge/current accounting |
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
