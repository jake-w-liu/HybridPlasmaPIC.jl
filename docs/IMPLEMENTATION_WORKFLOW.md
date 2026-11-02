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

The periodic-shear candidate `c53ef77`, branch
`implement/sheared-periodic-hybrid-20261102`, contains Hall and slab. Its
[full](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41455172923) and
[essential](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41455205769)
checks are submitted and still unaccepted. The implicit EM candidate `d665c47` builds on it; its
[full](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41457743663) and
[essential](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41457768797)
checks are submitted. Separate hybrid-clock fix `ae23c71` has
[full](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41458375525) and
[essential](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41458455218)
checks submitted. Static composite ES composes on that exact source. None is
runtime accepted; candidates remain off main while their gates are pending.

Hall `cd78339` contains slab `08e1bf2` and correctness candidate `bb7da93`.
Main remains `531f5ba`; none of those candidates is runtime accepted. Frozen
author worktrees remain intact. Root owns their integration and any subsequent
repairs; new workers do not edit those family files.

## Active assignments

| Track | Worker | Complete deliverable | Current stage |
| --- | --- | --- | --- |
| Static composite-grid electrostatic PIC | Root integration of frozen `cpu_model_factor_admission` source | Nested hierarchy, global Neumann field solve, physical-cloud crossings, diagnostics and complete physical persistence | Seven-file v2 source reviewed/integrated; phase/preparation findings repaired; essential GitHub checks prepared |
| First resident CUDA hybrid | `self_consistent_gyrokinetic_design` | Complete device field/particle lifecycle and physical persistence with actual CUDA/cuFFT APIs | Complete v2 contract and exact callable seams approved; nine new family files implementing |
| Shared CUDA scalar dependencies | `cpu_model_factor_admission` | Canonical arithmetic/coupling/scale/native Ohm/RK4 extractions with unchanged default CPU behavior | Eleven disjoint paths implementing; root owns upstream SO and registrations |
| Compatible implicit EM current/gather | `yee_longitudinal_coupled_implementation` | Paired orbit-averaged current/gather, midpoint field/particle solve and physical persistence | Complete contract and 29 source pins frozen for root review; no production edits |
| Integration and CI | Root | Shared registrations/docs, canonical source review, dependency publication, commits/pushes and exact-revision CI closure | Active alongside the independent tracks; existing clock-admission audit frozen for follow-up repair |

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

Completed Hall, slab, periodic-shear, synchronous implicit EM and static composite
ES author patches are frozen. Their owners transfer to root for integration and
later repairs. Static ES includes the byte-identical `_yl_source_scalar` move
from `models/yee_longitudinal_pic.jl` to `utils/validation.jl`; canonical cloud,
force and plate owners remain in place. Broader clock findings retain a separate
frozen audit; do not edit those owners without a new root assignment.

**CUDA family** owns only these nine new paths:

- `src/models/cuda_hybrid.jl`
- `ext/cuda_hybrid/types.jl`
- `ext/cuda_hybrid/validation.jl`
- `ext/cuda_hybrid/kernels.jl`
- `ext/cuda_hybrid/fourier.jl`
- `ext/cuda_hybrid/evolution.jl`
- `ext/cuda_hybrid/snapshot.jl`
- `ext/cuda_hybrid/checkpoint.jl`
- `docs/CUDA_HYBRID.md`

**Shared scalar dependencies** owns new `src/utils/scalar_ops.jl` and these ten
existing paths, including across worktrees:

- `src/particles/boris.jl`
- `src/particles/kinematics_range.jl`
- `src/coupling/shapes.jl`
- `src/coupling/deposit.jl`
- `src/coupling/gather.jl`
- `src/coupling/moments.jl`
- `src/models/hybrid_pic.jl`
- `src/electrons/closures.jl`
- `src/electrons/ohms_law.jl`
- `src/integrators/extrapolated_leapfrog.jl`

Its default expansion must retain current CPU arithmetic, generic branches,
validation, exact cold behavior and read/write/reduction order. CUDA consumes
frozen signatures and record fields; neither worker duplicates a physics owner.
Actual device lowering, semantics and performance still need remote acceptance.

**Compatible implicit EM** currently owns design evidence only. Root reviews
its complete method/geometry/source/phase/IO contract before production ownership.

**Root** owns `src/HybridPlasmaPIC.jl`, existing extension includes, shared
checkpoint/archive registrations, `Project.toml`, workflows, focused essentials,
README, CHANGELOG, global trackers/documentation map and CI repair/observation
files. Workers supply registration hunks. SpectralOperators canonical mode-row
and D1 curl-descriptor changes, all upstream CI/release/registration and the
actual accepted consumer compat floor remain root-owned. Every path has one
active writer; isolated worktrees do not permit competing edits.

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
| Current | Static composite ES integration; resident CUDA family and shared scalar implementation; compatible implicit EM contract review; known CI/dependency repairs | Complete physical lifecycle and persistence; no shared-file writers |
| Next | Compatible implicit EM implementation after review; broader geometry; completed-phase admission repairs | CUDA needs reviewed shared scalar seams and real device acceptance; every new method keeps its own field/coupling/phase contract |
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
