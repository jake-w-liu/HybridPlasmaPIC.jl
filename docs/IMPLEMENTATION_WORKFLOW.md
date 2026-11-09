# Parallel implementation workflow

Updated 2026-11-09. The full simulator goal remains active and incomplete.
The [completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md) owns the full
CAP-001–011 scope; the [audit tracker](DEEP_DEBUG_TRACKER.md#current-progress)
owns checked revisions, failures and acceptance. This page owns scheduling,
exclusive file writers and handoffs.

## Working rules

1. Give each actual path one writer, including across worktrees. Root owns
   registrations, exports, dependencies, workflows, current global docs and CI.
2. Complete real evolution, diagnostics, failure behavior, physical phase and
   persistence. A type, primitive or prescribed-field surrogate is not a complete
   self-consistent simulator family.
3. Preserve canonical numerical owners where equations agree. Freeze a small
   dependency interface early; never copy physics to bypass its owner.
4. Keep independent source work moving while GitHub runs. During this resumption
   Julia/parser/compiler/formatter/test/device workloads stay on remote runners.
   Preserve unrelated local edits and inherited processes/CI.
5. Keep every existing numerical, compiler, allocation and cost gate. Push each
   reviewed candidate to a validation branch; fast-forward main only after that
   exact revision passes CI. Record actual billing evidence for an exception.
6. Record source implementation, source review, CI, scientific and resource
   acceptance separately. Real backend claims require the actual device/ranks.

## Verified baseline and current work

Main is `9b4d363806722319c2241e4af04ed1ed1da00119` at this checked snapshot.
[Run 41813949753](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/41813949753)
passes all 26 jobs: minimum/latest Julia across eleven groups, formatting and
real MPI at 2/4/8 ranks. Root retained every raw job log and exact head.
The source stack, compatible implicit Yee EM and wave-1 species-schedule,
implicit-hybrid and slab-gyrokinetic controls are on main. Green CI is separate
from the remaining scientific/resource gates and a fresh source audit.

SpectralOperators 0.1.9 is registered, annotated-tagged and
[released](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.9)
at `76898af`; its candidate/main five-job matrices pass. HPIC already requires
that version. Root completed the previously missing tag/release and reconciled
its publication record.

| Track | Worker / owner | Deliverable | Stage |
| --- | --- | --- | --- |
| Boris certificate premises | Root integrates frozen `post_devin_numeric_audit` repair | Underflow-safe original quotient equality; DD range/mode/domain and outward-distance gates; symmetric zero products | Two canonical files reviewed/adopted; compact remote regressions prepared; executable acceptance pending |
| Dynamic composite ES | Root adopts final five-file handoff | Real completed-boundary refine/coarsen, canonical charge/field reconstruction, mesh-energy exchange, physical steps and complete restart | Complete bounded source/caller/phase/IO review and registrations; essential remote acceptance pending |
| Implicit Yee transverse acceptance | Root adopts frozen handoff | Independent true-Yee vacuum and unequal-species cold-plasma midpoint controls | Both references source-reviewed and registered; exact-head remote acceptance pending |
| Physical-pencil model acceptance | Root adopts frozen lifecycle checks | Actual charged distributed ES/EM/hybrid and nonconstant Maxwell evolution, restore and continuation controls | Complete source review, independent initial density/Maxwell references and minimum/latest2-rank registration; remote execution pending |
| Resident CUDA hybrid | Root holds frozen reconciled source | Complete resident D1/3V isothermal lifecycle, diagnostics and physical persistence | Nine files reconciled against current main; registration and actual NVIDIA acceptance remain open |
| Integration and CI | Root | Source review, shared registrations, current Markdown, commits/main publication and exact-head CI | Active alongside the three independent workers |

## Exclusive file writers

**Boris integration transfers to root:** `src/particles/boris_accuracy.jl` and
`src/particles/boris.jl`. Other arithmetic owners, formulas and tolerances remain
outside that repair. Root owns the existing accuracy test file and remote gates.

**Dynamic ES worker owns only five new paths:**

- `src/meshes/dynamic_composite_adaptation.jl`
- `src/models/dynamic_composite_electrostatic.jl`
- `src/diagnostics/dynamic_composite_electrostatic.jl`
- `src/io/dynamic_composite_electrostatic_checkpoint.jl`
- `docs/DYNAMIC_COMPOSITE_ELECTROSTATIC.md`

The static cloud, Neumann solve, KDK and physical-checkpoint owners remain
canonical. Regrid retains particles and global time, reconstructs the constrained
field on the new mesh and records actual signed mesh-energy exchange. Prior
physical-step history keeps its own detached mesh/epoch. Same-topology adaptation
is a validated no-op. Root owns includes/exports/writer unions/tests/workflows.

**Implicit Yee acceptance transfers to root:** new
`test/test_implicit_yee_transverse.jl`, `docs/IMPLICIT_YEE_EM.md` and canonical
full/focused registrations. Reference matrices use direct cyclic incidence and
original species normalization, without production residual/current/curl helpers.

**Physical-pencil acceptance transfers to root:** new
`test/test_distributed_models_core.jl` and `tools/ci_distributed_fft.jl` in its
isolated worktree. Root owns workflow/docs integration. The completed audit
covers actual model/checkpoint paths in `ext/pencil_*.jl`; any confirmed
production defect requires explicit canonical file ownership before editing.

**CUDA worker reviews the existing nine-file freeze:** `src/models/cuda_hybrid.jl`, seven files under
`ext/cuda_hybrid/` (`types`, `validation`, `kernels`, `fourier`, `evolution`,
`snapshot`, `checkpoint`) and `docs/CUDA_HYBRID.md`. Root retains extension/module,
Project compat and shared IO registrations. The canonical Boris repair is a
shared dependency; no CUDA-specific quotient workaround is added. The repository
Actions API currently lists no self-hosted runner; actual NVIDIA execution remains
open and no hardware-unavailability conclusion is inferred.

## Handoff and verification

A handoff includes complete owned files/patch, pre-edit contract, checked callers,
source/hash scope and unresolved runtime/resource gates. Existing-owner repairs
must retain the intended formulas, errors, ranges, publication and persistence.
Incoming alias checks precede detachment; prepare all replacement state before
publication. Do not reinitialize or reconstruct saved physical phase silently.
Root reads and composes the complete result, keeps essential independent controls
compact and submits the exact revision to GitHub. No source-only pass is called
runtime, scientific or performance acceptance.

## Full implementation queue

| Work | Required completion boundary |
| --- | --- |
| Current certificate repairs and implicit/pencil acceptance | Original numerical contracts, actual public lifecycle and complete continuation; no suppressed findings |
| Resident GPU families | Real device lowering, semantics, failure/restore and measured resources; broader precision/shape/model/backend combinations remain |
| Dynamic AMR and geometry | Current dynamic 1D ES family, then multidimensional/periodic/EM hierarchy, subcycling and full toroidal/variable-J geometry |
| Sources, boundaries and electron physics | General particle/EM boundaries, multiple events, broader transport/anisotropy/inertia and electron-ion exchange |
| Composition | Explicit supported species/collision/source/schedule/geometry/MPI/device/restart combinations |
| Scientific validation and optimization | Independent physical/code comparisons, space/time/particle refinement and fixed resource targets measured remotely |

The plan's **23/79 checked (29.1%)** count remains unweighted and does not measure
full simulator completion. All eleven capability areas remain open. No check is
raised because code exists, jobs queue or an unrelated primitive passes.
