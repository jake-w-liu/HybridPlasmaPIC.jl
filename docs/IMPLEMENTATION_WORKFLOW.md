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
| Boris certificate premises | Root integrates frozen `post_devin_numeric_audit` repair | Underflow-safe original quotient equality; DD range/mode/domain and outward-distance gates; symmetric zero products | Both focused jobs pass 6,142 assertions each at `a15be1c`; formatter correction passes at `62f985f`; full integrated/compiler/resource acceptance open |
| Dynamic composite ES | Root adopts final five-file handoff | Real completed-boundary refine/coarsen, canonical charge/field reconstruction, mesh-energy exchange, physical steps and complete restart | Complete source/registration review; minimum/latest focused jobs pass 136 each at `16321a9`; full integrated/scientific/resource acceptance open |
| Implicit Yee transverse acceptance | Root adopts frozen handoff | Independent true-Yee vacuum and unequal-species cold-plasma midpoint controls | Both references pass minimum/latest focused jobs with 366 assertions each at `62f985f`; full integrated/scientific acceptance open |
| Physical-pencil model acceptance | Root adopts frozen lifecycle checks | Actual charged distributed ES/EM/hybrid and nonconstant Maxwell evolution, restore and continuation controls | Charged latest controls passed before a Maxwell empty-population helper failure; helper corrected in CPU candidate `b8a5e7d`, complete remote rerun pending |
| Resident CUDA hybrid | Root integrates reviewed fourteen-file handoff | Complete resident D1/3V isothermal lifecycle, diagnostics and physical persistence | Latest host job at `38d1cd1` passes 38 checks on Julia 1.13.3; minimum host, full integrated and actual NVIDIA acceptance remain open |
| Toroidal annulus | Root integrated final handoff | Full coupled logical particles, metric magnetic field, adiabatic pressure, contact work and complete restart | Reviewed source/registration pushed at `7ea6ae8`; full and focused remote acceptance pending |
| Implicit evolving pressure | Root integrated final handoff | Canonical joint velocity/B/pe solve, third residual gate, signed physical ledger and schema2 persistence | Reviewed source/registration pushed at `5703c95`, including pinned old-schema raw/archive continuation; exact-head remote acceptance pending |
| Multidimensional adaptive ES | `dynamic_amr_es_implementation` | Genuine coupled D2/3 hierarchy, particle/field adaptation, work and restart | Complete grounded-Q1/tensor-cloud/KDK/regrid/history contract reviewed; eight source paths and compact essentials approved |
| Nondividing species schedules | `annulus_loader_guide` | Canonical EMPIC event/current/phase/schedule persistence | Persistent deadline/current/clock/schema3 contract reviewed; two source owners, two guides and compact essentials approved |
| Field-aligned electron conduction | Root integrates frozen `post_devin_numeric_audit` handoff | Varying-B tensor transport in canonical B/pressure lifecycle and persistence | Eight-owner source and essential file reviewed/adopted on `5703c95`; module/guide/registration wired; all remote acceptance pending |
| Persistent schedule IO | `post_devin_numeric_audit`, reconnaissance until interface transfer | Complete schema3 pending calendar/controller/current/phase persistence | Reads approved contract; edits require exact frozen model seam and exclusive IO-path transfer from schedule owner |
| Integration and CI | Root | Source review, shared registrations, current Markdown, commits/main publication and exact-head CI | Active alongside three parallel implementation tracks |

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

**CUDA integration transfers to root:** `src/models/cuda_hybrid.jl`, seven files under
`ext/cuda_hybrid/` (`types`, `validation`, `kernels`, `fourier`, `evolution`,
`snapshot`, `checkpoint`) and `docs/CUDA_HYBRID.md`. Root retains extension/module,
Project compat and shared IO registrations. The canonical Boris repair is a
shared dependency; no CUDA-specific quotient workaround is added. The repository
Actions API currently lists no self-hosted runner; actual NVIDIA execution remains
open and no hardware-unavailability conclusion is inferred.

**Toroidal annulus final handoff transfers eight paths to root:**
`src/meshes/toroidal_annulus_mesh.jl`,
`src/coupling/toroidal_annulus_hybrid_coupling.jl`,
`src/fields/toroidal_annulus_hybrid_operators.jl`,
`src/models/toroidal_annulus_hybrid_pic.jl`,
`src/diagnostics/toroidal_annulus_hybrid.jl`,
`src/io/toroidal_annulus_hybrid_checkpoint.jl`,
`src/initial_conditions/toroidal_annulus_hybrid.jl` and
`docs/TOROIDAL_ANNULUS_HYBRID.md`. Root owns the canonical scalar map extraction,
includes/exports/shared writer unions, essential remote controls and global docs.
The former parent/child writers have frozen their source; root now owns the integrated eight.

**Implicit-pressure final handoff transfers five existing paths to root:**
`src/integrators/implicit_hybrid.jl`,
`src/integrators/implicit_pic_nonlinear.jl` (scale dispatch only),
`src/electrons/energy_equation.jl`,
`src/integrators/extrapolated_leapfrog.jl` (literal shared-pressure extraction) and
`src/io/implicit_hybrid_checkpoint.jl`, plus new
`test/test_implicit_hybrid_pressure.jl`. Root adopted the exact frozen V2 bytes.
The shared nonlinear algorithm and isothermal report/schema1 remain; root owns
registrations, workflows, global docs and closure/model docstrings. No overlapping
writer is assigned.

**Approved next implementations:** the former annulus parent owns eight new
multidimensional adaptive-ES source/guide paths and one essential file. Its child
owns the two existing EMPIC collection/schedule source/IO files, their two species
guides and one new essential file. The former pressure worker owns the eight
field-aligned-conduction production paths and one essential file from its frozen
contract, rebased onto `5703c95` before edits. These include shared closure,
pressure, explicit-stage, diagnostic, collisional and generic checkpoint owners.
The conduction worker's final eight-owner handoff has transferred to root, so
root again owns the generic checkpoint file. The schedule owner retains model,
guides and essentials; its IO path transfers to the independent IO worker only
after an exact interface handoff. Module, exports, workflows and global docs
remain root-owned. No two workers edit the same canonical file. Accepted
source/runtime/scientific/resource gates remain separate.

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
