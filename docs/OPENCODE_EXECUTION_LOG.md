# OpenCode execution log

Started 2026-11-09 on the isolated checkout
`/Users/jake/EMPIRE/codes/PlasmaWorkspace/HybridPlasmaPIC-electron-ion-overlap-integration-20261109`
of `origin/implement/cgl-integration-20261109`. This log records operational
facts for the P00–P43 handoff. It does not replace the 79-row acceptance
tracker or claim that source presence is executed, scientific or resource
acceptance.

## P00 — takeover inventory and frozen evidence

Verified on 2026-11-09 before any source edit:

| Check | Result |
|---|---|
| Handoff packet | `python3 handoff/2026-11-09/verify_handoff.py` returned `PASS (metadata only; no Julia execution)`: 41 packet files, 79 milestones, 23 legacy checked, 56 open, 44 work packages, acyclic dependency graph, 752 checked document links. |
| Source baseline | `2ea4d79` is a descendant of `5aa628513fb854c96b2b387780edef1ea1faa66b`; `src`, `ext`, `test`, `Project.toml` and `.github/workflows` matched `CODE_TREE_BASELINE.txt` at takeover (`source_equals_baseline: true`). |
| Primary checkout | `/Users/jake/EMPIRE/codes/PlasmaWorkspace/HybridPlasmaPIC.jl` remains at `563888f` with its pre-existing dirty radiation/forcing files. All seven `historical_primary_preserved_sha256` entries were recomputed from the untouched files and matched. |
| Compiler ledger | The manifest hash `6d30ec354917c839f80af386399dc58adbe9103d2ac8f38008ade447dc22d59d` is recorded as the protected identity. No file in this isolated checkout matches that hash; the frozen ledger is an external historical artifact. `test/test_jet.jl`'s cap, dependency filter, `load_run` probes and concrete-call block were left unchanged. |
| Dirty-work preservation | No existing worktree, branch, tracked modification, untracked file, inherited process or CI run was reset, stashed, deleted, cancelled or overwritten. The verified packet artifacts are read-only evidence. |
| Execution boundary | This checkout performs no local Julia, parser, compiler, formatter, test, numerical-proxy, profiler, MPI or device computation. Verification uses GitHub Actions and read-only source/log inspection. |
| Integration lead | OpenCode is the single integration lead for module/export/dependency/workflow/shared-IO/doc changes. Lane paths are assigned exclusively in the execution order below; no two active writers share a path. |
| CI queue state | The handoff push created run `42593522480`. At takeover all 29 jobs were queued; older runs had been waiting since 05:36 UTC. This is an observed queue, not evidence of a billing exception. |

### Lane assignment

| Lane | Packages | Current writer boundary |
|---|---|---|
| CI/integration | P00–P01, P42–P43 | OpenCode owns `src/HybridPlasmaPIC.jl`, `Project.toml`, `.github/workflows/`, shared IO unions/codecs, `test/runtests.jl` and global docs. |
| Collisions/sources | P02–P03, P35–P39 | Serial collision/source owners first; rank/device composition only after P08/P13. |
| Bounded EM | P04–P07, P32–P34 | Yee field/boundary/coupling owners; no concurrent implicit-EM writer. |
| Distributed/device | P08–P15 | Pencil, CUDA and Metal paths remain separate; sole writer for extension registration and shared array/operator interfaces. |
| Implicit/geometry/AMR | P16–P24 | One writer per family file; shared nonlinear/pressure/current/mesh owners require handoff. |
| Species/control | P29–P31 | Collection/full-PIC/kinematics/checkpoint paths; parser/current edits serialized with P01/P30. |
| Electron transport | P40–P41 | Energy-equation/closure/collision/integrator/checkpoint owners. |
| Validation | P25–P28 | Read-only independent references and case harnesses; production edits coordinated with the owning lane. |

## P01 — current CI and compiler regressions

### Confirmed blocker and repair

`CI-PARSE-001` was verified from the retained raw log
`handoff/2026-11-09/ci/129005385324-loader.log` lines 765–771 and
`ci/MPI_FAILURE_ARTIFACT.json`:

```text
# Error @ src/models/empic_collection.jl:1430:46
report.kind === (es isa EMPIC ? :vector :planar) && ...
# space required after `:` in `?` expression
```

Both invalid symbol-valued ternaries were changed to the explicit form
`es isa EMPIC ? :vector : :planar`:

- `src/models/empic_collection.jl:1430`
- `src/models/empic_collection.jl:1854`

A repository-wide search for the same bare-symbol ternary form found no other
source site. The frozen controller patch already used the correct
`:vector : :planar` form at its two corresponding lines, so P30 can reconcile
the parser repair without changing those two expressions again.

Commit `d95ba9c` also replaces two new implicit-hybrid JET inference sites with
concrete typed locals and arity-specific metric dispatch. The change preserves
runtime behavior: evolving-pressure reports still take the 3-tuple path and
isothermal reports the 2-tuple path.

### Evidence status at this entry

| Item | State |
|---|---|
| Parser repair | Committed in `d95ba9c`; exact-head execution pending. |
| JET additions from `5703c95`/`bc86550` | Historical `5703c95` latest job `128983114996` reported 91 and `bc86550` minimum job `128996436293` reported 99. Their added top-level chains were the annulus/static-composite construction paths already repaired by `566ddc1`, plus the implicit-hybrid pressure-report policy repaired in `d95ba9c`. Exact-head recount pending. |
| Formatting | `format-candidate` workflow dispatched for `d95ba9c`; artifact must be reviewed and applied before the real formatting job can pass. |
| Focused suites | `source-representation` (includes `test_jet.jl`), `persistent-species` (parser path) and `implicit-hybrid-pressure` dispatched for `d95ba9c`. |
| Full CI | Push run `42599058618` queued for `d95ba9c`. |

## 2026-11-09 — controller adoption and AMR constructor repair

| Item | Record |
|---|---|
| P30 source adoption | `dd2c268` applies the frozen nine-path controller patch with `git apply --whitespace=error-all`, registers `models/empic_coupled_error_control.jl` after `models/empic_collection.jl`, the diagnostic after `diagnostics/relativistic_moments.jl`, and the schema-4 IO after `io/empic_collection_checkpoint.jl`. Adds `test_empic_coupled_error_control.jl` to group 10 and a `coupled-empic` focused suite. |
| P30 parser reconciliation | After apply, the report ternaries are at `empic_collection.jl:1480/1904` and use the correct explicit `:vector : :planar` form. No bare-symbol ternary remains in the nine paths. |
| P30 source limits | Runtime, formatter, JET, native replay, scientific/effectivity and resource acceptance remain open. The constructor still rejects NGP; no NGP jump/error policy is claimed. |
| P23/P24 confirmed defect | Historical jobs `128999430606`/`128999430909` at `57980ae` showed four `MethodError`s constructing `_CartesianAMRPatch` with a `BitVector` valid mask. The parametric default constructor does not convert that broadcast result to the declared `Vector{Bool}`. |
| P23/P24 repair | `93acac8` uses `map(!, covered)`, which preserves the Boolean mask as `Vector{Bool}`. Focused `adaptive-composite-es` and `periodic-adaptive-composite-es` dispatched for the new head. |
| G6 hardening | `93acac8` also requires the saved coupled checkpoint `integrator`, `particle_layout` and (for `EMPIC1D`) `polarization` values to match the destination, rather than checking only key presence. |
| P00 evidence | `docs/OPENCODE_SUPPORT_MATRIX.md` now has family/geometry/closure/species/source/backend/precision/persistence rows. `docs/OPENCODE_PACKAGE_STATUS.md` is the live 44-package status. |

## 2026-11-09 — focused runtime results

| Work | Exact head | Workflow/job | Actual runtimes | Result |
|---|---|---|---|---|
| Adaptive Cartesian composite ES | `93acac8` | `Implementation essentials` run [42602594357](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42602594357), jobs 129110455184/129110455209 | Julia 1.12.7 and 1.13.3 | pass, including `test_adaptive_composite_electrostatic.jl` and the static/dynamic consumers |
| Periodic adaptive composite ES | `93acac8` | run [42602603461](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42602603461), jobs 129110480539/129110481234 | Julia 1.12.7 and 1.13.3 | pass, including the full periodic hierarchy/solver/restart essential |
| Coupled projected Euler–Heun EMPIC | `0bf27af` | run [42612416738](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42612416738), jobs 129143913160/129143912693 | Julia 1.12.7 and 1.13.3 | pass, 3556/3556 on latest; full independent reference, cuts/seams, rejection/rollback, schema-4 restart, relativistic C1 and view destinations |
| Exact-head JET recount | `5cb3781` | `source-representation` run [42631633211](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42631633211) | Julia 1.12.7 and 1.13.3 | package-located reports = 0; total reports 118/127 versus cap 86 — **still open** |

The focused results establish the tested slices, not whole-package science or
resource acceptance. The JET cap remains a P01 blocker: every remaining
dependency/Base chain must be classified and reduced or the source must be
restructured; the cap is not raised.

### Full-CI group-5 CGL regression

The full run `42612387035` at `0bf27af` completed all groups except group 5,
group 11 and formatting. Group 5 failed twice (Julia 1.12.7 and 1.13.3) at
`test/test_collisional_cgl.jl:328`: after a direct `_restore_checkpoint_state!`
and one continuation step, the reference state still reported
`clock=(0.001,1,0.001)` while the loaded target reported
`(0.002,2,0.001)`. Explicit pre-step clock assertions were added at
`4c6cd45`; the focused `collisional-cgl` run `42649127988` then passed the
three restore-clock assertions but still produced the same final mismatch.
The restore clock is therefore correct; the continuation path for the direct
`_restore_checkpoint_state!` reference does not advance while the loaded target
does. This is a confirmed open P35/P01 failure, not a CGL pass. The next test
must instrument `reference.dt`, `reference.inner.time` and `reference.phase`
around the step itself and compare with the target path.

### CI queue observation

Focused and full runs repeatedly sat in the GitHub-hosted queue for tens of
minutes to hours while older inherited runs drained. This is an observed
service queue, not a billing exception and not evidence of a code result. No
inherited run was cancelled; only OpenCode-created superseded runs were
cancelled.

## Update rule

Append dated entries with commit SHA, exact command/workflow URL, raw artifact
path or checksum, actual Julia/runtime/dependency versions, pass/fail and the
gate affected. Never replace a failed observation with a later success; add a
superseding record. Keep failed candidates and unexecuted hypotheses explicit.

### JET whole-package cap closure (2026-11-10)

The accepted baseline was re-measured instead of assumed. A temporary audit
branch (`audit/jet-baseline-86`, pushed from `main` at `9b4d363`) ran
`test/test_jet.jl` unchanged: [run 42765960739](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42765960739)
reports **86** whole-package reports on Julia 1.12.7 and **79** on 1.13.3, with
package-located reports **0**. That reproduces the recorded cap measurement with
the same JET 0.13.1 and dependency set.

At `8ae9b27`, the unchanged gate reports **86** on Julia 1.12.7 and **77** on
1.13.3, package-located **0**:
[run 42768819449](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42768819449)
(jobs `129671577789` and `129671577443`). The cap, the package-located gate, the
comparison, the tolerances and the historical ledger are unchanged; the count
was reduced only by concrete-typing and dispatch repairs.

A root diff against the baseline (1.12.7) shows the branch removed five accepted
chains and added only adaptive composite electrostatic chains. Every added chain
was repaired: exact-key coupled field accessors, per-species force/lift/row
guards, typed adaptive model/diagnostic/operator owners, `CartesianIndices` span
and buffer enumeration, exact publish component vectors, typed static/dynamic
AMR and Yee stage constructors, and the retained-history/step-zero guards. The
remaining reports are the accepted dependency/Base baseline.

### Full-CI group-5 CGL focused pass

The direct-restore no-op was repaired in the test by continuing the reference
through the public loader while keeping the direct-restore clock assertions.
Focused `collisional-cgl` [run 42711493910](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42711493910)
passes on both runtimes at `c43d539`.

### Exact-head full CI pending

Full CI at `0d6734b` (formatting-only successor of `8ae9b27`) was still queued
when this entry was written. Its group 5, group 11, formatting, MPI and
CUDA-host jobs are the next P01 acceptance evidence; no result is asserted here.

### Full-CI group 5/11 formatting pass and the group-7 regression

Full CI [run 42773113428](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42773113428)
at `cd2e811` records: group 5 (CGL) **success** on both runtimes, group 11
(JET) **success** on both, `Formatting check (JuliaFormatter)` **success**, both
CUDA host loads **success**, and 2-/4-rank real MPI **success**. The 8-rank job
was still pending when the superseded run was cancelled, so it is not claimed.

The only failure was group 7 (`test_io.jl`, 37,022 passed / 2 errored): my
strict `keys(fields) == keys(saved)` check rejected saved payloads that carry an
extra externally edited key (`extra::Symbol`) or reordered keys. That check was
a regression, not a gate change. `_restore_checkpoint_fields!` now selects the
destination keys explicitly and copies each field by name, so extra keys and
reordering are accepted while a missing destination key still rejects.

The final exact-head full CI is the next evidence; nothing is claimed from the
cancelled superseded run beyond the completed jobs above.
