# OpenCode package status

Updated 2026-11-09. This is the live operator status for the 44 work packages.
It does not rewrite the frozen `handoff/2026-11-09/WORK_PACKAGES.json`; that
packet remains the historical handoff. `open` still means at least one
applicable G1–G10 gate is missing.

Status values: `open`, `in progress`, `source reviewed`, `runtime accepted`,
`science accepted`, `resource accepted`, `complete`, `blocked`.

| Package | Status | Current evidence / next gate |
|---|---|---|
| P00 Freeze and reconcile | in progress | Packet verified PASS; seven primary hashes match; dirty primary untouched; support matrix created. Dependency/runtime versions and final owner sign-off still need exact-head records. |
| P01 Close CI/compiler regressions | in progress | Parser repair `d95ba9c`; AMR constructor repair `93acac8`; package-located JET now 0 at `5cb3781`. Remaining blocker: whole-tree JET total 118/127 vs cap 86, formatter, and full exact-head CI. |
| P02 Serial collision resources | open | Source + focused tests; resource/occupancy campaign not run. |
| P03 Rank/device collisions | open | Depends P02/P08/P13; no real-rank/device acceptance. |
| P04 Bounded Yee families | open | Source + focused tests; independent cavity/trajectory/work and event scope pending. |
| P05 General material boundaries | open | D2/3 crossed-cell/face/edge/corner and repeated events still missing. |
| P06 Absorbers/open coupling | open | PML/conductive source exists; reflection spectra, corners and coupled plasma exits pending. |
| P07 Nonperiodic production composition | open | Depends P05/P06/P10/P14/P31. |
| P08 Distributed layout | open | Primitive operators focused; broader layout/conditioning and real-rank contract pending. |
| P09 Pencil solver lifecycles | open | Some 2-rank model controls; production families/restart pending. |
| P10 Broader distributed families | open | Depends P09/P35. |
| P11 Repartition and scaling | open | Changed layout/rank and multi-node scaling pending. |
| P12 GPU primitive contracts | open | Host admission only; actual device kernels/range/failure pending. |
| P13 Resident first solvers | open | CUDA source + host controls; actual device lifecycle and Metal resident solver pending. |
| P14 Resident physics families | open | Depends P13/P35. |
| P15 Multi-GPU | open | Not implemented. |
| P16 Implicit audit/common owner | open | Source + focused; redundant FFT/Krylov optimization and true-residual acceptance pending. |
| P17 Coupled implicit hybrid | open | Source + focused; independent nonlinear/order/restart/science pending. |
| P18 Coupled implicit EM | open | Source + focused; D2/3/boundary/relativity/backend and science pending. |
| P19 Shear/toroidal annulus | open | Source + focused; weighted adjoint, boundary work and resource pending. |
| P20 Full curvilinear/toroidal | open | Axis/material/equilibria/imports/backends missing. |
| P21 Slab gyrokinetics | open | Source + focused; collective mode/FLR/free-energy/science pending. |
| P22 Broader gyrokinetics | open | Geometry/collision/kinetic-EM/backends missing. |
| P23 Adaptive ES boundaries | in progress | Constructor defect fixed in `93acac8`; adaptive and periodic-adaptive focused suites passed both Julia runtimes (runs 42602594357/42602603461). Material/electrode/mixed-axis, independent refinement and resource gates remain open. |
| P24 Adaptive hybrid/EM/subcycling | open | Hybrid/EM AMR, subcycling and device hierarchy not implemented. |
| P25 Pinned references/wave matrix | open | Existing campaigns; production oblique/kinetic/thermal/EM matrix and uncertainty pending. |
| P26 Nonlinear profiles/datasets | open | Matched datasets/builds and full profiles pending. |
| P27 Cross-family validation | open | Independent bounded EM/collision/GK/geometry/AMR/reconnection campaigns pending. |
| P28 Error-matched effectiveness | open | No frozen comparative campaign. |
| P29 Arbitrary species/schedules/IO | open | Source + focused; coupled controller adoption is P30. Full pending-deadline/failure/raw-archive acceptance pending. |
| P30 Coupled controller | in progress | Nine-path source adopted in `dd2c268`; schema-4 tag hardening `93acac8`; focused essential passed both Julia runtimes at `0bf27af` (3556/3556 latest, run 42612416738) including reference, rollback, restart, relativistic C1 and view destinations. NGP jump, true-error/effectivity, replay portability, science and resource gates remain open. |
| P31 Schedule/source/boundary/backend intersections | open | Depends P30/P35/P05/P09/P14. |
| P32 Full-Hall open shock | open | Source + focused; physical boundary regime and sweeps pending. |
| P33 Open pressure/transport/inertia | open | Depends P32/P40/P41. |
| P34 Open shocks distributed/device | open | Depends P33/P10/P14. |
| P35 Composed sources/thermal/CGL | open | Source + focused; full CI group 5 exposed a CGL direct-restore/continuation clock mismatch at `test_collisional_cgl.jl:328` (`0bf27af`), with diagnostic assertions pending a focused rerun (`4c6cd45`). Integrated exact-head, native energy/entropy and resource gates remain open. |
| P36 Reaction growth/namespaces | open | Primitives exist; complete reaction/neutral/growth composition pending. |
| P37 Variable/higher-order source | open | Post-field stochastic split remains first order; order/retry/restart pending. |
| P38 Accounted RF/ray/radiation | open | Standalone components; accounted reservoir/coupling/restart pending. |
| P39 Rank/device sources/RNG growth | open | Depends P03/P10/P14/P36–P38. |
| P40 Limited/state-dependent/anisotropic conduction | open | Frozen proposal unapproved; design/domain decisions and implementation absent. |
| P41 Full drag/inertia/tensors | open | Reduced exchange/algebraic CGL/transverse filter are not these models. |
| P42 Whole-code audit/optimization | open | Depends all packages. |
| P43 Compatibility/docs/release | open | Depends P42. |

## Confirmed defect ledger for current work

| Candidate | Dimension | Verification | Status | Fix | Reverify |
|---|---|---|---|---|---|
| `empic_collection.jl` bare-symbol ternaries | parser/correctness | Retained raw log lines 765–771 and MPI artifact | fixed | explicit `: :planar` in both sites (`d95ba9c`) | exact-head package load pending |
| Implicit-hybrid report policy field access | compiler/JET | `5703c95` job 128983114996 full chains | fixed | typed solver local and arity-specific metric dispatch (`d95ba9c`) | exact-head JET pending |
| `_CartesianAMRPatch` `valid` mask type | correctness/runtime | `57980ae` jobs 128999430606/909: four `MethodError` at hierarchy construction | fixed | `map(!, covered)` (`93acac8`) | adaptive focused suites queued |
| Coupled checkpoint tag acceptance by key presence | persistence/robustness | Source audit: `integrator`/`particle_layout`/`polarization` values were not checked | fixed | explicit tag/value admission (`93acac8`) | coupled schema-4 test pending |
