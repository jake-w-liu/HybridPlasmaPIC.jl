# Plasma simulation audit and implementation tracker

Updated: 2026-10-14. Goal remains **active**.

User objective: audit the full codebase until an independent pass finds no new
confirmed defects; optimize measured bottlenecks; plan and implement missing
capabilities toward a leading, comprehensive plasma simulation tool. Compatibility with the latest stable Julia is an explicit acceptance gate. Every
verified fix or implementation is committed and pushed separately to `main`.
Passing a subset of tests does not establish whole-project correctness or
completeness. Comparative leadership requires reproducible external benchmarks.

## Current evidence

- Starting revision: `905ea1b`; clean worktree; remote `main` matched that revision.
- Inspected package exports, manifests, CI, contribution rules, README, tests,
  validation catalog, hybrid integrators, coupling kernels, and checkpoint paths.
- Installed Julia versions are 1.12.7 and 1.13.1. The official Julia downloads page
  [lists 1.13.1 (2026-09-25)](https://julialang.org/downloads/manual-downloads/)
  as the current stable release. CI now includes both the
  minimum supported minor (1.12) and the moving latest stable version (`1`).
- Initial full suites were started in the live worktree on Julia 1.13.1 and
  1.12.7. Concurrent test/code edits mean these runs cannot prove the final
  revision; both were stopped through their own launch sessions after this was confirmed.
  Final full-suite checks use immutable checkouts. The initial 1.13.1 SHA
  compatibility finding is fixed; HDF5 and PencilFFT extensions have now
  precompiled successfully in both fresh version runs.
- Full package suite at fixed `1138317` on Julia 1.13.1: **106,662/106,662 tests
  passed**. Julia 1.12.7: **106,661 passed, 1 failed**; the sole failure was
  the JET report-count gate (31 findings against a ceiling of 30). The checkpoint
  inference repair passes focused JET checks on both versions without increasing
  the ceiling. Full-suite runs at immutable `c76e39c` passed **117,438/117,438
  tests on each of Julia 1.12.7 and 1.13.1** (48m56.7s and 37m37.6s respectively;
  both processes exited successfully). Later fixes have separate focused
  evidence and still require a final integrated run. Logs:
  `/tmp/hybrid-full-1138317-julia112.log` and `/tmp/hybrid-full-1138317-julia113.log`.
  New-run logs: `/tmp/hybrid-full-c76e39c-julia112.log` and
  `/tmp/hybrid-full-c76e39c-julia113.log`.
- Julia 1.13.1 quick validation at fixed revision `b200c76`: **150/150 metrics
  passed across 20 cases**, no skipped results. Command:
  `julia +1.13.1 --project=validation validation/run_validation.jl --quick --no-plots`.
  Evidence: `/tmp/hybrid-validation-b200c76.log` and
  `/tmp/hybrid-validation-b200c76-results/artifacts/validation_summary.csv`.
- Julia 1.13.1 full electromagnetic suites at fixed `b62da7b`: **7,707/7,707
  assertions passed** across 40 test sets, including the new subcycling and
  Float32 cases. Log: `/tmp/empic-subcycling-b62da7b-julia113.log`.
- Julia 1.13.1 combined I/O, reproducibility and MPI checks at fixed `9856971`:
  **4,820/4,820 assertions passed** across 44 test sets. This includes the new
  two-species EM restart implementation. Log: `/tmp/hybrid-io-9856971-julia113.log`.
- CI run `38414172669` passed formatting and real MPI smoke at 2, 4, and 8 ranks.
  Its package suites were superseded by subsequent required pushes. Local full
  suites continue independently of push-triggered CI cancellation.
- Aqua passed all 9 quality checks on Julia 1.13.1 at `b62da7b`;
  the full JuliaFormatter 1.0.62 scan passed too.
- Baseline kernel harness: derivative 2.67 µs / 0 B; CIC deposit 305.29 µs for
  12,800 particles; hybrid step 1.54 ms for 3,200 particles. These are local
  measurements, not portable performance targets.
- Fused vector gather: all 304 deposition/gather assertions passed on
  Julia 1.12.7, including exact scalar equivalence and zero allocation. The repository benchmark
  measured 2.41–3.28× speedups against three scalar gathers with zero bytes allocated in
  every 1D–3D/NGP/CIC/TSC case. Hybrid and CAM-CL regression suites also passed.

## Audit coverage matrix

Each row requires domain/algorithm, implementation/error-path, resource, and
performance review. “Pending” means no current-run completion claim.

| Surface | Current work | Evidence needed to close |
| --- | --- | --- |
| Grids and spectral operators | Inventory complete; detailed pass pending | Manufactured solutions, boundary/Nyquist contracts, dependency tests |
| Particle storage/loading/pushers/boundaries | Coupling and step timing traced; remaining pass pending | Conservation, order, boundary maps, particle identity, invalid state |
| Deposition/gather/moments | Vector gather optimization verified; moment audit remains open | Adjointness, all shapes/dimensions/types, scalar equivalence, allocations |
| Electron closures and Ohm law | Integration traced; detailed pass pending | Tensor/force oracles, conservation budgets, dissipative terms |
| Hybrid/CAM-CL integration | Variable-step repair verified, including coupled-field convergence | Analytic timing, temporal order, subcycles, restart, MPI equivalence |
| Electrostatic/electromagnetic PIC | Timing, subcycling and restart verified; remaining audit open | Constant-force trajectory, Maxwell waves, charge continuity, mobile species |
| Hall-MHD, gyrokinetics, AMR | Hall-MHD domain/transaction/performance pass completed with no further confirmed defect; gyrokinetic position repairs active; AMR detailed pass pending | Integrate latest Hall fixes into final full suites; complete gyrokinetic/AMR passes |
| Collisions and radiation | Prior fixes exist; fresh pass pending | Conservation/rates, sampling, unequal weights, numerical robustness |
| Shocks/injection/campaigns | Inventory complete; detailed pass pending | RH laws, inflow flux, multidimensional/restart and published comparisons |
| Diagnostics and ray/Raycon solvers | Inventory complete; detailed pass pending | Independent references, units, transport/conversion invariants |
| Checkpoints/archives/field I/O | Independent audit active; defects confirmed | Configuration identity, exact bytes/checksums, atomic failure, restart |
| Threads/MPI/distributed FFT/GPU | Interfaces inspected; full pass pending | Real multirank tests, backend equivalence, ownership, hardware benchmarks |
| Packaging/quality/examples/validation | Full suites pass 117,438 assertions on each Julia version at `c76e39c`; later changes have focused checks | CI-version tests, Aqua/JET findings, formatting, examples, validation cases |

## Verified findings and fixes

| ID | Finding / acceptance criterion | State | Evidence / next action |
| --- | --- | --- | --- |
| PERF-001 | Vector gather repeats the same stencil three times | Pushed (`58fa769`) | 304 focused assertions plus hybrid/CAM-CL suites pass; reproducible speed/allocation comparison in `benchmark/benchmarks.jl` |
| IO-001 | `save_run` rejects its own valid dictionary output | Pushed (`d1db13d`) | Exact serialized bytes hashed under schema 4; 159 metadata/archive/field/async assertions pass, including cross-process and corruption cases |
| TIME-001 | ES/EM PIC accepts changed `dt` but uses the wrong half-step interval | ES pushed (`def6548`); EM pushed (`66ae6f3`) | ES suite passes including 180 new assertions; full EM suites pass 5,530 assertions; variable-step Maxwell/gyro convergence and charge/Gauss checks pass |
| IO-002 | Same-grid checkpoint silently accepts different model parameters | Pushed (`236d5d9`) | Restore saved model/shape parameters, reject incompatible concrete types; custom mutable configuration and alias tests pass |
| IO-003 | Truncated field file allocates declared payload before checking available bytes | Pushed (`c7b1c6f`) | 76 field/async/metadata assertions pass, including malformed extents, bounded allocation, empty and 0D arrays |
| IO-004 | HDF5 silently writes incorrect values from column-gapped views; range-backed arrays throw | Pushed (`b228b3f`) | 61 extension assertions pass, covering exact data/type/shape, contiguous-buffer identity and error recovery |
| TIME-002 | Hybrid/CAM-CL/MPI use incorrect kick interval after changing dt | Pushed (`236d5d9`) | 277 phase/convergence/MPI assertions; hybrid/CAM-CL/MPI regressions and real 2-rank suite pass; fixed-step comparison with original code is bitwise equal |
| IO-005 | Checkpoint restore fails on matching-length particle views and can partially mutate destination | Pushed (`236d5d9`) | 271 checkpoint/archive assertions pass, including 60 view/unsupported-resize assertions |
| PKG-001 | SHA compatibility excluded Julia 1.13’s bundled SHA 1.0 | Pushed (`4b5611d`) | Installed-version compatibility check plus 60 metadata assertions pass on each of Julia 1.12.7 and 1.13.1 |
| RESTART-001 | CAM-CL lacked checkpoint/archive methods and requires persistent magnetic state | Pushed (`1138317`) | 710 checkpoint/archive/MPI assertions pass, including 1D–3D, Float32/64, phase resets, edited fields, timestep/subcycle changes, malformed-state rejection, and exact continuation |
| CI-001 | Only Julia 1.12 was exercised; the test pipeline could mask failures | Pushed (`68b854d`) | YAML/matrix/artifact wiring checks pass; explicit CI Bash returns 1 for failed pipelines and 0 for successful ones; full version suites running separately |
| FORMAT-001 | Three files failed CI's JuliaFormatter 1 formatting rules | Pushed (`cf53bc4`) | Read-only scan with JuliaFormatter 1.0.62 identified exactly three Julia files; formatting-only repairs are idempotent and their diff was inspected. Contributor instructions now specify the CI formatter major |
| CFL-001 | Timestep recommendation lost representable results through intermediate overflow/underflow | Pushed (`27b4c8d`) | Valid 4-cell grid of length 1e-80 returned zero instead of the 512-bit oracle's 1.4184965709927286e-162. Logarithmic evaluation passes 3,802 helper/range/coupled-stability assertions on each of Julia 1.12.7 and 1.13.1 (`test/test_phase_helpers.jl`) |
| RESTART-002 | Electrostatic PIC lacked checkpoint/archive support | Pushed (`b200c76`) | Julia 1.13 combined IO/repro/MPI: 1,194 assertions; Julia 1.12 final focused suite: 540 assertions. Covers Float32/64, 1D–3D, exact variable-step continuation, changed fields, invalid destinations, views, empty populations, and archives |
| SUBCYCLE-001 | General EMPIC lacked electron subcycling | Pushed (`c5fef77`) | Julia 1.12.7 full EM suites: 7,257 assertions; added permanent Float32 coverage: 450 assertions in a separate final-file run. Includes independent constant-field oracles, 1D equivalence, 2D/3D charge/Gauss checks, changing dt/counts, mobile-ion coupled convergence. Default n_sub=1 arrays match prior code; allocations unchanged and scale linearly with subcycles |
| GYRO-001 | Guiding-centre drift loses finite results by reconstructing perpendicular speed squared from magnetic moment | Pushed (`2eeecd6`) | Direct moment-based grad-B drift; independent 2,048-bit equations cover extreme and ordinary Float32/64 inputs, charge signs, mass, zero moment, and failed-update rollback. Fusion suite: 256 assertions passed on each of Julia 1.12.7 and 1.13.1, including 76 new assertions; warmed push remains allocation-free |
| RESTART-003 | EM PIC lacked restart/archive support for two species and their staggered phases | Pushed (`9856971`) | Julia 1.12.7 combined I/O/reproducibility: 4,692 assertions, including 3,620 EM assertions. Covers legacy/general 1D–3D, Float32/64, all shapes, mobile/immobile and relativistic/nonrelativistic runs, exact continuation, diagnostics, changing dt/subcycles, invalid-state rollback, views and empty populations. Log: `/tmp/hybrid-io-empic-final112.log` |
| HALL-001 | Final Hall-MHD projection/Ohm evaluation can throw after changing public fields | Pushed (`f1f946e`) | Candidate fields are fully evaluated before copying to the public state. 114 Hall assertions pass on each of Julia 1.12.7 and 1.13.1; 108 independent before/after comparisons on 1.12 cover 1D–3D, Float32/64, three closures and projection on/off, with bitwise-equal successful steps and unchanged allocations. Failure tests compare every public field and preserve array identities. Logs: `/tmp/hall_transaction_final_tests.log`, `/tmp/hall_transaction_final_comparison.log`, `/tmp/hall-transaction-julia113.log` |
| TIME-003 | Finite solver steps overflow simulation time or wrap the integer step counter | Pushed (`f437600`) | Shared preflight covers hybrid, CAM-CL, Hall-MHD, both EM families and MPI. Clock and variable-step regressions pass 953 assertions on each of Julia 1.12.7 and 1.13.1; covers 1D–3D, Float32/64, moving particles, evolving fields, whole-state rollback, zero steps, Int boundary and EM substep rounding. Logs: `/tmp/hybrid-clock-final112.log`, `/tmp/hybrid-clock-final113.log` |
| JET-001 | Dynamic checkpoint field names widen inference to unrelated fields; Julia 1.12 full-suite JET gate fails | Pushed (`c76e39c`) | Typed field bundles and explicit ion-presence branches retain the existing schema. Frozen repair: JET reports 24 on Julia 1.12 and 22 on Julia 1.13, both below the unchanged ceiling of 30; no checkpoint reports. Concrete capture/validate/restore checks report zero across six families on both versions. I/O/repro: 4,692 assertions pass; eight serialized fixtures remain byte-identical. Logs: `/tmp/hybrid-io-typed-fields-*.log` |
| HALL-002 | Hall-MHD rejects stationary density exactly at the shared Ohm-law floor | Pushed (`d950dde`) | Inclusive bound matches the floor operation. Hall suite: 402 assertions pass on each of Julia 1.12.7 and 1.13.1; 288 new assertions cover exact/adjacent representable densities, three closures, 1D–3D, Float32/64, unchanged stationary fields and rejected-state rollback. Logs: `/tmp/hall_inclusive_floor_julia112.log`, `/tmp/hall_inclusive_floor_julia113.log` |
| HALL-003 | Final magnetic projection can produce an infinite electric field that Hall-MHD accepts | Pushed (`f96502c`) | Validate final E before committing candidate fields. Hall suite: 558 assertions pass on each of Julia 1.12.7 and 1.13.1; 156 new checks include unchanged rejected state, an unprojected retry and a finite projected neighbor. Independent comparison: 108 successful steps match prior fields/clock exactly, with unchanged allocations. Projected nonfinite B is already rejected by the spectral curl input guard (16 additional checks). Logs: `/tmp/hall_final_e_final_julia112.log`, `/tmp/hall_final_e_final_julia113.log`, `/tmp/hall_final_e_final_comparison.log`, `/tmp/hall_projected_b_guard.log` |
| RESTART-004 | Hall-MHD has no checkpoint/archive interface and carries inactive stage history | Pushed (`5c47d75`) | Add no-particle methods preserving all public fields and the retained stage `lapJ`. Focused Hall tests: 2,418 assertions on each Julia version; combined I/O/repro on 1.12: 7,110 assertions. Covers Float32/64, 1D–3D, three closures, initialization/live states, model changes, mutable closure ownership, array identity, rejection rollback and exact continuation. JET remains 24 reports on 1.12 and 22 on 1.13; Hall capture/validate/restore concrete analysis has zero reports. Logs: `/tmp/hybrid-io-hall-regressions112.log`, `/tmp/hybrid-io-hall-focused113.log` |
| GYRO-002 | Parallel velocity loses representable results by rounding field projection, force, acceleration and final addition separately | Pushed (`e9991e9`) | Compensated binary factors retain the complete expression; exact rational evaluation handles uncertain cancellation/range boundaries. Fusion: 550 assertions pass on each Julia version; 59,623 independent certificate/range/physical assertions pass on 1.12, including Float16 residual-underflow witnesses. Measured nonstationary and certified-equilibrium paths allocate 0 B; an uncertified mass=3 equilibrium still allocates 8,296 B. Local nonstationary latency increased about 20%; certified-equilibrium latency increased 36–44%. Logs: `/tmp/gyro_certified_final_112.log`, `/tmp/gyro_certified_final_latest.log`, `/tmp/gyro_certified_interleaved_performance.log` |
| PERF-002 | Hall-MHD repeats forward transforms for each velocity-gradient axis | Pushed (`870b544`) | Shared gradient in 2D/3D removes 12/24 forward FFTs per RK4 step. Hall tests: 627 assertions pass on each Julia version, including independent helical-wave, fourth-order convergence and nonlinear advection oracles. 108 normal-grid trajectory comparisons are exact; extreme anisotropic scaling can round differently. Benchmark `benchmark/hall_mhd_gradients.jl`: median baseline/current CPU ratios 1.085 at 128² and 1.133 at 64³; wall results vary, including mixed 2D repeats. Allocations increase 384 B/step in 2D/3D, unchanged in 1D. Logs: `/tmp/hall_gradient_final_julia112.log`, `/tmp/hall_gradient_final_julia113.log`, `/tmp/hall_gradient_final_comparison.log`, `/tmp/hall_gradient_repo_benchmark.log` |
| MSPEC-001 | Hybrid/CAM-CL accept one population despite existing multispecies moment deposition | Verified; commit in this change | Serial collection initialization/evolution preserves charge-weighted shared moments and independent q/m. Both Julia versions pass 2,286 integrated assertions, including 918 new collection checks. Independent uniform orbit, Boris force balance and wave-matrix convergence pass; 180 old/new scalar states agree exactly. Latest-Julia package JET remains 22 reports; all 24 concrete collection calls have zero reports; Aqua passes 9 checks. Example passes both versions with maximum velocity error 2.162e-8. Benchmark reports 18,992–35,456 B/two-ion step on Julia 1.13; alias checks account for 9,216–19,200 B and are a measured optimization target. Scalar extraction adds 16 B in 1D and 32 B in 2D/3D on Julia 1.12; latency samples vary. Logs: `/tmp/hybrid-multispecies-reviewed112.log`, `/tmp/hybrid-multispecies-reviewed113.log`, `/tmp/hybrid-multispecies-core-jet113.log`, `/tmp/hybrid-multispecies-benchmark113.log`, `/tmp/hybrid-multispecies-scalar-comparison.log` |

## Completeness roadmap

These are capability gaps or assessment tasks, not all classified as bugs.
Inspect implementation contracts before choosing the design; add an acceptance
benchmark before implementation. The roadmap remains open to new findings.

| Priority | Capability / gap | Plan and completion evidence | State |
| --- | --- | --- | --- |
| 1 | Reliable variable-step evolution | Correct phase state in ES/EM, hybrid/CAM-CL and MPI; preserve it in checkpoints; verify analytic trajectories and convergence | Active |
| 1 | Restart across production solvers | Hybrid, CAM-CL, electrostatic/electromagnetic PIC and Hall-MHD share capture/validation/restore infrastructure; exact continuation includes EM species and Hall retained stage history | Existing solver families verified; multispecies hybrid extension pending |
| 1 | Reproducible archives and bounded field I/O | Fix checksum representation and malformed-file allocation; cross-process and corruption tests | Active |
| 2 | Multispecies hybrid evolution | Collection initialization/stepping absence reproduced; [implementation plan](MULTISPECIES_HYBRID_PLAN.md) defines shared charge moments, independent species push, restart and conservation/dispersion acceptance | Serial core verified; collection restart and diagnostics active |
| 2 | Multispecies conservation diagnostics | Runtime probe confirms collection methods are absent for kinetic energy and energy/momentum budgets, and existing collection total momentum errors on a typed empty vector. Extend canonical reductions and count field reservoirs once; verify unequal masses/weights, empty collections and closure energies | Next core slice |
| 2 | Exactly conservative coupled hybrid integration | Existing explicit predictor has second-order momentum error verified against a discrete Boris oracle. Assess a coupled midpoint field/current solve with convergence, conservation and nonlinear solve-residual acceptance | Assessment |
| 2 | Evolving electron pressure integrated into steppers | Standalone `advance_electron_pressure!` exists; built-in periodic steppers use algebraic closures. Define transport/heating coupling and energy-budget acceptance | Planned |
| 2 | Multidimensional EM electron subcycling | Common EM path extended; charge-continuity, Gauss-law, species phase, coupled convergence and Float32/64 checks pass | Verified |
| 2 | Nonperiodic EM-PIC boundaries | README declares periodic EM-PIC only; design compatible field/particle boundaries with reflection/transmission and charge tests | Planned |
| 3 | Production distributed stepping | MPI reference stepping replicates fields; design domain-owned fields/moments with migration and distributed operators; verify serial equivalence and scaling | Planned |
| 3 | Production GPU kernels | README claims array plumbing only; implement coupled gather/push/deposit/field paths with backend equivalence, memory and throughput measurements | Planned |
| 2 | Spatially local, interspecies Coulomb collisions | Current public Coulomb operator takes one ParticleSet; design cell-local multispecies pairing with momentum/energy and equilibration-rate oracles | Planned |
| 2 | Integrated implicit field advancement | Current `run_whistler` CN implementation evolves a linear reference packet; assess nonlinear implicit hybrid/EM integration with stability, accuracy and solve-residual benchmarks | Assessment |
| 3 | Model/geometry completeness | Assess Hall-MHD, gyrokinetic, AMR and fusion building blocks for missing end-to-end coupling before planning extensions | Assessment |
| 2 | Guiding-centre coordinate range and cancellation cost | Confirmed coordinate update loses finite results through normalized-field underflow and displacement overflow; repair the complete position expression. Assess allocations for uncertified exact force balance without weakening numerical guarantees | Active |
| 3 | Published validation and external comparisons | Run local validation; inspect full Preisser profile replay and external solver cases; add reproducible accuracy/cost comparisons | Planned |
| 4 | Independent final audit | Revisit every matrix row from another entry point, resolve all confirmed defects, rerun final acceptance gates | Pending |

## Verification commands

- Full package/quality suite in immutable checkouts, on both Julia 1.12 and latest
  stable: `julia +1.12 --project=. -e 'using Pkg; Pkg.test()'` and
  `julia +1.13.1 --project=. -e 'using Pkg; Pkg.test()'`.
- Coupling: `julia +1.12 --project=. test/test_deposit.jl`.
- Kernel timing/allocation: `julia +1.12 --project=. benchmark/benchmarks.jl`.
- Validation: `julia +1.12 --project=validation validation/run_validation.jl --quick`
  followed by the complete case selection after inspecting runtime/dependencies.
- MPI: CI invokes `test/test_mpi_multirank.jl` at 2, 4 and 8 ranks using
  `MPI.mpiexec()`; record each observed result.
- Final diffs: `git diff --check`; format changed Julia files with the configured
  JuliaFormatter version; inspect the exact staged files before each commit.

No full-codebase, complete-physics, maximum-performance, GPU-production, or
cluster-scaling claim has been established by this run yet.
