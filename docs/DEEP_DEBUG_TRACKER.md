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
- Full package suite at immutable `67f01b7` on Julia 1.13.1 finished with
  **131,350 passed and 1 failed** in 81m34.8s. The collision rollback assertion
  compared uninitialized position storage containing `NaN` with `==`;
  `test/test_collisions.jl:653` and `src/particles/particle_set.jl:59` identify
  the fixture and allocation contract. The deterministic fixture repair passes all 512 collision assertions on
  each Julia version (`/tmp/hybrid-collision-fixture112.log` and
  `/tmp/hybrid-collision-fixture113.log`). The Julia 1.12.7 run remains active. Logs:
  `/tmp/hybrid-full-67f01b7-julia112.log` and
  `/tmp/hybrid-full-67f01b7-julia113.log`. These snapshots include serial
  multispecies core/restart/diagnostics but exclude later AMR geometry,
  guiding-centre position, alias and MPI collection repairs.
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
| Hall-MHD, gyrokinetics, AMR | Hall-MHD domain/transaction/performance pass completed with no further confirmed defect; gyrokinetic position repair verified; AMR independent dyadic and finite-volume pass verified | Integrate latest Hall fixes into final full suites; complete remaining gyrokinetic public drift/ring checks and assess extreme BigFloat resource bounds |
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
| MSPEC-001 | Hybrid/CAM-CL accept one population despite existing multispecies moment deposition | Pushed (`b57d11e`) | Serial collection initialization/evolution preserves charge-weighted shared moments and independent q/m. Both Julia versions pass 2,286 integrated assertions, including 918 new collection checks. Independent uniform orbit, Boris force balance and wave-matrix convergence pass; 180 old/new scalar states agree exactly. Latest-Julia package JET remains 22 reports; all 24 concrete collection calls have zero reports; Aqua passes 9 checks. Example passes both versions with maximum velocity error 2.162e-8. Benchmark reports 18,992–35,456 B/two-ion step on Julia 1.13; alias checks account for 9,216–19,200 B and are a measured optimization target. Scalar extraction adds 16 B in 1D and 32 B in 2D/3D on Julia 1.12; latency samples vary. Logs: `/tmp/hybrid-multispecies-reviewed112.log`, `/tmp/hybrid-multispecies-reviewed113.log`, `/tmp/hybrid-multispecies-core-jet113.log`, `/tmp/hybrid-multispecies-benchmark113.log`, `/tmp/hybrid-multispecies-scalar-comparison.log` |
| AMR-001 | General scale division introduces nonzero curvature on exact affine profiles | Pushed (`d23dec5`) | Power-of-two normalization and a TwoSum rounding certificate preserve exact curvature; uncertain cases use the existing Rational fallback. BigFloat stored precision is retained, with directed rounding handled exactly. New tests: 7,223 assertions; existing AMR: 62, on each Julia version. Additional 166,013 independent oracle assertions pass on 1.12. Ordinary Float32/64 reductions allocate 0 B; measured local latency rises 19–26%. Logs: `/tmp/amr_curvature_final_julia112.log`, `/tmp/amr_curvature_final_julia113.log`, `/tmp/amr_curvature_final_resources.log`. Geometry and prolongation have separate confirmed findings; the AMR audit remains open. |
| RESTART-005 | Multispecies hybrid evolution requires ordered population restart and archives | Pushed (`ecb8c86`) | Shared capture/restore preserves all species and Hybrid/CAM-CL phase. Full I/O/repro on 1.12: 8,694 assertions; focused final-file collection suite on 1.13: 1,584 assertions, including 144 continuation configurations. Malformed later populations, shared empty storage and incompatible destinations reject before mutation; independent source/destination views and typed empties work. Final production JET remains 24/22 reports on 1.12/1.13, with zero reports in all 12 concrete scalar/collection I/O calls; eight scalar fixtures remain byte-identical. Logs: `/tmp/hybrid-io-collection-verified-regressions112.log`, `/tmp/hybrid-io-collection-handoff-focused113.log`, `/tmp/hybrid-io-collection-verified-jet112.log`, `/tmp/hybrid-io-collection-v4-jet113.log`. Hybrid clock and EM empty-alias repairs remain separate confirmed follow-ups. |
| MSPEC-002 | Collection energy/momentum budgets are absent; empty collection momentum fails | Pushed (`67f01b7`) | Extend canonical budgets and kinetic reduction; preserve heterogeneous nonempty momentum collections. Both Julia versions pass 1,098 diagnostics/regression assertions with independent 512-bit physical sums, all closures/inertia, unequal masses/weights, views and typed empties. Latest-Julia JET remains 22 package reports, all 24 concrete diagnostic calls have zero reports; Aqua passes 9 checks. Warmed energy/momentum particle reductions allocate 0 B for Float32/64 in 1D–3D. Logs: `/tmp/hybrid-multispecies-diagnostics112.log`, `/tmp/hybrid-multispecies-diagnostics113.log`, `/tmp/hybrid-multispecies-diagnostics-jet113.log`. |
| AMR-002 | Cell centers lose finite coordinates through intermediate overflow and rounded half-cell indices | Pushed (`06b637a`) | Complete-coordinate fused evaluation with exact half-index bounds and Rational fallback; Float16 widens to Float64 to avoid double rounding. Both Julia versions pass 9,388 AMR assertions, including 2,103 geometry checks. Additional 250,242 exact-oracle assertions pass on 1.12. Ordinary Float16/32/64 center calls allocate 0 B; local varying-index measurements show no regression. Logs: `/tmp/amr_center_final_julia112.log`, `/tmp/amr_center_final_julia113.log`, `/tmp/amr_center_resources.log`. Prolongation repair evidence follows. |
| AMR-003 | Premature slope rounding erases representable reconstructed children | Pushed (`c083043`) | Complete-child reconstruction shares the curvature rounding certificate; uncertified results use exact stored inputs. Julia 1.12.7 passes 710,102 assertions: 15,198 permanent/existing and 694,904 additional exact-oracle checks. Julia 1.13.1 passes the same 15,198 permanent/existing assertions. This includes both sides of the smallest-normal boundary. Ordinary measured profiles allocate 0 B; extreme cancellation allocates 4,096–8,664 B and boundary fallback 2,120–3,144 B. Logs: `/tmp/amr_prolong_boundary_final_julia112.log`, `/tmp/amr_prolong_boundary_final_julia113.log`, `/tmp/amr_prolong_boundary_final_cpu_resources.log`, `/tmp/amr_prolong_boundary_cold_resources.log`. Stored-precision restriction evidence follows. |
| AMR-004 | Restriction loses stored BigFloat bits before cancellation at lower working precision | Pushed (`ffb23ae`) | Complete dyadic averages retain both child operands. Each Julia version passes 15,468 permanent/existing AMR assertions; 1.12 additionally passes 186,072 independent exact averages and 101 public finite-volume/convergence/ownership/failure checks. Ordinary hardware restriction allocates 0 B; retained BigFloat samples allocate 3,064–12,008 B. Logs: `/tmp/amr_restrict_final_julia112.log`, `/tmp/amr_restrict_final_julia113.log`, `/tmp/amr_restrict_resources.log`, `/tmp/amr_final_public_audit.log`. Public finite-volume checks passed; a later dependency audit found AMR-005. Arbitrary BigFloat exponent-span resource use has not been validated. |
| AMR-005 | Exact-output conversion can round twice through ambient BigFloat precision | Pushed (`e5b082a`) | Shared integer quotient/remainder conversion rounds Float16/32/64 directly, independent of ambient BigFloat settings; BigFloat materialization preserves exact operands and caller rounding. AMR centers, curvature, prolongation and restriction use the shared owner. Corrected independent dyadic/neighbor oracles pass 1,196,587 assertions on Julia 1.12.7 and 65,268 on 1.13.1; 40 concrete JET calls per version have zero reports. The isolated 1.12 package analysis has 26 reports against the unchanged ceiling of 30 (`/tmp/shared_exact_jet_package112.log`); a count within the gate does not classify each report. Ordinary hardware AMR paths remain 0 B; extreme prolongation measures 4,096/4,768/10,272 B for Float16/32/64. Bounded stored-BigFloat samples measure 25,088–112,640 B; arbitrary exponent-span resources remain unverified. Logs: `/tmp/shared_exact_final_julia112.log`, `/tmp/shared_exact_final_julia113.log`, `/tmp/shared_exact_amr_jet112.log`, `/tmp/shared_exact_amr_jet113.log`. Gyro migration is a separate pending change. |
| IO-006 | Shared empty storage escapes alias checks and couples particle resize/write operations | Pushed (`9397d9b`) | Canonical alias checks trace data parents and mutable view indices for empty storage, retaining disjoint nonempty views. Final focused suites pass 1,752 assertions on each Julia version; broad 1.12 regressions pass 4,073 across 51 testsets. Original corruption repros pass 4/4. Package JET remains 24/22, with zero findings in all 48 concrete calls on each version. Two-population validation allocation falls from 9,216–19,200 B to 3,456–7,040 B across 1D–3D; cross-population helper remains 0 B. Logs: `/tmp/hybrid-empty-alias-v3-regressions112.log`, `/tmp/hybrid-empty-alias-handoff-focused113.log`, `/tmp/hybrid-empty-alias-v3-jet112.log`, `/tmp/hybrid-empty-alias-v3-jet113.log`, `/tmp/hybrid-empty-alias-v3-allocation.log`. Intra-population checkpoint aliasing and clock validation remain separate confirmed repairs. |
| MSPEC-003 | MPI hybrid APIs lack ordered ion collection evolution | Pushed (`ac88424`) | Canonical raw moment/push/predictor stages serve serial and MPI paths; collection preflight is collective. Final frozen v3: Julia 1.13.1 real 2-rank collection tests pass 48,049 assertions per rank; 1.12.7 real 8-rank suite passes, including 13,993 collection assertions per rank. Earlier v2 real 4-rank suite passed 25,329 collection assertions per rank; v2 integrated serial/MPI regressions pass 95,567. Scalar before/after comparisons preserve 180 serial and 72 MPI states exactly. Latest-Julia package JET at v2 has 20 reports and Aqua passes all 9 checks; final v3 concrete calls execute successfully with zero init/moment/nonmigrating reports and the same two Serialization reports as the scalar migration baseline. Density-floor narrowing and scratch alias regressions are covered. Two-rank benchmark smoke at 256 proton +128 alpha particles per rank measures 125,376 B maximum local allocation per two-ion step; this is not cluster scaling evidence. Logs: `/tmp/hybrid-multispecies-mpi-v3-2rank113.log`, `/tmp/hybrid-multispecies-mpi-v3-8rank112.log`, `/tmp/hybrid-multispecies-mpi-v2-integrated112.log`, `/tmp/hybrid-multispecies-mpi-v3-concrete-jet113.log`, `/tmp/hybrid-multispecies-mpi-v3-benchmark-2rank113.log`. Distributed collection restart remains pending. |
| MPI-001 | A destination validation failure on one rank leaves other checkpoint destinations restored | Repair under final verification | VERIFIED: the two-rank Julia 1.13.1 probe returns success and changes rank 0 while rank 1 rejects a non-resizable view and remains unchanged. Log: `/tmp/hybrid-mpi-restore-collective-probe.log`. Canonical restore preflight now completes collectively before mutation. Final focused two-rank tests pass on both Julia versions (216 assertions on the failing rank and 258 on the healthy rank); serial restore and compiler regressions are still running. Logs: `/tmp/hybrid-mpi-restore-v2-2rank112.log` and `/tmp/hybrid-mpi-restore-v2-2rank113.log`. The first test snapshot compared copied `Ref` objects by identity and failed only its 42 rollback comparisons per rank; the corrected snapshot compares their values. Details in the multispecies plan. |
| GYRO-003 | Position loses finite updates by rounding drift, displacement and coordinate terms separately | Pushed (`109fe7a`) | Complete scaled displacement retains streaming, electric, magnetic-gradient and curvature terms through the final coordinate addition; exact fallbacks use the shared conversion owner. Each Julia version passes 2,422 fusion assertions; independent oracles add 157,888 on 1.12.7 and 245,045 on 1.13.1. Each version passes 60 concrete JET calls with zero reports and 76 runtime checks; matched package analysis retains 24 reports. All 16 measured ordinary hardware cases remain 0 B. Alternating local batches measured mixed-field costs 30–76% higher than `e9991e9`; equilibrium costs range from 21% lower to 17% higher. These shared-machine timings are not portable limits; further optimization remains open. Bounded retained-BigFloat pushes allocate 88,120–121,800 B; arbitrary exponent-span resources are unverified. Logs: `/tmp/gyro_shared_frozen_regressions112.log`, `/tmp/gyro_shared_frozen_regressions_latest.log`, `/tmp/gyro_shared_final_jet112.log`, `/tmp/gyro_shared_final_jet_latest.log`, `/tmp/gyro_shared_final_performance100k_settled.log`, `/tmp/gyro_shared_final_bigfloat_allocations.log`. Public drift and ring-averaging follow-ups remain open. |
| TEST-001 | Collision rollback test compares uninitialized positions containing NaN | Verified; commit follows | VERIFIED: the 1.13.1 full milestone log shows identical position vectors containing NaN failing `==`; the ParticleSet constructor allocates positions with `undef`. Initialize the fixture explicitly without changing the rejection or unchanged-state assertions. The full collision file passes 512 assertions across 32 testsets on each Julia version. Logs: `/tmp/hybrid-full-67f01b7-julia113.log`, `/tmp/hybrid-collision-fixture112.log`, `/tmp/hybrid-collision-fixture113.log`. A new integrated run is still required. |
| MOM-001 | Particle weights used as moment scratch are silently overwritten | Confirmed; canonical moment repair planned | VERIFIED on Julia 1.13.1: `momentum!`, `compute_moments!` and `compute_moments_multi!` change weight 0.4 to 9.6; both bulk-velocity APIs return (2, 6, 24) for input (2, 3, 4). Log: `/tmp/hybrid-serial-moment-alias-probe.log`. [Repair plan](MOMENT_STORAGE_PLAN.md) covers shared output/scratch preflight before mutation, independent views, serial/MPI callers and existing allocation gates. |

## Completeness roadmap

These are capability gaps or assessment tasks, not all classified as bugs.
Inspect implementation contracts before choosing the design; add an acceptance
benchmark before implementation. The roadmap remains open to new findings.

| Priority | Capability / gap | Plan and completion evidence | State |
| --- | --- | --- | --- |
| 1 | Reliable variable-step evolution | Correct phase state in ES/EM, hybrid/CAM-CL and MPI; preserve it in checkpoints; verify analytic trajectories and convergence | Active |
| 1 | Restart across production solvers | Hybrid, CAM-CL, electrostatic/electromagnetic PIC and Hall-MHD share capture/validation/restore infrastructure; exact continuation includes EM species and Hall retained stage history | Existing solver families and serial hybrid collections verified; boundary validation follow-ups active |
| 1 | Reproducible archives and bounded field I/O | Fix checksum representation and malformed-file allocation; cross-process and corruption tests | Active |
| 2 | Multispecies hybrid evolution | Collection initialization/stepping absence reproduced; [implementation plan](MULTISPECIES_HYBRID_PLAN.md) defines shared charge moments, independent species push, restart and conservation/dispersion acceptance | Serial core, collection restart and conservation diagnostics verified; MPI collection evolution verified; distributed collection restart follows |
| 2 | Multispecies conservation diagnostics | Canonical reductions and budgets support typed collections and empty populations; independent particle/field/closure sums and scalar regressions pass on both Julia versions | Verified |
| 2 | Exactly conservative coupled hybrid integration | Existing explicit predictor has second-order momentum error verified against a discrete Boris oracle. Assess a coupled midpoint field/current solve with convergence, conservation and nonlinear solve-residual acceptance | Assessment |
| 2 | Evolving electron pressure integrated into steppers | Standalone `advance_electron_pressure!` exists; built-in periodic steppers use algebraic closures. Define transport/heating coupling and energy-budget acceptance | Planned |
| 2 | Multidimensional EM electron subcycling | Common EM path extended; charge-continuity, Gauss-law, species phase, coupled convergence and Float32/64 checks pass | Verified |
| 2 | Nonperiodic EM-PIC boundaries | README declares periodic EM-PIC only; design compatible field/particle boundaries with reflection/transmission and charge tests | Planned |
| 3 | Production distributed stepping | MPI reference stepping replicates fields; design domain-owned fields/moments with migration and distributed operators; verify serial equivalence and scaling | Planned |
| 3 | Production GPU kernels | README claims array plumbing only; implement coupled gather/push/deposit/field paths with backend equivalence, memory and throughput measurements | Planned |
| 2 | Spatially local, interspecies Coulomb collisions | Current public Coulomb operator takes one ParticleSet; design cell-local multispecies pairing with momentum/energy and equilibration-rate oracles | Planned |
| 2 | Integrated implicit field advancement | Current `run_whistler` CN implementation evolves a linear reference packet; assess nonlinear implicit hybrid/EM integration with stability, accuracy and solve-residual benchmarks | Assessment |
| 3 | Model/geometry completeness | Assess Hall-MHD, gyrokinetic, AMR and fusion building blocks for missing end-to-end coupling before planning extensions | Assessment |
| 2 | Guiding-centre coordinate range and cancellation cost | Complete-position expression and exact cancellation repair verified. All measured ordinary hardware cases remain allocation-free; mixed-field latency increased in local paired measurements, so further profiling and optimization remain open | Correctness slice verified; performance follow-up planned |
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
