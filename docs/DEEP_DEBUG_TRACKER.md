# Plasma simulation audit and implementation tracker

Updated: 2026-10-13. Goal remains **active**.

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
  Final full-suite checks will use an immutable checkout. The 1.13.1
  log also reports SHA compatibility and extension-precompilation issues to assess.
- Full package suites are now running from isolated, fixed `1138317` checkouts
  on Julia 1.12.7 and 1.13.1; no final result yet. Logs:
  `/tmp/hybrid-full-1138317-julia112.log` and `/tmp/hybrid-full-1138317-julia113.log`.
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
| Electrostatic/electromagnetic PIC | Independent audit active; variable-step defect confirmed | Constant-force trajectory, Maxwell waves, charge continuity, mobile species |
| Hall-MHD, gyrokinetics, AMR | Inventory complete; detailed pass pending | Model-specific invariants, convergence, integration contracts |
| Collisions and radiation | Prior fixes exist; fresh pass pending | Conservation/rates, sampling, unequal weights, numerical robustness |
| Shocks/injection/campaigns | Inventory complete; detailed pass pending | RH laws, inflow flux, multidimensional/restart and published comparisons |
| Diagnostics and ray/Raycon solvers | Inventory complete; detailed pass pending | Independent references, units, transport/conversion invariants |
| Checkpoints/archives/field I/O | Independent audit active; defects confirmed | Configuration identity, exact bytes/checksums, atomic failure, restart |
| Threads/MPI/distributed FFT/GPU | Interfaces inspected; full pass pending | Real multirank tests, backend equivalence, ownership, hardware benchmarks |
| Packaging/quality/examples/validation | Immutable full-suite check pending; focused suites pass | CI-version tests, Aqua/JET findings, formatting, examples, validation cases |

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
| CFL-001 | Timestep recommendation lost representable results through intermediate overflow/underflow | Verified; committing | Valid 4-cell grid of length 1e-80 returned zero instead of the 512-bit oracle's 1.4184965709927286e-162. Logarithmic evaluation passes 3,802 helper/range/coupled-stability assertions on each of Julia 1.12.7 and 1.13.1 (`test/test_phase_helpers.jl`) |

## Completeness roadmap

These are capability gaps or assessment tasks, not all classified as bugs.
Inspect implementation contracts before choosing the design; add an acceptance
benchmark before implementation. The roadmap remains open to new findings.

| Priority | Capability / gap | Plan and completion evidence | State |
| --- | --- | --- | --- |
| 1 | Reliable variable-step evolution | Correct phase state in ES/EM, hybrid/CAM-CL and MPI; preserve it in checkpoints; verify analytic trajectories and convergence | Active |
| 1 | Restart across production solvers | CAM-CL implemented through the shared checkpoint path, preserving Blead/Blag/Bout/cl_h/cl_live. Electrostatic restart is next; assess other models afterward | CAM-CL verified; electrostatic design active |
| 1 | Reproducible archives and bounded field I/O | Fix checksum representation and malformed-file allocation; cross-process and corruption tests | Active |
| 2 | Multispecies hybrid evolution | `compute_moments_multi!` exists; built-in hybrid steppers document one proton species. Design shared-field multispecies advancement and conservation/dispersion benchmarks | Planned |
| 2 | Evolving electron pressure integrated into steppers | Standalone `advance_electron_pressure!` exists; built-in periodic steppers use algebraic closures. Define transport/heating coupling and energy-budget acceptance | Planned |
| 2 | Multidimensional EM electron subcycling | Extend the common EM path with charge-continuity, Gauss-law, species phase, and coupled convergence checks | Implementation active |
| 2 | Nonperiodic EM-PIC boundaries | README declares periodic EM-PIC only; design compatible field/particle boundaries with reflection/transmission and charge tests | Planned |
| 3 | Production distributed stepping | MPI reference stepping replicates fields; design domain-owned fields/moments with migration and distributed operators; verify serial equivalence and scaling | Planned |
| 3 | Production GPU kernels | README claims array plumbing only; implement coupled gather/push/deposit/field paths with backend equivalence, memory and throughput measurements | Planned |
| 2 | Spatially local, interspecies Coulomb collisions | Current public Coulomb operator takes one ParticleSet; design cell-local multispecies pairing with momentum/energy and equilibration-rate oracles | Planned |
| 2 | Integrated implicit field advancement | Current `run_whistler` CN implementation evolves a linear reference packet; assess nonlinear implicit hybrid/EM integration with stability, accuracy and solve-residual benchmarks | Assessment |
| 3 | Model/geometry completeness | Assess Hall-MHD, gyrokinetic, AMR and fusion building blocks for missing end-to-end coupling before planning extensions | Assessment |
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
