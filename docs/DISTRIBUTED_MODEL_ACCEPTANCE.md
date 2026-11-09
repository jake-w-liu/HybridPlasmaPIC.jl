# Physical-pencil model acceptance

The actual rank-owned models have a separate lifecycle gate from the replicated
MPI reference and distributed FFT/operator/metadata tests. Source review and
registration are complete. The first latest-Julia2-rank run at `38d1cd1` passes
3,728/3,725 charged ES/EM/hybrid assertions per rank. Maxwell's mode controls run,
then its continuation helper fails on `copy(::Tuple{})` for the empty particle
collection ([job log](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/42556458252/job/128963184546)).
The helper now collects population references for both tuples and vectors;
complete minimum/latest acceptance requires a successful rerun. No convergence,
scaling, compiler-resource or device result follows from these controls.

`test/test_distributed_models_core.jl` constructs real charged electrostatic,
electromagnetic and isothermal hybrid models on physical PencilArrays, plus
nonconstant vacuum Maxwell fields. Its dyadic load starts all 64 markers on
rank0, with empty peers. Independent CIC density and ID/tag/weight inventories,
nonzero fields/acceleration, actual owner migration, all physical arrays and
phase/cache state are checked. A separate analytic single-mode KDK reference
checks nonconstant Maxwell evolution. One-rank invalid-duration rejection must
complete collectively and preserve public physical/phase state.

Every family saves an advanced boundary, advances a different positive duration,
loads over stale fields/clocks/counts, and repeats continuation. All 4 ES, 11 EM,
18 hybrid and 6 Maxwell arrays are compared with the saved/continued values.
This checks matching-layout carried-state restart; other layouts, closures,
precisions, shapes, backend combinations and scientific refinement remain open.

The main workflow runs the file with 2 actual ranks on minimum and latest Julia.
The minimum job reuses its existing distributed environment; a separate latest
job resolves the current package in its own runner environment. Existing
operator and metadata failures remain failures, with separate retained logs.
This MPI-only file is intentionally outside single-process test groups.

The canonical launcher prepares its external package environment and runs it:

```sh
julia --project=. tools/ci_distributed_fft.jl prepare "$RUNNER_TEMP/distributed-env" 2
julia --project=. tools/ci_distributed_fft.jl models "$RUNNER_TEMP/distributed-env" 2
```

Use GitHub Actions when the development machine is busy. Inspect the exact
commit, every required job and retained test log before accepting a revision.
The source contracts remain in the [ES](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md),
[EM](DISTRIBUTED_EMPIC_IMPLEMENTATION.md), [hybrid](DISTRIBUTED_HYBRID_IMPLEMENTATION.md)
and [Maxwell](DISTRIBUTED_MAXWELL_IMPLEMENTATION.md) guides.
