# Distributed FFT Evaluation

Evaluation date: 2026-06-26.

## Decision

Use `PencilFFTs.jl`/`PencilArrays.jl` as the first MPI-distributed FFT path for
larger fully periodic 3D HybridPlasmaPIC cases. Do not maintain a custom transpose
implementation unless focused benchmarks show a concrete gap that PencilFFTs
cannot cover.

## Evidence

- `PencilFFTs.jl` is a Julia package for fast Fourier transforms on
  MPI-distributed Julia arrays through `PencilArrays.jl`.
- Its documented 3D tutorial creates an MPI communicator, a `Pencil`, and a
  `PencilFFTPlan`; by default a 3D data set is distributed on a 2D MPI topology.
- The package supports distributed N-dimensional transforms, arbitrary
  combinations of FFT-related transforms along dimensions, in-place and
  out-of-place plans, and documented high scalability.
- The repository has a maintained release series and documented MPI scaling.
  Pin the exact release during implementation after rechecking upstream tags.

Primary sources:

- https://github.com/jipolanco/PencilFFTs.jl
- https://jipolanco.github.io/PencilFFTs.jl/dev/tutorial/

`DaggerFFT` is a relevant emerging Julia research direction for task-scheduled
distributed FFTs, including CPU and GPU results, but it should be tracked rather
than adopted as the baseline until the maintained package/API and HybridPlasmaPIC
integration path are verified.

Primary source:

- https://arxiv.org/html/2601.12209v1

## Integration Plan

1. Keep `SpectralOperators.jl` as the single-rank operator API and oracle.
2. Add an optional distributed FFT extension around `PencilFFTs.jl` instead of
   adding `PencilFFTs` as an unconditional core dependency.
3. Limit the first integration to fully periodic 3D Fourier operators. Shock
   configurations with nonperiodic shock-normal physics should continue to use
   slab/local finite-difference logic until a mixed distributed operator is
   designed and tested.
4. Validate the extension with real `mpiexec` tests:
   - forward/backward round-trip on 2, 4, and 8 ranks;
   - distributed spectral derivative parity against `SpectralOperators.deriv!`;
   - divergence/projection identities against analytic periodic fields;
   - plan allocation and scratch-buffer reuse checks.
5. Benchmark before optimizing:
   - compare PencilFFTs default communication against any custom transpose idea;
   - record strong and weak scaling separately from correctness gates;
   - keep GPU-distributed FFT adoption gated on CUDA/ROCm hardware tests.

## Non-Goals

- This evaluation does not implement pencil decomposition.
- This evaluation does not establish production scaling.
- This evaluation does not claim GPU-distributed FFT support for HybridPlasmaPIC.
