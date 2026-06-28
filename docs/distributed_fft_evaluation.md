# Distributed FFT Evaluation

Evaluation date: 2026-06-26.
Implementation update: 2026-06-28.

## Decision

Use `PencilFFTs.jl`/`PencilArrays.jl` as the first MPI-distributed FFT path for
larger fully periodic 3D HybridPlasmaPIC cases. Do not maintain a custom transpose
implementation unless focused benchmarks show a concrete gap that PencilFFTs
cannot cover.

## Implementation Status

Implemented as the optional `HybridPlasmaPICPencilFFTSExt` package extension.
`PencilFFTs` and `PencilArrays` are weak dependencies, so the core package does
not load them unless users explicitly import those packages. The exported API is:

- `distributed_fft_plan((nx, ny, nz); comm, transform, periodic)`
- `distributed_fft_input(plan)` / `distributed_fft_output(plan)`
- `distributed_fft_forward!(output, plan, input)`
- `distributed_fft_inverse!(input, plan, output)`
- `distributed_fft_roundtrip_error(plan, input; output=...)`

Verified locally on `MPI.COMM_SELF`: extension loading, input/output allocation,
validation failures for non-3D/nonperiodic requests, forward-transform parity
against `FFTW.fft` through `PencilArrays.gather`, inverse transform, and
round-trip error reduction.

## Evidence

- `PencilFFTs.jl` is a Julia package for fast Fourier transforms on
  MPI-distributed Julia arrays through `PencilArrays.jl`.
- Its documented 3D tutorial creates an MPI communicator, a `Pencil`, and a
  `PencilFFTPlan`; by default a 3D data set is distributed on a 2D MPI topology.
- The package supports distributed N-dimensional transforms, arbitrary
  combinations of FFT-related transforms along dimensions, in-place and
  out-of-place plans, and documented high scalability.
- The repository has a maintained release series and documented MPI scaling.
  The implemented compatibility is `PencilFFTs = 0.15` and
  `PencilArrays = 0.19`.

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
   adding `PencilFFTs` as an unconditional core dependency. **Done.**
3. Limit the first integration to fully periodic 3D Fourier operators. **Done.**
   Shock
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

- This evaluation does not establish production scaling.
- This evaluation does not claim GPU-distributed FFT support for HybridPlasmaPIC.
