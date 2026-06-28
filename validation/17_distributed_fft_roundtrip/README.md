# Distributed FFT Roundtrip

This case validates the PencilFFTs/PencilArrays package extension.

## What Is Tested

- Construction of a distributed 3D FFT plan on `MPI.COMM_SELF`.
- Forward distributed FFT output against a serial FFTW reference.
- Inverse distributed FFT roundtrip back to the original input.

## Reference

The external open-source reference is FFTW applied to the gathered serial array.
The extension path uses PencilFFTs and PencilArrays.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 17_distributed_fft_roundtrip
```
