# Parallel Backend Extension Validation

This case validates CPU backend, extension metadata, pencil decomposition, halo
exchange, and MPI serial-smoke contracts.

## What Is Tested

- Supported extension registry and extension dependency metadata.
- CPU particle backend storage, copy roundtrip, and memory-pressure arithmetic.
- Pencil decomposition rank coordinates, local bounds, owner lookup, and wrapper APIs.
- Local field halo exchange and MPI single-rank Cartesian setup.

## Reference

The references are exact API contracts, exact array copies, known pencil bounds,
and exact halo arrays for a constructed local layout.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 14_parallel_backend_extension_validation
```
