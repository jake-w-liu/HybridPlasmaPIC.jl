# Analytic Spectral Operators

This case validates the three-dimensional Fourier operator layer.

## What Is Tested

- Fourier derivative and Laplacian accuracy.
- Gradient, divergence, and curl vector identities.
- Divergence-free projection behavior and recovery of projected fields.

## Reference

The references are exact trigonometric fields with known derivatives and vector
identities. Metrics are relative L2 errors and residual norms.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 03_analytic_spectral_operators
```
