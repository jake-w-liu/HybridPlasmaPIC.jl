# Analytic Poisson 2D

This case validates the two-dimensional spectral electrostatic Poisson solve.

## What Is Tested

- Fourier-space Poisson inversion on a 2D periodic grid.
- Electric-field reconstruction from the solved potential.
- Agreement with an exact single Fourier-mode field.

## Reference

The reference is the analytic electric field for a manufactured sinusoidal
charge-density mode. The primary metric is the maximum absolute field error.

## Outputs

The runner writes the sampled field comparison CSV in this folder's `artifacts/`
directory. The plotter writes a PlotlySupply PDF of the measured and expected
slice.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 01_analytic_poisson_2d
```
