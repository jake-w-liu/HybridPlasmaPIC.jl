# Analytic Hall-MHD Continuity

This case validates the Hall-MHD continuity equation on a periodic Fourier grid.

## What Is Tested

- The continuity RHS produced by the Hall-MHD implementation.
- Spectral differentiation of an advected density wave.
- Consistency of the numerical RHS with the analytic `-u dot grad(n)` reference.

## Reference

The reference is an exact trigonometric density wave with constant advection
velocity. The measured metric is the maximum absolute RHS error.

## Outputs

The runner writes a metric CSV in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error to tolerance.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 02_analytic_hall_mhd_continuity
```
