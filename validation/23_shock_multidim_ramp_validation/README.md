# Shock Multidim Ramp Validation

This case validates constructed 2D/3D shock surfaces, transverse coherence, 3D
magnetic divergence, and initial ramp setup.

## What Is Tested

- 2D shock-surface position extraction.
- 2D shock-surface spectrum peak and transverse autocorrelation.
- 3D shock-surface position extraction.
- Uniform 3D magnetic-divergence residual.
- Initial tanh ramp magnetic-field profile.

## Reference

The references are constructed step-front positions, the known transverse mode,
the exact transverse autocorrelation of that front, zero divergence for a uniform
field, and the analytic tanh ramp profile.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 23_shock_multidim_ramp_validation
```
