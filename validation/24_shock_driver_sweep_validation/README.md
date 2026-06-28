# Shock Driver Sweep Validation

This case validates the perpendicular-shock driver wrappers and low-cost shock
campaign contracts.

## What Is Tested

- Direct 1D `run_perp_shock` frozen-in and finite-diagnostic contracts.
- `reproduce_established_shock` oracle behavior.
- `ramp_width_scan` and `box_length_scan` driver contracts.
- Direct `run_perp_shock3d` finite diagnostics and one-step four-spacecraft traces.
- `perp_shock_sweep`, `mach_sweep`, convergence, 3D campaign, restart, and dimension-comparison wrappers.

## Reference

The references are flux-freezing, mass-conservation shock-speed arithmetic,
compression-band checks, convergence thresholds, and exact driver/restart
contracts. The 3D campaign checks use deliberately small smoke parameters.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 24_shock_driver_sweep_validation
```
