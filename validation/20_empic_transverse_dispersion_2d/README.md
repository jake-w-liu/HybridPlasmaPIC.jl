# EMPIC Transverse Dispersion 2D

This case validates the 2D electromagnetic PIC transverse-wave dispersion.

## What Is Tested

- 2D EMPIC transverse wave evolution.
- Frequency extraction from the simulated field history.
- Charge-conservation residual during the run.

## Reference

The analytic cold electromagnetic plasma reference is
`omega^2 = omega_pe^2 + c^2 k^2`. The case compares measured frequency to that
dispersion relation.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 20_empic_transverse_dispersion_2d
```
