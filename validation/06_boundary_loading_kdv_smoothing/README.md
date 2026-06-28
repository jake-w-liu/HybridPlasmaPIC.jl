# Boundary Loading KdV Smoothing

This case validates boundary handling, particle loading helpers, smoothing, and
the KdV analytic soliton helper.

## What Is Tested

- Absorbing/open-boundary particle compaction.
- Quiet velocity loading and density-weight normalization.
- Cold flux speed, cold flux per density, and injection batch contracts.
- Binomial smoothing transfer function for a single mode.
- KdV soliton profile against the analytic `sech^2` solution.

## Reference

The references are exact particle-count/weight identities, analytic cold-flow
flux values, the binomial filter transfer function, and the analytic KdV soliton.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 06_boundary_loading_kdv_smoothing
```
