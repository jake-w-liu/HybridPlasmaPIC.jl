# Metrics Load-Balance Validation

This case validates runtime metrics, mixed-grid residuals, sorting, load balance,
and memory arithmetic.

## What Is Tested

- Particle work in a uniform electric field.
- Mixed SBP/Fourier `div(curl(.))` residual convergence behavior.
- Balanced tile partition coverage and minimax capacity.
- Cell-index oracle and particle sorting order.
- Exact memory-byte accounting.

## Reference

The references are exact uniform-field work arithmetic, expected residual
monotonicity, and exact integer partition/cell/memory calculations.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 08_metrics_loadbalance_validation
```
