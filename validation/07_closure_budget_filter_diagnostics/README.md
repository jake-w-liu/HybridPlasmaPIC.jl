# Closure Budget Filter Diagnostics

This case validates electron closures, diagnostic budgets, pressure diagnostics,
and spectral diagnostics.

## What Is Tested

- Electron pressure evolution for manufactured density histories.
- Energy-budget component arithmetic and isothermal-budget NaN handling.
- `J dot E` density integration and resistive dissipation.
- Parallel/perpendicular temperature diagnostics and pressure-strain response.
- Power-spectrum peak detection for a manufactured mode.

## Reference

The references are exact manufactured fields and direct budget arithmetic. The
power-spectrum reference is the known input mode index.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 07_closure_budget_filter_diagnostics
```
