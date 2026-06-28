# Threaded Backend API Validation

This case validates threaded deposition/density helpers and backend extension
API contracts.

## What Is Tested

- Threaded scalar deposition against serial deposition for multiple dimensions and shapes.
- Threaded density deposition against serial density deposition.
- Conservation, offset-value, and zero-particle threaded edge cases.
- Binomial smoothing workspace equivalence.
- CPU backend storage/copy/memory telemetry and extension registry/fallback contracts.
- CPU backend preparation and abstract backend API type contracts.

## Reference

The references are serial deposition/density outputs, exact conservation sums,
exact CPU copy roundtrips, and deterministic extension API contracts.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 13_threaded_backend_api_validation
```
