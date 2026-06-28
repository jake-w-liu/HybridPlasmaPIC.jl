# API Contract Regression Validation

This case validates deterministic contracts across exported helper APIs.

## What Is Tested

- Rank-seeded RNG and global particle-id contracts.
- SBP derivative identities and mixed SBP/Fourier derivatives.
- Spectral utilities, electron closures, gather/current/temperature helpers.
- Rankine-Hugoniot residuals, diagnostic reductions, load-balance helpers.
- MPI host-buffer contracts, 3D shock checkpoint roundtrip, and reference comparison.

## Reference

The references are exact algebraic contracts, manufactured derivative identities,
known conservation residuals, and exact roundtrip comparisons.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 11_api_contract_regression_validation
```
