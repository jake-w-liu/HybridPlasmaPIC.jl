# Normalization Migration IO

This case validates unit normalization, logical rank layout, particle migration,
particle append, and checkpoint/restart invariants.

## What Is Tested

- Unit conversion roundtrips and normalized scalar identities.
- Periodic rank wrapping, rank indexing, rank counts, and rank bounds.
- Logical particle migration count and ID preservation.
- `append_particles!` ID behavior.
- Checkpoint restart state equality.

## Reference

The references are exact inverse unit maps, exact logical-rank geometry, exact
particle ID arrays, and bitwise/equality restart comparisons.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 09_normalization_migration_io
```
