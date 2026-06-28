# IO Metadata Archive Validation

This case validates metadata capture, field IO, async save, archive, and sample
dump helpers.

## What Is Tested

- Raw field write/read roundtrip.
- Async save snapshot behavior.
- `RunMetadata`, `save_run`, and `load_run` contracts.
- `archive_run` and `load_archive` state recovery.
- Particle sampling stride and operator-match checks.

## Reference

The references are exact serialization/raw-array roundtrips, preserved metadata
fields, and exact sampled particle indices.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 10_io_metadata_archive_validation
```
