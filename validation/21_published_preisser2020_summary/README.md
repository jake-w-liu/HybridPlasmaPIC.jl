# Published Preisser 2020 Summary

This case validates the bundled published-reference summary oracle.

## What Is Tested

- Lookup of the packaged Preisser et al. hybrid-shock reference summary.
- Metadata and scalar-summary comparison through the reference-comparison helper.
- Max relative summary error against the bundled published oracle.

## Reference

The reference is the bundled scalar summary derived from the Zenodo dataset with
DOI `10.5281/zenodo.3697360`. The full external HDF5 source is not vendored.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 21_published_preisser2020_summary
```
