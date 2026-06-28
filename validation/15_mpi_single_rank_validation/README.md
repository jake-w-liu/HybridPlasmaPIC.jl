# MPI Single-Rank Validation

This case validates MPI.jl integration on `MPI.COMM_SELF`.

## What Is Tested

- MPI initialization, rank/size, and Cartesian communicator layout.
- Nested allreduce diagnostics.
- GPU-aware status query and host-staging copy-back contracts.
- MPI moment computation against the serial reference.
- Single-rank particle migration, halo exchange, and checkpoint roundtrip.

## Reference

The references are serial moment/deposition results, exact one-rank MPI
transport invariants, and exact checkpoint field/particle roundtrips.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 15_mpi_single_rank_validation
```
