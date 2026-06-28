# Model Runtime Smoke Validation

This case validates finite, deterministic behavior across core model runtime
paths.

## What Is Tested

- Ohm-law Hall term and Faraday RHS identities.
- Divergence projection contract.
- Single-species and multi-species moment equivalence.
- Energy/momentum budget arithmetic and gathered-push equivalence.
- Electrostatic PIC and Hall-MHD uniform-step invariants.
- 1D/2D/3D shock deposition, electric-field computation, init, and zero-step paths.

## Reference

The references are exact manufactured field identities, uniform-equilibrium
invariants, and exact mass-conservation checks for shock moment deposition.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 12_model_runtime_smoke_validation
```
