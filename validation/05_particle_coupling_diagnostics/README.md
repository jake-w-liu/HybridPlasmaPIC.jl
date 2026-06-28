# Particle Coupling Diagnostics

This case validates particle deposition, gather, pressure, and histogram
diagnostics.

## What Is Tested

- Density deposition integral conservation.
- CIC gather reproduction of a linear field.
- Deposit/gather adjoint identity.
- Pressure tensor for exact centered velocities.
- Velocity histogram and phase-space histogram weight conservation.

## Reference

The references are partition-of-unity identities, exact linear interpolation,
an adjoint inner-product identity, and exact weighted sums.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 05_particle_coupling_diagnostics
```
