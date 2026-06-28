# Shock Spacecraft Diagnostics

This case validates shock diagnostics and synthetic-spacecraft geometry helpers.

## What Is Tested

- Shock-front position/width extraction.
- Particle crossing logger count and energy gain.
- Shock-frame and normal-incidence frame transforms.
- Boundary reflection fraction and boundary energy flux.
- Synthetic probe gather/sample/advance and crossing-time interpolation.
- Four-spacecraft timing, de Hoffmann-Teller velocity, and reflected-particle classification.

## Reference

The references are constructed shock profiles, exact particle trajectories,
analytic frame transforms, exact boundary flux arithmetic, and manufactured
spacecraft timing geometry.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 22_shock_spacecraft_diagnostics
```
