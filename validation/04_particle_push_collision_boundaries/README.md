# Particle Push Collision Boundaries

This case validates particle pushing, boundary maps, and BGK collision
conservation.

## What Is Tested

- Boris pusher speed conservation in a uniform magnetic field.
- Periodic and reflecting particle boundary behavior.
- BGK collision momentum conservation.
- BGK collision kinetic-energy conservation.

## Reference

The references are exact uniform-B speed invariance, exact boundary-coordinate
maps, and exact total momentum/energy before and after collision relaxation.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 04_particle_push_collision_boundaries
```
