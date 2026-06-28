# Integrator CAM-CL Semiimplicit

This case validates CAM-CL and semi-implicit integrator behavior.

## What Is Tested

- CAM-CL ion-acoustic frequency against the analytic reference.
- CAM-CL frequency agreement with the existing hybrid stepper.
- Crank-Nicolson whistler amplification and energy behavior.
- Explicit Euler stiff-mode growth behavior.
- `compare_integrators_whistler` resolved and stiff-regime contracts.

## Reference

The references are the ion-acoustic dispersion relation and analytic linear
amplification properties of the CN/Euler whistler test equations.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 19_integrator_camcl_semiimplicit
```
