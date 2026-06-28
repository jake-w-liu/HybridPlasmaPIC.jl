# Hybrid Ion-Acoustic Dispersion

This case validates the hybrid PIC ion-acoustic wave frequency.

## What Is Tested

- Hybrid PIC linear ion-acoustic wave evolution.
- Frequency extraction from the simulated density/field signal.
- Agreement between measured frequency and analytic ion-acoustic theory.

## Reference

The analytic reference is `omega = k * sqrt(Te)` for the configured normalized
ion-acoustic setup.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 18_hybrid_ion_acoustic_dispersion
```
