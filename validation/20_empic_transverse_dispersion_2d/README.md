# EMPIC Transverse Dispersion 2D

This case validates the 2D electromagnetic PIC transverse-wave dispersion.

## What Is Tested

- 2D EMPIC transverse wave evolution.
- Frequency extraction from the simulated field history.
- Charge-conservation residual during the run.

The shared frequency estimator selects the strongest nonzero Fourier peak in
the complex field history, including both rotation directions, and returns its
unsigned angular frequency. It removes DC and fits a parabola to the peak and
its two neighbors. Constant, nonfinite or fewer-than-five-sample histories have
no measured frequency and fail validation. Sampling intervals must be finite
and positive. This finite-history estimate does not measure damping or separate
unresolved neighboring modes; its off-bin interpolation can be biased.

## Reference

The analytic cold electromagnetic plasma reference is
`omega^2 = omega_pe^2 + c^2 k^2`. The case compares measured frequency to that
dispersion relation.

## Outputs

The runner writes metric CSV rows in this folder's `artifacts/` directory. The
plotter writes a PlotlySupply PDF comparing error ratios to tolerances.

## Run

```sh
julia --project=validation validation/run_validation.jl --case 20_empic_transverse_dispersion_2d
```
