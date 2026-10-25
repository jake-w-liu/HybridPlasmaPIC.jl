# Production parallel-wave dispersion

This case evolves cold circularly polarized waves through the production
`HallMHDState` and `HybridStepper` solvers. It fits the measured complex Fourier
coefficient to `A exp((gamma - i*omega)*t)`, retaining signed frequency, growth
or damping, residuals, conditional slope errors and window sensitivity.

Run from the repository root after instantiating the validation environment:

```sh
julia --project=validation validation/run_validation.jl --case 34_production_parallel_dispersion --no-plots
julia --project=@v#.# validation/34_production_parallel_dispersion/plot.jl
```

The plotter uses the repository's required global PlotlySupply 1.8.0 environment.
If your shared environment uses another PlotlySupply version, the standalone
plotter also accepts Julia's `--project=/path/to/environment` pointing to an
isolated environment with PlotlySupply 1.8.0, DelimitedFiles and TOML installed.
Both commands accept `--artifact-dir`; the suite runner places case artifacts
under `<root>/34_production_parallel_dispersion/artifacts/`.

## Matrix and references

There are 40 periodic CPU Float64 profiles, all with density and background
magnetic field equal to one, zero ion/electron temperature, domain length
`2pi`, amplitude `1e-4` and duration 40. Both the whistler and ion-cyclotron
branches run at modes 1 and 2. Hall resistivity is 0 or 0.03; hybrid resistivity
is zero. Baseline resolution is 16 cells, timestep 0.01 and 16 particles/cell.
Separate profiles double the grid, halve the timestep or, for hybrid, double
the particle count. Hybrid runs use eight magnetic subcycles.

The independent cold continuum equations supply the reference eigenvalues.
The hybrid initializer uses the analytical CIC interpolation/deposition matrix
to select a single discrete eigenmode; the comparison remains against the
continuum reference. Neither reference enters the history fit. With the stored
helical coefficient `B_y + i B_z` at positive spatial mode, the whistler has
negative signed `omega` and the ion-cyclotron branch has positive signed `omega`.

Each fit uses the final three quarters of its raw history and must satisfy:

- Relative frequency error at most 3%; absolute growth/damping error at most 0.01.
- Phase and log-amplitude residual RMS at most 0.03.
- Absolute differences between first-half and second-half fitted frequency/rate
  at most 0.01, and at least two fitted periods.
- Maximum density departure from one at most `1e-6`.

Six further histories measure temporal order using the mode-2 whistler,
16 cells, amplitude `1e-6`, duration 1 and timesteps 0.04, 0.02 and 0.01.
The maximum complex-history error is compared against the independently
initialized linear eigenmode. Both successive orders must be at least 3.5 for
Hall RK4 or 1.5 for hybrid. Limits are fixed before measurement; no measured
failure is removed as a presumed noise floor.

## Artifacts and interpretation

- Every evolved history is saved as time, real/imaginary coefficient and magnitude.
- `*_fits.csv` contains all configurations, reference values, measured rates,
  residuals, conditional OLS standard errors and half-window differences.
  `duration` is the complete run; `fit_duration` is the fitted window.
- `*_orders.csv` records all six temporal error measurements.
- `*_provenance.toml` hashes the loaded package source, validation owners,
  active project/manifest and generated data. It records runtime and thread
  counts. `running` means incomplete; only a completed run records `pass` or `fail`.
- The standard summary contains 324 individual acceptance metrics.
- Plots compare all fitted rates with their references, show both temporal
  error curves and display one representative hybrid ion-cyclotron history.
  The plot-label CSV maps profile indices back to complete configurations.

OLS standard errors assume independent residuals of constant variance. They
are conditional regression diagnostics, not confidence intervals for correlated
PIC noise or discretization error. The deterministic cold lattice does not
replace independent-seed uncertainty studies. This case does not validate
oblique propagation, kinetic damping, finite temperature, mass-ratio scans,
electromagnetic PIC, nonperiodic boundaries or GPU/distributed evolution.
Those remain in the [production dispersion plan](../../docs/PRODUCTION_DISPERSION_PLAN.md).
