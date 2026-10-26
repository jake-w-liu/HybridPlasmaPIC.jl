# Production oblique Hall-MHD histories

This case advances the actual periodic Hall-MHD RK4 solver at angles of 30 and
60 degrees between the spatial wavevector and the uniform magnetic field. It
measures all density, velocity and magnetic Fourier components, then fits one
measured complex component to `A exp((gamma-i*omega)*t)`. Signed frequency and
numerical growth/damping are compared with an independent linear fluid oracle.

```sh
julia --project=validation validation/run_validation.jl --case 35_production_oblique_hall --no-plots --artifact-dir /path/to/new-run
julia --project=/path/to/plot-environment validation/35_production_oblique_hall/plot.jl --artifact-dir /path/to/new-run/35_production_oblique_hall/artifacts
```

The plot environment requires PlotlySupply 1.8.0, DelimitedFiles and TOML, as in
the existing validation plotters. Use a fresh artifact directory: this case
refuses to overwrite an existing receipt or raw history. The suite runner's
existing artifact-cleaning policy is unchanged; do not reuse evidence roots.

The 60-profile matrix uses a one-dimensional periodic grid with three vector
components and a tilted background field. There are two cold electromagnetic
branches and three warm branches continued from parallel whistler,
ion-cyclotron and acoustic seeds. Both frequency signs are tested. Temperature
pairs are `(Ti,Te)=(0,0)` and `(0.3,0.2)`, giving the actual isothermal total
pressure coefficient `Ti+Te`. Baseline n16/dt0.01, independent n32 grid
refinement and dt0.005 time refinement use mode 1, amplitude1e-5 and duration80.
Twelve additional histories check fourth-order time convergence for mode-2
whistlers with both signs and temperatures.

Branch labels identify continuous eigenvectors, not a re-sorted frequency at
each angle. Tracking starts from distinct closed-form parallel seeds and uses
fixed 0.1-degree steps with one-to-one overlap assignments. Degenerate seeds,
near-zero roots and ambiguous assignments raise errors. No cold acoustic zero
mode is presented as a propagating branch. The same nonzero-mode continuity
relation recovers density at cold and warm temperatures without dividing by
zero pressure.

The final three quarters of each measured history must contain at least two
periods. Relative signed frequency error must be at most 1e-3 and absolute rate
error at most 1e-4. Residual, window-sensitivity, full-vector history,
polarization, mass and divergence constraints apply as specified in the
[implementation plan](../../docs/OBLIQUE_HALL_VALIDATION_PLAN.md). Independent
refinement differences and both consecutive time orders are checked explicitly.
The oracle supplies initialization and comparison only; it is never passed to
the complex-history fitter. The standard summary has 748 metrics.

Each raw CSV contains actual time and real/imaginary coefficients of all seven
components, normalized by `2/n` to account for the conjugate mode in real fields.
Branch CSVs retain every continuation node; the seed's assignment gap is NaN
because seed identity is determined by closed-form frequencies, not overlap.
The fit CSV records signed rates, residuals, conditional OLS standard errors,
half-window sensitivity, exact fit indices, density bounds and complete-vector
errors. The order CSV retains all temporal errors. The receipt hashes source,
package environment and artifacts, and records initial eigenvectors, completed
profiles and failure causes. Successful raw samples are streamed so a thrown
evolution error retains the completed prefix. `running` is incomplete; a
completed receipt says `pass`, `fail` or `error`.

The OLS errors are conditional regression diagnostics, not confidence intervals
for correlated simulation noise. These ideal-fluid modes have zero physical
growth/damping; their measured gamma checks numerical behavior. This case does
not establish kinetic damping, oblique hybrid PIC, electron inertia/mass-ratio
scans, anisotropic/evolving pressure, multidimensional spatial convergence,
nonperiodic boundaries or GPU/distributed behavior. CAP-007B3 remains open for
those separate regimes.
