# Validation Suite

This folder contains reproducible validation cases for `HybridPlasmaPIC.jl`.
The default runner compares package output against analytical references and a
bundled published external-reference summary, then writes ignored CSV/PDF
artifacts for review.

Run the quick suite:

```sh
julia --project=validation validation/run_validation.jl --quick
```

Run every local case:

```sh
julia --project=validation validation/run_validation.jl --all
```

List cases:

```sh
julia --project=validation validation/run_validation.jl --list
```

Each validation case lives in its own folder:
`validation/<case_id>/README.md` explains what is tested and validated,
`validation/<case_id>/run.jl` contains the case runner,
`validation/<case_id>/plot.jl` contains the case plotter, and
`validation/<case_id>/artifacts/` contains that case's generated CSV/PDF files.
The top-level scripts only discover and orchestrate those case folders.

Suite summaries and plot metadata are written to `validation/artifacts/`. All
generated validation CSV/PDF files and artifact folders are intentionally
gitignored. The runner cleans generated CSV/PDF files for the selected case
folders and the summary folder at the start of each run so plots correspond to
the current case selection.

## Implemented Cases

| id | model path | comparison reference | default |
| --- | --- | --- | --- |
| `03_analytic_spectral_operators` | Fourier derivative, gradient, divergence, curl, projection, Laplacian | exact trigonometric derivatives and vector identities | quick |
| `05_particle_coupling_diagnostics` | particle deposition/gather, density, pressure tensor, velocity histogram | partition of unity, CIC linear reproduction, adjoint identity, exact centered pressure, weight conservation | quick |
| `04_particle_push_collision_boundaries` | Boris pusher, periodic/reflecting boundaries, BGK collisions | magnetic-speed invariant, exact boundary maps, momentum/energy conservation | quick |
| `09_normalization_migration_io` | unit normalization, logical rank layout/migration, checkpoint/restart | inverse unit maps, exact rank geometry, migration invariants, bitwise restart | quick |
| `22_shock_spacecraft_diagnostics` | shock-front, crossing, boundary-flux, frame-transform, and synthetic-spacecraft diagnostics | constructed shock/particle/spacecraft references | quick |
| `10_io_metadata_archive_validation` | metadata, portable field dumps, async save, run/archive containers, sampled particle dumps | exact roundtrip and checksum/oracle equality | quick |
| `08_metrics_loadbalance_validation` | particle work, mixed SBP/Fourier residual, sorting, load balance, memory estimate | uniform-field work identity, residual convergence, exact cell/load arithmetic | quick |
| `06_boundary_loading_kdv_smoothing` | open absorbing boundary, loading/weights, flux injection, binomial smoothing, KdV | exact compaction/cold flux, transfer function, analytic `sech^2` soliton | quick |
| `07_closure_budget_filter_diagnostics` | electron pressure evolution, energy/work budgets, pressure tensor/strain, power spectrum | exact manufactured fields and budget arithmetic | quick |
| `16_pic1d_hdf5_extension_validation` | `Electrostatic1D`, `EMPIC1D`, mobile/subcycled PIC, HDF5 extension | 1D Poisson exact field, discrete charge conservation, HDF5 dense-array roundtrip | quick |
| `14_parallel_backend_extension_validation` | extension metadata, CPU backend, pencil decomposition, local halo exchange, MPI single-rank smoke | exact contracts/ranges/halo arrays, MPI `COMM_WORLD` size-1 when loadable | quick |
| `13_threaded_backend_api_validation` | threaded deposition/density, CPU backend storage, memory telemetry, extension API fallbacks | serial deposition/density references, exact CPU particle roundtrip, extension contract behavior | quick |
| `15_mpi_single_rank_validation` | MPI.jl `COMM_SELF` reductions, staging, moments, migration, halos, checkpoint/restart | serial moment/reduction references and exact one-rank roundtrip invariants | quick |
| `23_shock_multidim_ramp_validation` | 2D/3D constructed shock surfaces, 3D `div B`, initial tanh ramp | constructed front positions/spectrum, uniform `div B`, analytic tanh field | quick |
| `11_api_contract_regression_validation` | exported API contracts, FD/spectral helpers, provenance IDs, RH solver, diagnostic reductions, restart/reference helpers | deterministic contracts, SBP identity, Fourier derivative, RH conservation residuals, exact roundtrips | quick |
| `12_model_runtime_smoke_validation` | Ohm/Faraday, moment construction, pusher equivalence, electrostatic PIC, Hall-MHD, 1D/2D/3D shock init/step | exact field identities, uniform equilibria, single/multi-species equality, finite zero-step shock contracts | quick |
| `01_analytic_poisson_2d` | 2D electrostatic spectral Poisson solve | exact Fourier-mode electric field | quick |
| `02_analytic_hall_mhd_continuity` | Hall-MHD continuity RHS | exact advected density wave RHS | quick |
| `21_published_preisser2020_summary` | bundled hybrid shock reference oracle | scalar summary from Preisser et al. Zenodo DOI `10.5281/zenodo.3697360` | quick |
| `17_distributed_fft_roundtrip` | PencilFFTs/PencilArrays extension | serial FFTW reference and inverse roundtrip | quick |
| `19_integrator_camcl_semiimplicit` | CAM-CL hybrid integrator and semi-implicit whistler integrator | ion-acoustic dispersion, hybrid-stepper frequency, CN/Euler amplification factors | all |
| `18_hybrid_ion_acoustic_dispersion` | hybrid PIC linear ion-acoustic wave | analytic `omega = k sqrt(Te)` | all |
| `20_empic_transverse_dispersion_2d` | 2D electromagnetic PIC transverse wave | analytic `omega^2 = omega_pe^2 + c^2 k^2` | all |
| `24_shock_driver_sweep_validation` | 1D shock sweep/convergence drivers and 3D campaign/dimension-comparison wrappers | flux-freezing, mass conservation, convergence threshold, and low-cost driver contract checks | all |

The quick cases are deterministic and fast. The dispersion cases use fixed
loading/seeding and tolerances above the measured FFT-bin and particle-noise
error. `validation/Project.toml` installs `HDF5`, `PencilArrays`, and
`PencilFFTs`, so the extension validations are expected to run and pass when the
suite is launched with `--project=validation`.

## Plotting

`run_validation.jl` calls `validation/plot_validation.jl` unless `--no-plots` is
passed. Plot generation intentionally runs in Julia's shared environment:

```sh
julia --project=@v#.# validation/plot_validation.jl --artifact-dir validation/artifacts
```

The plotter requires global `PlotlySupply` `1.8.0` and uses PlotlySupply's
high-level constructors without overriding the package template or Cartesian
axis styling. The current verified default template is `plotly_white`.

## External Open-Source Cross-Checks

The local suite does not vendor or launch heavyweight third-party solvers. The
CSV outputs are designed so external adapters can compare the same observables
against open-source tools when those tools are installed in a separate
environment:

- PlasmaPy: use its analytical plasma formula APIs as an independent Python-side
  check for simple dispersion/reference values.
- WarpX or Smilei: reproduce the cold electromagnetic PIC transverse mode and
  compare measured frequency against `20_empic_transverse_dispersion_2d`.
- Gkeyll: reproduce fluid or kinetic linear waves and compare frequency/growth
  observables against the corresponding CSV outputs.

Until those external solvers are installed and run, the executed comparisons are
the analytical references, FFTW serial reference, and the bundled published
hybrid-code summary.
