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
gitignored. The runner cleans generated CSV/PDF files for every case folder
under the chosen artifact root and the summary folder at the start of each run
so outputs correspond to the current case selection.

## Implemented Cases

| id | model path | comparison reference | default |
| --- | --- | --- | --- |
| `01_analytic_poisson_2d` | 2D electrostatic spectral Poisson solve | exact Fourier-mode electric field | quick |
| `02_analytic_hall_mhd_continuity` | Hall-MHD continuity RHS | exact advected density wave RHS | quick |
| `03_analytic_spectral_operators` | Fourier derivative, gradient, divergence, curl, projection, Laplacian | exact trigonometric derivatives and vector identities | quick |
| `04_particle_push_collision_boundaries` | Boris pusher, periodic/reflecting boundaries, BGK collisions | magnetic-speed invariant, exact boundary maps, momentum/energy conservation | quick |
| `05_particle_coupling_diagnostics` | particle deposition/gather, density, pressure tensor, velocity histogram | partition of unity, CIC linear reproduction, adjoint identity, exact centered pressure, weight conservation | quick |
| `06_boundary_loading_kdv_smoothing` | open absorbing boundary, loading/weights, flux injection, binomial smoothing, KdV | exact compaction/cold flux, transfer function, analytic `sech^2` soliton | quick |
| `07_closure_budget_filter_diagnostics` | electron pressure evolution, energy/work budgets, pressure tensor/strain, power spectrum | exact manufactured fields and budget arithmetic | quick |
| `08_metrics_loadbalance_validation` | particle work, mixed SBP/Fourier residual, sorting, load balance, memory estimate | uniform-field work identity, residual convergence, exact cell/load arithmetic | quick |
| `09_normalization_migration_io` | unit normalization, logical rank layout/migration, checkpoint/restart | inverse unit maps, exact rank geometry, migration invariants, bitwise restart | quick |
| `10_io_metadata_archive_validation` | metadata, portable field dumps, async save, run/archive containers, sampled particle dumps | exact roundtrip and checksum/oracle equality | quick |
| `11_api_contract_regression_validation` | exported API contracts, FD/spectral helpers, provenance IDs, RH solver, diagnostic reductions, restart/reference helpers | deterministic contracts, SBP identity, Fourier derivative, RH conservation residuals, exact roundtrips | quick |
| `12_model_runtime_smoke_validation` | Ohm/Faraday, moment construction, pusher equivalence, electrostatic PIC, Hall-MHD, 1D/2D/3D shock init/step | exact field identities, uniform equilibria, single/multi-species equality, finite zero-step shock contracts | quick |
| `13_threaded_backend_api_validation` | threaded deposition/density, CPU backend storage, memory telemetry, extension API fallbacks | serial deposition/density references, exact CPU particle roundtrip, extension contract behavior | quick |
| `14_parallel_backend_extension_validation` | extension metadata, CPU backend, pencil decomposition, local halo exchange, MPI single-rank smoke | exact contracts/ranges/halo arrays, MPI `COMM_WORLD` size-1 when loadable | quick |
| `15_mpi_single_rank_validation` | MPI.jl `COMM_SELF` reductions, staging, moments, migration, halos, checkpoint/restart | serial moment/reduction references and exact one-rank roundtrip invariants | quick |
| `16_pic1d_hdf5_extension_validation` | `Electrostatic1D`, `EMPIC1D`, mobile/subcycled PIC, HDF5 extension | 1D Poisson exact field, discrete charge conservation, HDF5 dense-array roundtrip | quick |
| `17_distributed_fft_roundtrip` | PencilFFTs/PencilArrays extension | serial FFTW reference and inverse roundtrip | quick |
| `18_hybrid_ion_acoustic_dispersion` | hybrid PIC linear ion-acoustic wave | analytic `omega = k sqrt(Te)` | all |
| `19_integrator_camcl_semiimplicit` | CAM-CL hybrid integrator and semi-implicit whistler integrator | ion-acoustic dispersion, hybrid-stepper frequency, CN/Euler amplification factors | all |
| `20_empic_transverse_dispersion_2d` | 2D electromagnetic PIC transverse wave | analytic `omega^2 = omega_pe^2 + c^2 k^2` | all |
| `21_published_preisser2020_summary` | bundled hybrid shock reference oracle | scalar summary from Preisser et al. Zenodo DOI `10.5281/zenodo.3697360` | quick |
| `22_shock_spacecraft_diagnostics` | shock-front, crossing, boundary-flux, frame-transform, and synthetic-spacecraft diagnostics | constructed shock/particle/spacecraft references | quick |
| `23_shock_multidim_ramp_validation` | 2D/3D constructed shock surfaces, 3D `div B`, initial tanh ramp | constructed front positions/spectrum, uniform `div B`, analytic tanh field | quick |
| `24_shock_driver_sweep_validation` | 1D shock sweep/convergence drivers and 3D campaign/dimension-comparison wrappers | flux-freezing, mass conservation, convergence threshold, and low-cost driver contract checks | all |
| `25_leroy1982_perp_shock` | sustained two-state perpendicular shock | published overshoot, compression and reflected-fraction comparisons; pointwise limitations in case README | all |
| `26_hellinger2002_resolution` | perpendicular-shock resolution sweep | published maximum-gradient resolution trends | all |
| `27_nhds_dispersion_comparison` | independent Hall-MHD eigenvalue oracle | pinned NHDS kinetic dispersion; does not evolve the production solver | all |
| `28_hybridvpic_perp_shock` | production hybrid shock | pinned Hybrid-VPIC scalar shock comparisons | all |
| `29_firehose_instability` | anisotropic hybrid ions | unstable/stable threshold and mode-amplification controls | all |
| `30_ion_cyclotron_instability` | anisotropic hybrid ions | ion-cyclotron threshold and amplification controls | all |
| `31_weibel_instability` | counterstreaming electromagnetic PIC | field-energy growth versus a nonstreaming control | all |
| `32_magnetic_reconnection` | hybrid Harris sheet | coherent-mode growth versus uniform-field control; mechanism/rate validation remains open | all |
| `33_raycon_paper_validation` | ray propagation and mode conversion | exact Tracy identities and figure-derived Jaun targets | all |
| `34_production_parallel_dispersion` | evolved cold parallel hybrid/Hall waves | independent continuum/discrete equations, signed frequency/rate fits and refinement | all |

The default cases use deterministic inputs. Each non-default case documents
its own initialization, reference assumptions and fixed acceptance limits.
Case 34 uses production complex histories, signed frequency/rate fits and
independent refinement; conditional fit errors are separate from discretization
and correlated-particle uncertainty. `validation/Project.toml` installs `HDF5`,
`PencilArrays`, and
`PencilFFTs`, so the extension validations are expected to run and pass when the
suite is launched with `--project=validation`.

## Plotting

`run_validation.jl` calls `validation/plot_validation.jl` unless `--no-plots` is
passed. Plot generation intentionally runs in Julia's shared environment:

```sh
julia --project=@v#.# validation/plot_validation.jl --artifact-dir validation/artifacts
```

To plot a subset from existing artifacts, repeat `--case`:

```sh
julia --project=@v#.# validation/plot_validation.jl --artifact-dir validation/artifacts --case 01_analytic_poisson_2d
```

The plotter requires global `PlotlySupply` `1.8.0` and uses PlotlySupply's
high-level constructors without overriding the package template or Cartesian
axis styling. The current verified default template is `plotly_white`.

## External Open-Source Cross-Checks

Cases 27 and 28 build pinned NHDS and Hybrid-VPIC sources through the
[shared provenance owner](../docs/EXTERNAL_REFERENCE_PROVENANCE_PLAN.md).
Required toolchains must already be installed; missing/build-failed references
produce failure metrics. Cache reuse checks source, tools, inputs and outputs,
and failed build attempts remain available for inspection.

Additional comparisons remain planned; CSV outputs can support matched setups
in separately installed tools:

- PlasmaPy: use its analytical plasma formula APIs as an independent Python-side
  check for simple dispersion/reference values.
- WarpX or Smilei: reproduce the cold electromagnetic PIC transverse mode and
  compare measured frequency against `20_empic_transverse_dispersion_2d`.
- Gkeyll: reproduce fluid or kinetic linear waves and compare frequency/growth
  observables against the corresponding CSV outputs.

These proposed adapters do not constitute executed cross-code validation.
Case 27 compares eigenvalues; case 34 separately measures production solver
histories. Full-profile and matched cost/accuracy comparisons remain open in
the [completeness plan](../docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md).
