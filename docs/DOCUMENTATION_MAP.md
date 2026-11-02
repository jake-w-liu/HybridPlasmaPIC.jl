# Documentation map

Start with the [README](../README.md) for supported models and entry points. The
[audit tracker](DEEP_DEBUG_TRACKER.md) owns current integration and exact-revision
acceptance; the [completeness plan](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)
owns missing capabilities and their gates. Source presence, a passing primitive
test and an end-to-end simulator are separate milestones.

The inventory below covers all 128 Markdown files tracked at audit baseline
`cd8b94d8a9737f68ed63fe0a64cb4c8d3196d70b`; this map is the additional audit entry
point. It is a dated inventory, not a claim that every historical experiment
was rerun. Later contracts and work orders are linked below as they land;
the current tracked set is available with `git ls-files '*.md'`.
Versioned plans retain baseline failures, candidate checklists and
older runtime results as evidence history. Read their named revision and current
tracker before treating a pending or completed gate as present-day acceptance.

At that baseline, local links and heading anchors were checked across the
inventory. External website availability, historical `/tmp` files and untracked generated artifacts
are separate evidence boundaries; commands and scientific results were not
rerun by this documentation-only audit.

## User and contributor entry points

- [CONTRIBUTING.md](../CONTRIBUTING.md)
- [README.md](../README.md)

## Current implementation trackers

- [docs/DEEP_DEBUG_TRACKER.md](DEEP_DEBUG_TRACKER.md)
- [docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md](PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md)

## Physics and API guides

- [docs/MIXED_HALL_OPERATORS.md](MIXED_HALL_OPERATORS.md) — real mixed operators,
  raw full-Hall induction, represented transpose, exact face-work ledger and
  constrained field/rate projector; full shock and runtime acceptance remain open.

- [docs/IMPLICIT_HYBRID.md](IMPLICIT_HYBRID.md) — synchronous joint velocity/B
  residual, bounded Newton–GMRES, represented work defects, distinct full
  checkpoint/archive and deferred runtime acceptance.

- [docs/COLLISIONAL_HYBRID.md](COLLISIONAL_HYBRID.md)
- [docs/COLLISIONAL_ELECTROSTATIC.md](COLLISIONAL_ELECTROSTATIC.md) — fixed-step
  electrostatic PIC/local-Coulomb source composition, sole clock and RNG restart.
- [docs/DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md](DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md)
  — actual-pencil cell-local ES/Coulomb composition, sole inner clock/source
  ledger and complete per-rank RNG/physical persistence; runtime acceptance open.
- [docs/DISTRIBUTED_COLLISIONAL_EMPIC.md](DISTRIBUTED_COLLISIONAL_EMPIC.md)
  — actual-pencil EM/Coulomb composition, carried field/current phase and full
  per-rank RNG/physical persistence; runtime acceptance open.
- [docs/COLLISIONAL_EMPIC.md](COLLISIONAL_EMPIC.md) — field-first electromagnetic
  PIC/local-Coulomb source composition, fixed normalization and full RNG restart.
- [docs/EVOLVING_ELECTRONS.md](EVOLVING_ELECTRONS.md) — stored pressure, transport
  and the bounded original-quotient conduction source candidate (policy 2).
- [docs/LOCAL_COULOMB.md](LOCAL_COULOMB.md)
- [docs/BOUNDED_PARTICLE_COUPLING_PLAN.md](BOUNDED_PARTICLE_COUPLING_PLAN.md) — static
  cloud primitives and the event-free longitudinal/1D3V electromagnetic source
  candidates, including original-geometry force completion and schema-2 policy.
- [docs/SOURCED_YEE_FIELDS.md](SOURCED_YEE_FIELDS.md)
- [docs/YEE_CONDUCTIVE_LAYER.md](YEE_CONDUCTIVE_LAYER.md)
  — graded scalar E/B conductivity, exact loss/work/field-defect budgets, fixed
  staggering and full nested restart; runtime acceptance remains open.
- [docs/YEE_PML.md](YEE_PML.md) — causal directional ADE memories, midpoint KDK,
  signed exchange/field-defect budgets and complete restart; reflection/stability
  and runtime acceptance remain open.
- [docs/MAXWELL_VACUUM_CHECKPOINT.md](MAXWELL_VACUUM_CHECKPOINT.md)
  — synchronized serial Yee vacuum state, matching-grid restore and archives.
- [docs/YEE_IMPACT_LONGITUDINAL.md](YEE_IMPACT_LONGITUDINAL.md)
  — one-cluster finite-cloud reflection, carried-error and restart source contract.
- [docs/YEE_ABSORBING_LONGITUDINAL.md](YEE_ABSORBING_LONGITUDINAL.md)
  — one isolated passive full-cloud capture, force-jump certificate and restart.
- [docs/DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md)
  — actual-pencil ES particle ownership, complete CIC transport, migration, phase
  and charge ledgers; runtime acceptance deferred.
- [docs/DISTRIBUTED_ELECTROSTATIC_CHECKPOINT.md](DISTRIBUTED_ELECTROSTATIC_CHECKPOINT.md)
  — matching actual-layout ES generations, complete particles/fields/phase and
  detached all-rank restore; runtime acceptance deferred.
- [docs/DISTRIBUTED_POISSON_FIELD.md](DISTRIBUTED_POISSON_FIELD.md)
  — static rank-owned electrostatic field, exact charge ledger and neutrality policy.
- [docs/DISTRIBUTED_MAXWELL_IMPLEMENTATION.md](DISTRIBUTED_MAXWELL_IMPLEMENTATION.md)
  — rank-owned periodic vacuum field source, scalar/curl boundaries and deferred
  runtime acceptance; separate ES/EM coupled source models are available.
- [docs/DISTRIBUTED_EMPIC_IMPLEMENTATION.md](DISTRIBUTED_EMPIC_IMPLEMENTATION.md)
  — actual-pencil all-mobile nonrelativistic CIC/Boris, corrected current,
  carried magnetic/velocity phase and scalar diagnostics; runtime acceptance
  and restart acceptance remain open.
- [docs/DISTRIBUTED_EMPIC_CHECKPOINT.md](DISTRIBUTED_EMPIC_CHECKPOINT.md)
  — complete eleven-field/current/endpoint phase, checksummed rank generations
  and detached stale-destination restore; runtime acceptance remains open.
- [docs/DISTRIBUTED_HYBRID_IMPLEMENTATION.md](DISTRIBUTED_HYBRID_IMPLEMENTATION.md)
  — actual-pencil isothermal ion collections, frozen-midpoint RK4, endpoint
  predictor and scalar diagnostics; runtime acceptance remains open.
- [docs/DISTRIBUTED_HYBRID_CHECKPOINT.md](DISTRIBUTED_HYBRID_CHECKPOINT.md)
  — complete eighteen-field/particle/phase generations, matching layout and
  detached stale-destination restore; runtime acceptance remains open.
- [docs/DISTRIBUTED_MAXWELL_CHECKPOINT.md](DISTRIBUTED_MAXWELL_CHECKPOINT.md)
  — matching-layout physical-pencil restart, checksums, publication and restore scope.

## Implementation, audit and design records

- [docs/ARCHIVE_INFERENCE_PLAN.md](ARCHIVE_INFERENCE_PLAN.md)
- [docs/BACKEND_MEMORY_CONSTRUCTOR_PLAN.md](BACKEND_MEMORY_CONSTRUCTOR_PLAN.md)
- [docs/BIGFLOAT_TIMESTEP_CONTEXT_PLAN.md](BIGFLOAT_TIMESTEP_CONTEXT_PLAN.md)
- [docs/BOUNDED_MAXWELL_IMPLEMENTATION.md](BOUNDED_MAXWELL_IMPLEMENTATION.md)
- [docs/BOUNDED_PARTICLE_COUPLING_PLAN.md](BOUNDED_PARTICLE_COUPLING_PLAN.md)
- [docs/BORIS_ACCURACY_IMPLEMENTATION.md](BORIS_ACCURACY_IMPLEMENTATION.md)
- [docs/CHARGE_DEPOSITION_CONSERVATION_PLAN.md](CHARGE_DEPOSITION_CONSERVATION_PLAN.md)
- [docs/CHECKPOINT_CONFIGURATION_PLAN.md](CHECKPOINT_CONFIGURATION_PLAN.md)
- [docs/CI_RELIABILITY_PLAN.md](CI_RELIABILITY_PLAN.md)
- [docs/COLLISION_COMPOSITION_RESTART_PLAN.md](COLLISION_COMPOSITION_RESTART_PLAN.md)
- [docs/CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md](CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md)
- [docs/COULOMB_FINITE_DRAW_PLAN.md](COULOMB_FINITE_DRAW_PLAN.md)
- [docs/COULOMB_NULL_ROTATION_PLAN.md](COULOMB_NULL_ROTATION_PLAN.md)
- [docs/COULOMB_SMALL_ANGLE_PLAN.md](COULOMB_SMALL_ANGLE_PLAN.md)
- [docs/COULOMB_WEIGHT_MEAN_PLAN.md](COULOMB_WEIGHT_MEAN_PLAN.md)
- [docs/DISTRIBUTED_DERIVATIVE_PLAN.md](DISTRIBUTED_DERIVATIVE_PLAN.md)
- [docs/DISTRIBUTED_FFT_COLLECTIVE_PLAN.md](DISTRIBUTED_FFT_COLLECTIVE_PLAN.md)
- [docs/DISTRIBUTED_FFT_DIAGNOSTIC_PLAN.md](DISTRIBUTED_FFT_DIAGNOSTIC_PLAN.md)
- [docs/DISTRIBUTED_FFT_PRECISION_PLAN.md](DISTRIBUTED_FFT_PRECISION_PLAN.md)
- [docs/DISTRIBUTED_INVERSE_HEADROOM_REPAIR.md](DISTRIBUTED_INVERSE_HEADROOM_REPAIR.md)
- [docs/DISTRIBUTED_METADATA_ADMISSION_REPAIR.md](DISTRIBUTED_METADATA_ADMISSION_REPAIR.md)
- [docs/DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md](DISTRIBUTED_OPERATOR_IMPLEMENTATION_PLAN.md)
- [docs/DISTRIBUTED_OPERATOR_TAKEOVER_AUDIT.md](DISTRIBUTED_OPERATOR_TAKEOVER_AUDIT.md)
- [docs/ELECTROSTATIC_SPECIES_IMPLEMENTATION.md](ELECTROSTATIC_SPECIES_IMPLEMENTATION.md)
- [docs/EMPIC_SPECIES_IMPLEMENTATION.md](EMPIC_SPECIES_IMPLEMENTATION.md)
  — arbitrary full-vector collections, opt-in generic D1–3 dividing schedules,
  phase/history, checkpoint migration and deferred acceptance.
- [docs/EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md](EMPIC1D_SPECIES_IMPLEMENTATION_PLAN.md)
- [docs/EMPIC_CHARGE_COMBINATION_PLAN.md](EMPIC_CHARGE_COMBINATION_PLAN.md)
- [docs/EMPIC_CURRENT_COMBINATION_PLAN.md](EMPIC_CURRENT_COMBINATION_PLAN.md)
- [docs/EMPIC_LEGACY_CURRENT_MODE_PLAN.md](EMPIC_LEGACY_CURRENT_MODE_PLAN.md)
- [docs/EMPIC_NYQUIST_RANGE_PLAN.md](EMPIC_NYQUIST_RANGE_PLAN.md)
- [docs/EMPIC_RANGE_PLAN.md](EMPIC_RANGE_PLAN.md)
- [docs/ENERGY_RANGE_PLAN.md](ENERGY_RANGE_PLAN.md)
- [docs/EXACT_BIGFLOAT_CONVERSION_PLAN.md](EXACT_BIGFLOAT_CONVERSION_PLAN.md)
- [docs/EXACT_IEEE_CONVERSION.md](EXACT_IEEE_CONVERSION.md)
- [docs/EXTERNAL_REFERENCE_PROVENANCE_PLAN.md](EXTERNAL_REFERENCE_PROVENANCE_PLAN.md)
- [docs/GATHERED_FIELD_STORAGE_PLAN.md](GATHERED_FIELD_STORAGE_PLAN.md)
- [docs/GATHERED_MIDPOINT_STORAGE_PLAN.md](GATHERED_MIDPOINT_STORAGE_PLAN.md)
- [docs/GYRO_RING_PLAN.md](GYRO_RING_PLAN.md)
- [docs/HYBRID_DIAGNOSTIC_RANGE_PLAN.md](HYBRID_DIAGNOSTIC_RANGE_PLAN.md)
  — current original-reservoir complete-total source and deferred acceptance,
  with rejected V6 cost history retained separately.
- [docs/HYBRID_ENERGY_MODEL_PLAN.md](HYBRID_ENERGY_MODEL_PLAN.md)
- [docs/IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md](IDEAL_FORCE_SOURCE_RANGE_WORK_ORDER.md)
- [docs/INTERVAL_VALIDATION_ALLOCATION_PLAN.md](INTERVAL_VALIDATION_ALLOCATION_PLAN.md)
- [docs/IONIZATION_ARITHMETIC_PLAN.md](IONIZATION_ARITHMETIC_PLAN.md)
- [docs/ISOTHERMAL_ENERGY_ACCURACY_PLAN.md](ISOTHERMAL_ENERGY_ACCURACY_PLAN.md)
- [docs/KINEMATICS_AUDIT_20261027.md](KINEMATICS_AUDIT_20261027.md)
- [docs/KINEMATICS_V5_INTEGRATION_20261027.md](KINEMATICS_V5_INTEGRATION_20261027.md)
- [docs/LINEAR_CN_RANGE_PLAN.md](LINEAR_CN_RANGE_PLAN.md)
- [docs/LOCAL_COULOMB_RATE_ALLOCATION_PLAN.md](LOCAL_COULOMB_RATE_ALLOCATION_PLAN.md)
- [docs/LOCAL_COULOMB_RATE_PLAN.md](LOCAL_COULOMB_RATE_PLAN.md)
- [docs/LOCAL_MULTISPECIES_COULOMB_PLAN.md](LOCAL_MULTISPECIES_COULOMB_PLAN.md)
- [docs/MAGNETIC_ENERGY_ACCURACY_PLAN.md](MAGNETIC_ENERGY_ACCURACY_PLAN.md)
- [docs/METAL_GATHERED_VALIDATION_PLAN.md](METAL_GATHERED_VALIDATION_PLAN.md)
- [docs/METAL_UNIFORM_PUSH_PLAN.md](METAL_UNIFORM_PUSH_PLAN.md)
- [docs/MIGRATION_STORAGE_PLAN.md](MIGRATION_STORAGE_PLAN.md)
- [docs/MODE_HISTORY_CENTERING_PLAN.md](MODE_HISTORY_CENTERING_PLAN.md)
- [docs/MOMENT_NORMALIZATION_PLAN.md](MOMENT_NORMALIZATION_PLAN.md)
- [docs/MOMENT_STORAGE_PLAN.md](MOMENT_STORAGE_PLAN.md)
- [docs/MPFR_CONCURRENCY_PLAN.md](MPFR_CONCURRENCY_PLAN.md)
- [docs/MPFR_UPSTREAM_REPORT.md](MPFR_UPSTREAM_REPORT.md)
- [docs/MPI_CHECKPOINT_PLAN.md](MPI_CHECKPOINT_PLAN.md)
- [docs/MULTISPECIES_HYBRID_PLAN.md](MULTISPECIES_HYBRID_PLAN.md)
- [docs/NGP_NEAREST_NODE_PLAN.md](NGP_NEAREST_NODE_PLAN.md)
- [docs/NORMALIZATION_RANGE_PLAN.md](NORMALIZATION_RANGE_PLAN.md)
- [docs/OBLIQUE_HALL_VALIDATION_PLAN.md](OBLIQUE_HALL_VALIDATION_PLAN.md)
- [docs/PARTICLE_CHARGE_MASS_RANGE_PLAN.md](PARTICLE_CHARGE_MASS_RANGE_PLAN.md)
- [docs/RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md](RELATIVISTIC_OUTPUT_PUBLICATION_PLAN.md)
- [docs/RELATIVISTIC_NORMAL_ACCURACY_PLAN.md](RELATIVISTIC_NORMAL_ACCURACY_PLAN.md)
- [docs/PARTICLE_APPEND_PLAN.md](PARTICLE_APPEND_PLAN.md)
- [docs/PARTICLE_KINEMATICS_RANGE_PLAN.md](PARTICLE_KINEMATICS_RANGE_PLAN.md)
- [docs/PARTICLE_PRIMING_INTERVALS.md](PARTICLE_PRIMING_INTERVALS.md)
- [docs/PERIODIC_MIGRATION_PLAN.md](PERIODIC_MIGRATION_PLAN.md)
- [docs/POLYTROPIC_ENERGY_RANGE_PLAN.md](POLYTROPIC_ENERGY_RANGE_PLAN.md)
- [docs/PRODUCTION_DISPERSION_PLAN.md](PRODUCTION_DISPERSION_PLAN.md)
- [docs/RADIATION_FORCE_SPAN_PLAN.md](RADIATION_FORCE_SPAN_PLAN.md)
- [docs/RADIATION_REACTION_PLAN.md](RADIATION_REACTION_PLAN.md)
- [docs/RANK_CLASSIFICATION_ALLOCATION_PLAN.md](RANK_CLASSIFICATION_ALLOCATION_PLAN.md)
- [docs/RELATIVISTIC_VALIDATION_PLAN.md](RELATIVISTIC_VALIDATION_PLAN.md)
- [docs/RESISTIVE_DIAGNOSTIC_RANGE_AUDIT.md](RESISTIVE_DIAGNOSTIC_RANGE_AUDIT.md)
- [docs/RK4_MIXED_STABILITY_PLAN.md](RK4_MIXED_STABILITY_PLAN.md)
- [docs/SERIAL_PUBLICATION_PLAN.md](SERIAL_PUBLICATION_PLAN.md)
- [docs/SCALAR_PRESSURE_RATE_RANGE_PLAN.md](SCALAR_PRESSURE_RATE_RANGE_PLAN.md)
- [docs/SHOCK_CN_SOLVE_PLAN.md](SHOCK_CN_SOLVE_PLAN.md)
- [docs/SHOCK_ENERGY_FLUX_RANGE_PLAN.md](SHOCK_ENERGY_FLUX_RANGE_PLAN.md)
- [docs/TEST_SUBPROCESS_ENVIRONMENT_PLAN.md](TEST_SUBPROCESS_ENVIRONMENT_PLAN.md)
- [docs/THREADED_DEPOSITION_PLAN.md](THREADED_DEPOSITION_PLAN.md)
- [docs/TRANSPORT_HYPERRESISTIVITY_PLAN.md](TRANSPORT_HYPERRESISTIVITY_PLAN.md)
- [docs/TRANSPORT_SOURCE_RANGE_PLAN.md](TRANSPORT_SOURCE_RANGE_PLAN.md)
- [docs/VALIDATION_FREQUENCY_PLAN.md](VALIDATION_FREQUENCY_PLAN.md)
- [docs/WHISTLER_ALLOCATION_CI_DIAGNOSIS.md](WHISTLER_ALLOCATION_CI_DIAGNOSIS.md)
- [docs/WHISTLER_BIGINT_COUNTER_PLAN.md](WHISTLER_BIGINT_COUNTER_PLAN.md)
- [docs/WHISTLER_PACKET_ALLOCATION_PLAN.md](WHISTLER_PACKET_ALLOCATION_PLAN.md)
- [docs/distributed_fft_evaluation.md](distributed_fft_evaluation.md)

## Historical design specifications

- [docs/superpowers/specs/2026-07-05-ray-tracing-design.md](superpowers/specs/2026-07-05-ray-tracing-design.md)
- [docs/superpowers/specs/2026-07-05-raycon-port-notes.md](superpowers/specs/2026-07-05-raycon-port-notes.md)

## Release history

- [CHANGELOG.md](../CHANGELOG.md)

## Validation catalog and cases

- [validation/01_analytic_poisson_2d/README.md](../validation/01_analytic_poisson_2d/README.md)
- [validation/02_analytic_hall_mhd_continuity/README.md](../validation/02_analytic_hall_mhd_continuity/README.md)
- [validation/03_analytic_spectral_operators/README.md](../validation/03_analytic_spectral_operators/README.md)
- [validation/04_particle_push_collision_boundaries/README.md](../validation/04_particle_push_collision_boundaries/README.md)
- [validation/05_particle_coupling_diagnostics/README.md](../validation/05_particle_coupling_diagnostics/README.md)
- [validation/06_boundary_loading_kdv_smoothing/README.md](../validation/06_boundary_loading_kdv_smoothing/README.md)
- [validation/07_closure_budget_filter_diagnostics/README.md](../validation/07_closure_budget_filter_diagnostics/README.md)
- [validation/08_metrics_loadbalance_validation/README.md](../validation/08_metrics_loadbalance_validation/README.md)
- [validation/09_normalization_migration_io/README.md](../validation/09_normalization_migration_io/README.md)
- [validation/10_io_metadata_archive_validation/README.md](../validation/10_io_metadata_archive_validation/README.md)
- [validation/11_api_contract_regression_validation/README.md](../validation/11_api_contract_regression_validation/README.md)
- [validation/12_model_runtime_smoke_validation/README.md](../validation/12_model_runtime_smoke_validation/README.md)
- [validation/13_threaded_backend_api_validation/README.md](../validation/13_threaded_backend_api_validation/README.md)
- [validation/14_parallel_backend_extension_validation/README.md](../validation/14_parallel_backend_extension_validation/README.md)
- [validation/15_mpi_single_rank_validation/README.md](../validation/15_mpi_single_rank_validation/README.md)
- [validation/16_pic1d_hdf5_extension_validation/README.md](../validation/16_pic1d_hdf5_extension_validation/README.md)
- [validation/17_distributed_fft_roundtrip/README.md](../validation/17_distributed_fft_roundtrip/README.md)
- [validation/18_hybrid_ion_acoustic_dispersion/README.md](../validation/18_hybrid_ion_acoustic_dispersion/README.md)
- [validation/19_integrator_camcl_semiimplicit/README.md](../validation/19_integrator_camcl_semiimplicit/README.md)
- [validation/20_empic_transverse_dispersion_2d/README.md](../validation/20_empic_transverse_dispersion_2d/README.md)
- [validation/21_published_preisser2020_summary/README.md](../validation/21_published_preisser2020_summary/README.md)
- [validation/22_shock_spacecraft_diagnostics/README.md](../validation/22_shock_spacecraft_diagnostics/README.md)
- [validation/23_shock_multidim_ramp_validation/README.md](../validation/23_shock_multidim_ramp_validation/README.md)
- [validation/24_shock_driver_sweep_validation/README.md](../validation/24_shock_driver_sweep_validation/README.md)
- [validation/25_leroy1982_perp_shock/README.md](../validation/25_leroy1982_perp_shock/README.md)
- [validation/26_hellinger2002_resolution/README.md](../validation/26_hellinger2002_resolution/README.md)
- [validation/27_nhds_dispersion_comparison/README.md](../validation/27_nhds_dispersion_comparison/README.md)
- [validation/28_hybridvpic_perp_shock/README.md](../validation/28_hybridvpic_perp_shock/README.md)
- [validation/29_firehose_instability/README.md](../validation/29_firehose_instability/README.md)
- [validation/30_ion_cyclotron_instability/README.md](../validation/30_ion_cyclotron_instability/README.md)
- [validation/31_weibel_instability/README.md](../validation/31_weibel_instability/README.md)
- [validation/32_magnetic_reconnection/README.md](../validation/32_magnetic_reconnection/README.md)
- [validation/33_raycon_paper_validation/README.md](../validation/33_raycon_paper_validation/README.md)
- [validation/34_production_parallel_dispersion/README.md](../validation/34_production_parallel_dispersion/README.md)
- [validation/35_production_oblique_hall/README.md](../validation/35_production_oblique_hall/README.md)
- [validation/README.md](../validation/README.md)

## Reference tool guide

- [tools/matlab_shadow/README.md](../tools/matlab_shadow/README.md)
