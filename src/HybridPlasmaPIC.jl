"""
    HybridPlasmaPIC

Dimension-parametric (1D3V/2D3V/3D3V) hybrid particle-in-cell plasma solver:
kinetic ions, massless fluid electrons, spectral fields. Built on the standalone
SpectralOperators.jl package (no global mutable state; configuration is passed
explicitly). See the implementation checklist in the repository root.
"""
module HybridPlasmaPIC

using FFTW
using LinearAlgebra
using Random
import SpectralOperators
using SpectralOperators:
    FourierGrid,
    deriv!,
    deriv,
    gradient!,
    divergence!,
    curl!,
    laplacian!,
    project_divfree!,
    exp_filter!,
    dealias_two_thirds!,
    fft_friendly_size,
    with_fftw_wisdom,
    smoothing_transfer,
    BinomialSmoothWorkspace,
    binomial_smooth!,
    SBP1D,
    sbp_deriv!,
    sbp_deriv,
    sbp_deriv_x!,
    FourierDerivYWorkspace,
    fourier_deriv_y!

include("extensions.jl")
include("utils/validation.jl")

# §5.3 source tree. Grouped by subsystem directory; include order is the original
# verified dependency order. Spectral and mixed SBP/Fourier operators live in the
# sibling SpectralOperators.jl package and are imported above.

# --- Meshes ---
include("meshes/cartesian_mesh.jl")
include("meshes/local_finite_difference.jl")

# --- Particles ---
include("models/abstract_model.jl")
include("particles/particle_set.jl")
include("particles/loading.jl")
include("particles/boris.jl")
include("particles/sorting.jl")
include("particles/collisions.jl")
include("parallel/domain_decomposition.jl")
include("parallel/pencil_decomposition.jl")
include("particles/provenance.jl")
include("particles/migration.jl")

# --- Particle boundaries ---
include("boundaries/periodic.jl")
include("boundaries/reflecting_wall.jl")
include("boundaries/open_field_boundary.jl")
include("boundaries/sponge_layer.jl")

# --- Coupling: deposition / gather / moments ---
include("coupling/shapes.jl")
include("coupling/deposit.jl")
include("coupling/gather.jl")
include("coupling/moments.jl")

# --- Electrons + hybrid model + field state ---
include("electrons/closures.jl")
include("fields/field_state.jl")              # HybridFields
include("models/hybrid_pic.jl")               # HybridModel + compute_moments!
include("electrons/ohms_law.jl")
include("electrons/energy_equation.jl")

# --- Integrators ---
include("integrators/abstract_integrator.jl")
include("integrators/extrapolated_leapfrog.jl")
include("integrators/camcl.jl")
include("integrators/field_subcycling.jl")
include("integrators/semi_implicit.jl")

# --- Initial conditions: waves + shocks ---
include("initial_conditions/rankine_hugoniot_shock.jl")
include("initial_conditions/linear_waves.jl")
include("initial_conditions/reflecting_wall_shock.jl")
include("boundaries/particle_reservoir.jl")
include("initial_conditions/reflecting_wall_shock2d.jl")
include("initial_conditions/reflecting_wall_shock3d.jl")
include("initial_conditions/shock_ramp.jl")
include("initial_conditions/uniform_plasma.jl")
include("initial_conditions/instabilities.jl")

# --- Models: full PIC ---
include("models/electrostatic.jl")
include("models/full_pic.jl")
include("models/hall_mhd.jl")

# --- Diagnostics ---
include("diagnostics/pressure_tensor.jl")
include("diagnostics/phase_space.jl")
include("diagnostics/spectra.jl")
include("diagnostics/conservation.jl")        # energy/momentum budgets (needs kinetic_energy)
include("diagnostics/shock_surface.jl")       # shock diag (needs PerpShock2D/3D)
include("diagnostics/synthetic_spacecraft.jl")
include("diagnostics/particle_history.jl")

# --- Parallel ---
include("parallel/threads.jl")
include("parallel/gpu.jl")
include("parallel/mpi.jl")

# --- IO ---
include("io/checkpoint.jl")
include("io/metadata.jl")
include("io/hdf5_output.jl")
include("io/restart.jl")

# --- Verification ---
include("verification/normalization.jl")
include("verification/shock_sweep.jl")        # run_perp_shock (used by Oracles/Campaign)
include("verification/shock_campaign.jl")
include("verification/shock_convergence.jl")
include("verification/oracles.jl")
include("verification/metrics.jl")
include("verification/dispersion.jl")

export FourierGrid, deriv!, deriv, gradient!, divergence!, curl!, project_divfree!, laplacian!
export ParticleSet,
    nparticles,
    load_uniform!,
    load_lattice!,
    load_lattice_1d!,
    load_maxwellian!,
    load_quiet_velocities!,
    set_density_weight!,
    rank_seed,
    rank_rng,
    global_particle_id,
    assign_global_particle_ids!,
    boris_kick,
    push_uniform!,
    push_gathered!,
    apply_periodic!,
    apply_reflecting!,
    apply_absorbing!
export ShapeFunction,
    NGP,
    CIC,
    TSC,
    deposit_scalar!,
    gather_scalar!,
    gather_vector!,
    density!,
    momentum!,
    current!,
    pressure_tensor!,
    temperature_components
export ElectronClosure,
    IsothermalElectrons,
    PolytropicElectrons,
    electron_pressure!,
    closure_gamma,
    HybridModel,
    HybridFields,
    compute_moments!,
    ohms_law!,
    faraday_rhs!,
    magnetic_divergence!,
    project_b!,
    compute_moments_multi!,
    advance_electron_pressure!
export HybridStepper,
    init!, step!, kinetic_energy, magnetic_energy, electron_internal_energy, mode_amplitude
export MHDState, rankine_hugoniot
export save_checkpoint, load_checkpoint!
export kdv_soliton, kdv_solve
export SBP1D, sbp_deriv!, sbp_deriv, sbp_deriv_x!, FourierDerivYWorkspace, fourier_deriv_y!
export PerpShock, init_shock!, step_shock!, deposit_moments!, compute_E!, shock_density_weight
export Electrostatic1D, ElectrostaticPIC, init_espic!, step_espic!, poisson_E!, field_energy
export flux_speed, flux_per_density, inject_face_1d!
export total_momentum,
    electric_work,
    temperatures_par_perp,
    velocity_histogram,
    phase_space_histogram,
    power_spectrum,
    pressure_strain,
    shock_front
export cell_index, sort_particles!, particles_per_cell, memory_bytes
export LogicalRankLayout,
    nranks,
    rank_coords,
    rank_index,
    rank_bounds,
    rank_of_position,
    PencilDecomposition3D,
    pencil_decomposition,
    pencil_nranks,
    pencil_rank_coords,
    pencil_rank_index,
    pencil_bounds,
    pencil_local_size,
    pencil_owner,
    pencil_orientation_axes,
    GPUAwareMPIStatus,
    gpu_aware_mpi_status,
    mpi_buffer_uses_host_staging,
    host_staging_buffer,
    copy_from_host_staging!,
    MPIBufferPlan,
    prepare_mpi_buffer,
    finish_mpi_buffer!,
    ensure_mpi_initialized!,
    mpi_initialized,
    mpi_comm_size,
    mpi_comm_rank,
    mpi_dims_create,
    MPICartesianCommunicator,
    create_cartesian_communicator,
    free_mpi_communicator!,
    mpi_cartesian_neighbor,
    mpi_rank_layout_description,
    MPI_CHECKPOINT_SCHEMA_VERSION,
    save_mpi_checkpoint,
    load_mpi_checkpoint!,
    mpi_allreduce_diagnostics,
    mpi_compute_moments!,
    mpi_exchange_field_halos!,
    mpi_exchange_ghost_moments!,
    mpi_init!,
    mpi_step!,
    mpi_migrate_particles!,
    exchange_field_halos!,
    exchange_ghost_moments!,
    reduce_diagnostics,
    sum_diagnostics,
    min_diagnostics,
    max_diagnostics,
    append_particles!,
    migrate_particles!
export particle_array_backend,
    particle_storage_backend,
    BackendMemoryStatus,
    backend_memory_status,
    memory_pressure,
    reclaim_backend_memory!,
    copy_particles_to_backend,
    copy_particles_to_host,
    prepare_gpu_backend!
export binomial_smooth!, smoothing_transfer, BinomialSmoothWorkspace
export supported_extensions,
    extension_name,
    extension_dependency_name,
    extension_loaded,
    loaded_extensions,
    require_extension,
    extension_dependency_module,
    extension_device_array_type,
    disallow_scalar_indexing!,
    write_field_hdf5,
    read_field_hdf5
export gather_at,
    SyntheticProbe, sample!, advance!, shock_frame, dehoffmann_teller_velocity, classify_reflected
export RunMetadata, capture_metadata, CHECKPOINT_SCHEMA_VERSION, save_run, load_run
export PlasmaUnits, alfven_speed, gyrofrequency, inertial_length, to_SI, to_normalized
export rh_branches
export SpectralOperators, exp_filter!, dealias_two_thirds!, fft_friendly_size, with_fftw_wisdom
export CAMCLStepper, init_camcl!, step_camcl!
export EMPIC, EMPIC1D, init_empic!, step_empic!, em_field_energy, charge_conservation_residual
export HallMHDModel, HallMHDState, hall_mhd_ohms_law!, hall_mhd_rhs!, step_hall_mhd!
export run_perp_shock, perp_shock_sweep
export collide_bgk!
export deposit_scalar_threaded!, density_threaded!
export energy_budget, momentum_budget, jdotE_density, resistive_dissipation
export archive_run, load_archive, sample_particles, operators_match
export initial_ramp!, ramp_width_scan, box_length_scan
export PerpShock2D,
    init_shock2d!,
    step_shock2d!,
    deposit_moments2d!,
    compute_E2d!,
    shock2d_density_weight,
    shock_surface
export PerpShock3D,
    init_shock3d!,
    step_shock3d!,
    deposit_moments3d!,
    compute_E3d!,
    shock3d_density_weight,
    shock_surface3d,
    magnetic_divergence3d,
    run_perp_shock3d
export checkpoint_shock3d,
    restore_shock3d, production_3d_case, shock_campaign_3d, compare_dims_shock
export cn_multiplier, run_whistler, compare_integrators_whistler
export crossing_time, four_spacecraft_timing, four_spacecraft_traces
export load_imbalance,
    tile_loads,
    particle_load_imbalance,
    balanced_tile_ranges,
    balanced_tile_loads,
    particle_load_balance
export compare_to_reference,
    reproduce_established_shock,
    published_hybrid_reference_ids,
    published_hybrid_reference_metadata,
    published_hybrid_reference,
    compare_to_published_hybrid_reference
export boundary_energy_flux,
    shock_surface_spectrum,
    transverse_coherence,
    CrossingLogger,
    log_crossings!,
    crossing_count,
    energy_gain,
    boundary_reflection_fraction,
    normal_incidence_frame
export particle_work!, mixed_divcurl_residual
export write_field, read_field, async_save
export mach_sweep, convergence_study

end # module
