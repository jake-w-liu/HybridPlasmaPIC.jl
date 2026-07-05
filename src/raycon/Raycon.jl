# Raycon.jl — Julia port of the RAYCON ray-tracing / mode-conversion package.
#
# RAYCON (A. Jaun, KTH; A.N. Kaufman, LBL; E.R. Tracy, William & Mary; v7.0,
# 2006, with S. Richardson's modifications) computes linear mode conversion of
# RF waves in tokamak plasmas using ray tracing and ray splitting:
#
#   * E.R. Tracy, A.N. Kaufman, A. Jaun, Phys. Lett. A 290 (2001) 309.
#   * A. Jaun, E.R. Tracy, A.N. Kaufman, Plasma Phys. Control. Fusion 49 (2007) 43.
#   * E.R. Tracy, A.N. Kaufman, A. Jaun, Phys. Plasmas 14 (2007) 082102.
#
# This is a faithful port of /Users/jake/PlasmaWorkspace/raycon (MATLAB) — the
# mapping of every upstream .m file, all preserved upstream quirks, and every
# deliberate deviation are documented in
# docs/superpowers/specs/2026-07-05-raycon-port-notes.md. Reference data dumped
# from the original MATLAB code (tools/raycon_reference.m) pins the ported
# layers against the original in test/test_raycon.jl — at machine precision for
# the deterministic layers (equilibrium, geometry, dispersion, saddle,
# trajectory RHS) and at algorithm-tolerance level where upstream's own loose
# iterations limit reproducibility (conversion coefficients, ~1e-5).
#
# UNIFIED UNITS: the user-facing interface follows the package's Ω_ci
# normalization — every entry point has a method taking a `PlasmaUnits` first
# argument (raycon_normalized.jl) with lengths in d_i, wavenumbers in 1/d_i,
# frequencies in Ω_ci, B in B0, n in n0, T in m_i·v_A². Use those methods; the
# same functions without `PlasmaUnits` are the raw SI engine (meters, Tesla,
# rad/s, keV — the layer regression-pinned against the original MATLAB, which
# the normalized methods rescale exactly at the boundary). The traced state is
# z = (R, Z, kR, kZ) in the poloidal plane with the toroidal wavenumber k_φ
# held constant (upstream convention), evolved in the ray parameter σ (scale-
# invariant; the physical time direction is sign(∂U/∂ω)).
#
# Deliberate deviations from upstream (all behavior-preserving unless noted):
#   - no globals, no GUI/plots: explicit structs + programmatic driver;
#   - flux-surface root finding to 1e-12 (upstream fzero TolX 1e-4);
#   - our own Dormand–Prince 4(5) integrator with the upstream tolerances;
#     conversion/caustic monitors are recorded at accepted steps;
#   - cld3x3 mode-conversion analysis is implemented with corrected math
#     (upstream's 3x3 second-derivative expressions carry a known sign bug,
#     dispertok.m lines 476-477 "FIX THIS"): exact polynomial extraction of
#     the determinant derivatives plus Tracy–Kaufman near-null-subspace
#     coupling, reducing to the pinned 2x2 result when the electrostatic
#     branch decouples (RCN-014);
#   - the (upstream-disabled) amplitude-transport layer is completed rather
#     than skipped: integrate_ray_amplitude / trace_rays(amplitude=true)
#     evolve the Riccati focusing tensor, lnE² and eikonal phase with Maslov
#     caustic switching, cyclotron/Landau/TTMP damping with per-species
#     deposition, and conversion amplitude bookkeeping, verified against
#     symplectic tangent-map oracles (RCN-015).

module Raycon

using LinearAlgebra
using ..HybridPlasmaPIC: PlasmaUnits, alfven_speed, gyrofrequency, inertial_length

export RayconConstants,
    SolovevEquilibrium,
    RayconProblem,
    cmod_parameters,
    cmod_units,
    PlasmaUnits,
    rho_edge,
    solovev_flux,
    map_flux,
    flux_surface_mesh,
    magnetic_geometry,
    plasma_profiles,
    stix_elements,
    dispersion_U,
    trajectory_rhs,
    conversion_monitors,
    polarization,
    cyclotron_frequencies,
    dUdomega,
    msw_dispersion,
    adjust_to_dispersion,
    cgamma,
    analyze_conversion,
    is_valid,
    integrate_ray,
    integrate_ray_amplitude,
    antenna_focusing,
    launch_ray,
    trace_rays,
    RayconTrace,
    RayconConversion,
    AmplitudeTrace

include("raycon_types.jl")
include("raycon_solovev.jl")
include("raycon_magnetic.jl")
include("raycon_dispersion.jl")
include("raycon_conversion.jl")
include("raycon_trace.jl")
include("raycon_amplitude.jl")
include("raycon_driver.jl")
include("raycon_normalized.jl")

end # module Raycon
