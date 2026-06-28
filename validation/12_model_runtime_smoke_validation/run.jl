#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_12_model_runtime_smoke_validation(artifact_dir::AbstractString)
    id = "12_model_runtime_smoke_validation"

    n = 64
    l = 2π
    g = FourierGrid((n,), (l,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    k = 3.0
    f = HybridFields{1,Float64}((n,))
    fill!(f.n, 1.0)
    f.B[3] .= cos.(k .* x)
    ohms_law!(f, HybridModel(IsothermalElectrons(0.0)), g)
    ohm_error = maximum(
        (
            maximum(abs, f.J[2] .- (k .* sin.(k .* x))),
            maximum(abs, f.E[1] .- (@. k * sin(k * x) * cos(k * x))),
            maximum(abs, f.E[2]),
            maximum(abs, f.E[3]),
        ),
    )

    E = (zeros(Float64, n), zeros(Float64, n), cos.(k .* x))
    dB = ntuple(_ -> zeros(Float64, n), 3)
    faraday_rhs!(dB, E, g)
    faraday_error = maximum(
        (
            maximum(abs, dB[1]),
            maximum(abs, dB[2] .- (.-k .* sin.(k .* x))),
            maximum(abs, dB[3]),
        ),
    )

    fproj = HybridFields{1,Float64}((n,))
    fproj.B[1] .= sin.(x)
    divbuf = zeros(Float64, n)
    div0 = magnetic_divergence!(divbuf, fproj, g)
    project_b!(fproj, g)
    div1 = magnetic_divergence!(divbuf, fproj, g)
    projection_error = div0 > 0 && div1 / div0 < 1e-12 ? 0.0 : 1.0

    ps = ParticleSet{1,Float64}(n)
    load_lattice_1d!(ps, 0.0, l)
    set_density_weight!(ps, 1.0, g)
    ps.v[1] .= 0.1
    ps.v[2] .= -0.2
    ps.v[3] .= 0.3
    fsingle = HybridFields{1,Float64}((n,))
    fmulti = HybridFields{1,Float64}((n,))
    compute_moments!(fsingle, ps, g, CIC(), 1e-6)
    compute_moments_multi!(fmulti, [ps], g, CIC(), 1e-6)
    moments_error = maximum(
        (
            maximum(abs, fsingle.n .- fmulti.n),
            maximum(abs, fsingle.ui[1] .- fmulti.ui[1]),
            maximum(abs, fsingle.ui[2] .- fmulti.ui[2]),
            maximum(abs, fsingle.ui[3] .- fmulti.ui[3]),
        ),
    )

    Bconst = (zeros(Float64, n), zeros(Float64, n), ones(Float64, n))
    polytropic = PolytropicElectrons(0.5, 1.0, 5 / 3)
    energy_error = maximum(
        (
            abs(magnetic_energy(Bconst, g) - 0.5 * prod(g.L)),
            abs(electron_internal_energy(fill(1.0, n), polytropic, g) - 0.5 / (5 / 3 - 1) * prod(g.L)),
            maximum(abs.(momentum_budget(ps, Bconst, g).particle .- total_momentum(ps))),
        ),
    )

    p_uniform = ParticleSet{2,Float64}(3; q = 1.5, m = 2.0)
    p_gathered = ParticleSet{2,Float64}(3; q = p_uniform.q, m = p_uniform.m)
    p_uniform.x[1] .= [0.1, 0.2, 0.3]
    p_uniform.x[2] .= [0.4, 0.5, 0.6]
    p_uniform.v[1] .= [0.2, -0.1, 0.4]
    p_uniform.v[2] .= [-0.4, 0.3, -0.2]
    p_uniform.v[3] .= [0.7, -0.6, 0.5]
    for d = 1:2
        copyto!(p_gathered.x[d], p_uniform.x[d])
    end
    for c = 1:3
        copyto!(p_gathered.v[c], p_uniform.v[c])
    end
    E0 = (0.2, -0.1, 0.05)
    B0 = (0.03, 0.04, 0.2)
    dt = 0.03
    Epart = ntuple(c -> fill(E0[c], nparticles(p_uniform)), 3)
    Bpart = ntuple(c -> fill(B0[c], nparticles(p_uniform)), 3)
    push_uniform!(p_uniform, E0, B0, dt)
    push_gathered!(p_gathered, Epart, Bpart, dt)
    gathered_push_error = maximum(
        (
            maximum(abs, p_gathered.x[1] .- p_uniform.x[1]),
            maximum(abs, p_gathered.x[2] .- p_uniform.x[2]),
            maximum(abs, p_gathered.v[1] .- p_uniform.v[1]),
            maximum(abs, p_gathered.v[2] .- p_uniform.v[2]),
            maximum(abs, p_gathered.v[3] .- p_uniform.v[3]),
        ),
    )
    bkick = boris_kick(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.1)
    boris_error = abs(sum(abs2, bkick) - 1.0)

    ge = FourierGrid((4, 4), (2π, 2π))
    electrons = ParticleSet{2,Float64}(16; q = -1.0, m = 1.0)
    load_lattice!(electrons, (0.0, 0.0), ge.L, (4, 4))
    set_density_weight!(electrons, 1.0, ge)
    es = ElectrostaticPIC(ge, nparticles(electrons); n0 = 1.0)
    init_espic!(es, electrons)
    step_espic!(es, electrons, 0.1)
    espic_uniform_error = maximum(
        (
            maximum(abs, es.ne .- 1.0),
            maximum(abs, es.E[1]),
            maximum(abs, es.E[2]),
            maximum(abs, es.E[3]),
        ),
    )

    hall = HallMHDState(ge, HallMHDModel(IsothermalElectrons(0.5); Ti = 0.2))
    fill!(hall.fields.n, 1.0)
    fill!(hall.fields.B[3], 1.0)
    hall_mhd_ohms_law!(hall)
    step_hall_mhd!(hall, 0.05)
    hall_uniform_error = maximum(
        (
            maximum(abs, hall.fields.n .- 1.0),
            maximum(abs, hall.fields.B[3] .- 1.0),
            maximum(abs, hall.fields.E[1]),
            abs(hall.time[] - 0.05),
            abs(hall.step[] - 1),
        ),
    )

    sh1 = PerpShock(8, 1.0)
    ps1 = ParticleSet{1,Float64}(2)
    ps1.x[1] .= [0.25, 0.75]
    ps1.v[1] .= [-0.1, -0.2]
    ps1.weight .= shock_density_weight(1.0, 1.0, nparticles(ps1))
    deposit_moments!(sh1, ps1)
    compute_E!(sh1)
    shock1_deposit_error = abs(sum(sh1.n .* sh1.s.H) - sum(ps1.weight))
    init_shock!(sh1, ps1)
    step_shock!(sh1, ps1, 0.0; NB = 1)
    shock1_error =
        all(isfinite, sh1.Bz) && all(isfinite, sh1.Ex) && nparticles(ps1) == 2 ? 0.0 : 1.0

    sh2 = PerpShock2D(8, 4, 1.0, 1.0)
    ps2 = ParticleSet{2,Float64}(2)
    ps2.x[1] .= [0.25, 0.75]
    ps2.x[2] .= [0.1, 0.8]
    ps2.v[1] .= [-0.1, -0.2]
    ps2.weight .= shock2d_density_weight(1.0, 1.0, 1.0, nparticles(ps2))
    deposit_moments2d!(sh2, ps2)
    shock2_deposit_error =
        abs(sum(sh2.n .* reshape(sh2.sbp.H, :, 1)) * sh2.dy - sum(ps2.weight))
    init_shock2d!(sh2, ps2)
    compute_E2d!(sh2)
    step_shock2d!(sh2, ps2, 0.0; NB = 1)
    shock2_error =
        all(isfinite, sh2.Bz) && all(isfinite, sh2.Ex) && nparticles(ps2) == 2 ? 0.0 : 1.0

    sh3 = PerpShock3D(6, 4, 4, 1.0, 1.0, 1.0)
    ps3 = ParticleSet{3,Float64}(2)
    ps3.x[1] .= [0.25, 0.75]
    ps3.x[2] .= [0.1, 0.8]
    ps3.x[3] .= [0.2, 0.6]
    ps3.v[1] .= [-0.1, -0.2]
    ps3.weight .= shock3d_density_weight(1.0, 1.0, 1.0, 1.0, nparticles(ps3))
    deposit_moments3d!(sh3, ps3)
    shock3_deposit_error =
        abs(sum(sh3.n .* reshape(sh3.sbp.H, :, 1, 1)) * sh3.dy * sh3.dz - sum(ps3.weight))
    init_shock3d!(sh3, ps3)
    compute_E3d!(sh3)
    step_shock3d!(sh3, ps3, 0.0; NB = 1)
    shock3_error =
        all(isfinite, sh3.B[3]) && all(isfinite, sh3.E[1]) && nparticles(ps3) == 2 ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "12_model_runtime_smoke_validation.csv")
    rows = (
        ("ohms_law_hall_term_max_abs_error", ohm_error, 0.0, "absolute", ohm_error, 1e-10),
        ("faraday_rhs_max_abs_error", faraday_error, 0.0, "absolute", faraday_error, 1e-10),
        ("project_b_divergence_contract_error", projection_error, 0.0, "absolute", projection_error, 0.0),
        ("single_multi_species_moment_max_abs_error", moments_error, 0.0, "absolute", moments_error, 1e-12),
        ("energy_momentum_budget_max_abs_error", energy_error, 0.0, "absolute", energy_error, 1e-12),
        ("push_gathered_uniform_equivalence_error", gathered_push_error, 0.0, "absolute", gathered_push_error, 1e-12),
        ("boris_kick_speed_invariant_error", boris_error, 0.0, "absolute", boris_error, 1e-12),
        ("electrostatic_pic_uniform_step_error", espic_uniform_error, 0.0, "absolute", espic_uniform_error, 1e-12),
        ("hall_mhd_uniform_step_error", hall_uniform_error, 0.0, "absolute", hall_uniform_error, 1e-12),
        ("shock1d_deposit_compute_e_mass_error", shock1_deposit_error, 0.0, "absolute", shock1_deposit_error, 1e-12),
        ("shock1d_init_step_finite_error", shock1_error, 0.0, "absolute", shock1_error, 0.0),
        ("shock2d_deposit_mass_error", shock2_deposit_error, 0.0, "absolute", shock2_deposit_error, 1e-12),
        ("shock2d_init_step_finite_error", shock2_error, 0.0, "absolute", shock2_error, 0.0),
        ("shock3d_deposit_mass_error", shock3_deposit_error, 0.0, "absolute", shock3_deposit_error, 1e-12),
        ("shock3d_init_step_finite_error", shock3_error, 0.0, "absolute", shock3_error, 0.0),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "model_runtime_smoke",
        reference_kind = "analytic_or_contract",
        reference = "exact Ohm/Faraday identities, moment equivalence, uniform PIC/Hall-MHD equilibria, and finite shock init/step contracts",
        rows = rows,
        artifact = artifact,
    )
end
VALIDATION_CASE = ValidationCase(
    id = "12_model_runtime_smoke_validation",
    default = true,
    description = "Ohm/Faraday identities, moment equivalence, uniform PIC/Hall-MHD steps, and shock init/step smokes.",
    runner = case_12_model_runtime_smoke_validation,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
