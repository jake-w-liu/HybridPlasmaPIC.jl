#!/usr/bin/env julia
#
# SHK / Hellinger, Trávníček & Matsumoto (2002): the key quantitative result of
# 1-D perpendicular hybrid shock simulations (Geophys. Res. Lett. 29(24), 2234,
# doi:10.1029/2002GL015915; ref/Geophysical Research Letters - 2002 - Hellinger ...pdf):
#   Fig 1 — the maximum magnetic gradient max|dBy/dx|/B0 in the shock front vs the
#   grid spacing dx, for upstream β_p = 0.2, 0.5, 1.0 (β_e = 0.5, M_A ≈ 6.6):
#     • β_p = 1.0 (HOT): the gradient is RESOLUTION-INDEPENDENT (≈ const ~5) —
#       proton reflection stops the steepening (quasi-stationary).
#     • β_p = 0.2 (COLD): the gradient RISES as dx decreases (grid-determined,
#       nonstationary). Hellinger's headline: 1-D hybrid cannot describe the
#       nonstationary (cold) case.
# We reproduce the TREND (the physically meaningful claim): β_p=1.0 ~ flat,
# β_p=0.2 rising. (At dx≲0.08 with a fixed dt the whistler CFL is violated and the
# field blows up — recommended_dt avoids that; we use the resolvable dx=0.16,0.32.)

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC, Random

# time-averaged max |dBz/dx|/B0 in a sustained perpendicular shock at spacing dx
# and upstream ion thermal speed `vthi` (β_i = 2 vthi²); β_e = 0.5 fixed.
function _hellinger_maxgrad(dx::Float64, vthi::Float64; Lx = 80.0, U0 = 4.0, nppc = 150, seed = 1)
    T = Float64
    N = max(8, round(Int, Lx / dx))
    dx2 = Lx / N
    sh = PerpShock(N, T(Lx); Te = 0.25, γe = 5 / 3, η = 0.02, τ = U0, B0 = 1.0)  # β_e=0.5
    Np = nppc * N
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(seed)
    for p = 1:Np
        ps.x[1][p] = Lx * rand(rng)
        ps.v[1][p] = -U0 + vthi * randn(rng)
        ps.v[2][p] = vthi * randn(rng)
        ps.v[3][p] = vthi * randn(rng)
    end
    wp = shock_density_weight(1.0, T(Lx), Np)
    ps.weight .= wp
    init_shock!(sh, ps)
    inj = ShockInjector(
        MersenneTwister(seed + 1);
        n0 = 1.0,
        drift = U0,
        vthi = vthi,
        σt = vthi,
        weight = wp,
        first_id = Np + 1,
    )
    g = Float64[]
    for st = 1:900
        step_shock!(sh, ps, 0.02; NB = 2, injector = inj)
        if st * 0.02 >= 8 && st % 25 == 0
            xf, _ = shock_front(sh.Bz, sh.x)
            (isfinite(xf) && 8 < xf < Lx - 8) || continue
            push!(g, maximum(abs.(diff(sh.Bz))) / dx2)
        end
    end
    return isempty(g) ? NaN : sum(g) / length(g)
end

function case_26_hellinger2002_resolution(artifact_dir::AbstractString)
    id = "26_hellinger2002_resolution"
    vth_hot = sqrt(0.5)                # β_i = 1.0
    vth_cold = sqrt(0.1)               # β_i = 0.2
    g_hot_coarse = _hellinger_maxgrad(0.32, vth_hot)
    g_hot_fine = _hellinger_maxgrad(0.16, vth_hot)
    g_cold_coarse = _hellinger_maxgrad(0.32, vth_cold)
    g_cold_fine = _hellinger_maxgrad(0.16, vth_cold)

    # β_p=1.0: gradient ~ resolution-independent (relative change small)
    hot_rel_change = abs(g_hot_fine - g_hot_coarse) / g_hot_coarse
    # β_p=0.2: gradient RISES with refinement (fine clearly > coarse)
    cold_rises = (g_cold_fine > 1.3 * g_cold_coarse) ? 0.0 : 1.0
    # and the cold case steepens more strongly than the hot case (the discriminator)
    cold_steeper = (g_cold_fine / g_cold_coarse > g_hot_fine / g_hot_coarse) ? 0.0 : 1.0
    hot_magnitude_rel = abs(g_hot_fine - 5.0) / 5.0
    cold_magnitude_rel = abs(g_cold_fine - 13.0) / 13.0

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        (
            "beta1.0_resolution_independence_rel_change",
            hot_rel_change,
            0.0,
            "relative",
            hot_rel_change,
            0.5,
        ),
        ("beta0.2_gradient_rises_with_refinement", cold_rises, 0.0, "absolute", cold_rises, 0.0),
        ("cold_steepens_more_than_hot", cold_steeper, 0.0, "absolute", cold_steeper, 0.0),
        (
            "beta1.0_gradient_magnitude_vs_Hellinger_fig1_rel_error",
            g_hot_fine,
            5.0,
            "relative",
            hot_magnitude_rel,
            0.3,
        ),
        (
            "beta0.2_gradient_magnitude_vs_Hellinger_fig1_rel_error",
            g_cold_fine,
            13.0,
            "relative",
            cold_magnitude_rel,
            0.3,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "published_shock_benchmark",
        reference_kind = "literature_digitized",
        reference = "Hellinger et al. 2002 GRL 29(24):2234, doi:10.1029/2002GL015915 (Fig 1)",
        rows = rows,
        artifact = artifact,
        notes = "The case validates both Hellinger's qualitative result (hot shocks are nearly resolution-independent; cold shocks steepen with refinement) and the digitized order-of-magnitude Fig. 1 gradients at the resolvable dx values.",
    )
end

VALIDATION_CASE = ValidationCase(
    id = "26_hellinger2002_resolution",
    default = false,
    description = "Perpendicular-shock max-gradient resolution dependence vs Hellinger et al. 2002 (β_p=1.0 flat, β_p=0.2 rising).",
    runner = case_26_hellinger2002_resolution,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
