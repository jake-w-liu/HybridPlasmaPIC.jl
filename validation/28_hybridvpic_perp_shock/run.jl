#!/usr/bin/env julia
#
# Case 28 — live external HYBRID-CODE shock comparison vs Hybrid-VPIC.
#
# Hybrid-VPIC (LANL, github.com/lanl/vpic-kokkos@hybridVPIC) is a different, independent
# kinetic-ion/fluid-electron hybrid PIC code. build_vpic.sh runs it for a PERPENDICULAR
# shock matched to ours (θ=90°, β=1, Te=Ti, drift Mach 4 → shock-frame M_A≈6) and dumps
# Bz(x). Since Bz/n is frozen-in, the Bz profile gives compression (Bz₂/Bz₁) and overshoot
# (Bz_max/Bz₂); the shock-front speed gives the shock-frame M_A.
#
# This is the nonlinear analog of case 27 (NHDS dispersion): a live external plasma code,
# run from source, cross-checking our shock. VPIC source/binary/output are GITIGNORED.
# This all-validation case treats a missing or unparseable VPIC profile as a failed
# validation, not a skip.

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC, Printf

# compression, overshoot, shock-frame M_A from a wall-frame Bz(x) profile (cells, hx=1 d_i;
# wall/downstream at low x, upstream at high x; taui = run time in Ω_ci⁻¹; Vd = drift Mach).
function _vpic_shock_metrics(x::Vector{Float64}, bz::Vector{Float64}; taui, Vd)
    Bz1 = sum(bz[x.>0.92*maximum(x)]) / count(x .> 0.92 * maximum(x))   # upstream
    half = (maximum(bz) + Bz1) / 2
    xf = 0.0
    for i = length(bz)-1:-1:2
        if (bz[i] - half) * (bz[i+1] - half) <= 0 && bz[i] != bz[i+1]
            xf = x[i]
            break
        end
    end
    dn = (x .> 10) .& (x .< xf - 8)                                    # downstream plateau
    Bz2 = any(dn) ? sum(bz[dn]) / count(dn) : NaN
    over = maximum(bz) / Bz2
    Vs = xf / taui                                                    # shock speed (vA=1)
    return (compression = Bz2 / Bz1, overshoot = over, xf = xf, MA = Vd + Vs)
end

function case_28_hybridvpic_perp_shock(artifact_dir::AbstractString)
    id = "28_hybridvpic_perp_shock"
    cat = "external_plasma_code"
    refkind = "external_plasma_code"
    ref = "Hybrid-VPIC (lanl/vpic-kokkos@hybridVPIC): perpendicular shock, β=1, drift Mach 4"
    prof = joinpath(@__DIR__, "vpic_build", "run", "bz_profile.txt")

    if !isfile(prof)
        try
            run(`bash $(joinpath(@__DIR__, "build_vpic.sh"))`)
        catch
        end
    end
    if !isfile(prof)
        artifact = joinpath(artifact_dir, "$(id).csv")
        rows = (("hybridvpic_profile_available_error", 1.0, 0.0, "absolute", 1.0, 0.0),)
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = cat,
            reference_kind = refkind,
            reference = ref,
            rows = rows,
            artifact = artifact,
            notes = "Hybrid-VPIC not built (needs git+cmake+MPI; macOS installs OpenMPI via brew); run validation/28_hybridvpic_perp_shock/build_vpic.sh first.",
        )
    end

    x = Float64[]
    bz = Float64[]
    for ln in eachline(prof)
        f = split(strip(ln))
        length(f) < 2 && continue
        push!(x, parse(Float64, f[1]))
        push!(bz, parse(Float64, f[2]))
    end
    if length(x) < 20
        artifact = joinpath(artifact_dir, "$(id).csv")
        rows =
            (("hybridvpic_profile_parse_error", length(x), 20.0, "margin", 20.0 - length(x), 0.0),)
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = cat,
            reference_kind = refkind,
            reference = ref,
            rows = rows,
            artifact = artifact,
            notes = "Hybrid-VPIC profile unparseable or too short.",
        )
    end

    v = _vpic_shock_metrics(x, bz; taui = 30.0, Vd = 4.0)
    # ours at the matched shock-frame M_A
    o = run_perp_shock_rh(;
        MA = v.MA,
        β = 1.0,
        N = 320,
        nppc = 96,
        nsteps = 1000,
        Lx = 200.0,
        t_avg_start = 8.0,
        seed = 1,
    )
    Xrh = rankine_hugoniot(MHDState(1.0, v.MA, 0.0, 1.0, 0.0, 1.0), 5 / 3).X

    over_band = max(0.0, 1.1 - v.overshoot, v.overshoot - 1.5)            # VPIC overshoot in Leroy's band
    kinetic_lt_fluid = max(0.0, v.compression - Xrh)                       # VPIC compresses ≤ fluid RH
    ov_agree = abs(o.overshoot - v.overshoot) / v.overshoot                # ours vs VPIC overshoot
    comp_agree = abs(o.compression - v.compression) / v.compression
    ma6_rel = abs(v.MA - 6.0) / 6.0
    profile_count_error = max(0.0, 20.0 - length(x))

    # VPIC profile (gitignored) for plotting the external shock structure
    overlay = joinpath(artifact_dir, "$(id)_vpic_profile.csv")
    _write_csv(overlay, ("x_di", "Bz_vpic"), ((x[i], bz[i]) for i = 1:length(x)))

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        ("vpic_overshoot_in_Leroy_band_1.1_1.5", v.overshoot, 1.26, "margin", over_band, 0.0),
        (
            "vpic_compression_le_fluidRH_kinetic",
            v.compression,
            Xrh,
            "margin",
            kinetic_lt_fluid,
            0.0,
        ),
        ("ours_vs_vpic_overshoot_rel_error", o.overshoot, v.overshoot, "relative", ov_agree, 0.25),
        (
            "ours_vs_vpic_compression_rel_error",
            o.compression,
            v.compression,
            "relative",
            comp_agree,
            0.15,
        ),
        ("hybridvpic_shockframe_MA6_rel_error", v.MA, 6.0, "relative", ma6_rel, 0.05),
        (
            "hybridvpic_profile_point_count_error",
            length(x),
            20.0,
            "margin",
            profile_count_error,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = cat,
        reference_kind = refkind,
        reference = ref,
        rows = rows,
        artifact = artifact,
        notes = "Hybrid-VPIC provides an independent nonlinear perpendicular-shock profile. The validation checks VPIC's shock regime, kinetic compression bound, and ours-vs-VPIC compression/overshoot agreement.",
    )
end

VALIDATION_CASE = ValidationCase(
    id = "28_hybridvpic_perp_shock",
    default = false,
    description = "Perpendicular shock vs the external Hybrid-VPIC code (live hybrid-to-hybrid comparison).",
    runner = case_28_hybridvpic_perp_shock,
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
