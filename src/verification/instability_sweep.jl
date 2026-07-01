# instability_sweep.jl — driven kinetic-instability runners (Phase-1 physics that
# reuses the hybrid PIC engine: no new solver, only an initial-condition + a growth
# diagnostic). Analogous to shock_sweep.jl for shocks.

"""
    firehose_growth(; vth_par, vth_perp, β_e=1.0, N=128, Lx=25.6, nppc=200,
                      η=0.001, dt=0.02, nsteps=1600, NB=2, seed=1)

Drive the **parallel firehose instability** in the hybrid model and return the
transverse magnetic-fluctuation energy it produces. A uniform background field
`B₀ = x̂` is seeded with an anisotropic bi-Maxwellian ion distribution
(`vth_par` along B₀, `vth_perp` across); the fluid electrons are isotropic, so the
firehose threshold is set by the ions alone:

    unstable  ⇔  β_∥ − β_⊥ > 2  ⇔  T_∥ − T_⊥ > B₀²  ⇔  vth_par² − vth_perp² > 1
    (n₀ = 1, B₀ = 1, μ₀ = 1)

Above threshold the transverse field `δB_⊥ = (B_y, B_z)` grows exponentially; below
it, it stays at the particle-noise floor. Returns `(; wperp_max, wperp_final, wB0,
ratio_max, anisotropy, unstable_theory, nsamples)` where `ratio_max = wperp_max/wB0`
is the peak transverse energy as a fraction of the background magnetic energy.
"""
function firehose_growth(;
    vth_par::Real,
    vth_perp::Real,
    β_e::Real = 1.0,
    N::Integer = 128,
    Lx::Real = 25.6,
    nppc::Integer = 200,
    η::Real = 0.001,
    dt::Real = 0.02,
    nsteps::Integer = 1600,
    NB::Integer = 2,
    seed::Integer = 1,
)
    T = Float64
    N >= 8 || throw(ArgumentError("N must be ≥ 8"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    NB >= 1 || throw(ArgumentError("NB must be ≥ 1"))
    vpar = _require_finite_positive_real("vth_par", vth_par, T)
    vperp = _require_finite_positive_real("vth_perp", vth_perp, T)
    LxT = _require_finite_positive_real("Lx", Lx, T)
    βeT = _require_finite_nonnegative_real("β_e", β_e, T)
    dtT = _require_finite_positive_real("dt", dt, T)

    g = FourierGrid((Int(N),), (LxT,))
    dx = g.dx[1]
    # whistler CFL on the RK4 B-subcycle (d_i = 1): reject an unstable dt up front
    Kmax = T(π) / dx
    ω_w = T(0.5) * (sqrt(Kmax^4 + 4Kmax^2) + Kmax^2)
    dt_cfl = T(0.9) * Int(NB) * T(2.8) / ω_w
    dtT <= dt_cfl || throw(
        ArgumentError(
            "dt=$dt exceeds the whistler-CFL limit $(round(dt_cfl, sigdigits = 3)) at " *
            "N=$N (dx=$(round(dx, sigdigits = 3))); reduce dt or NB",
        ),
    )

    Te = βeT / 2                                   # β_e = 2 Te (n=1, B=1)
    model = HybridModel(IsothermalElectrons(Te); η = η)
    Np = Int(nppc) * Int(N)
    st = HybridStepper(g, model, CIC(), Np)
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(Int(seed))
    load_uniform!(ps, rng, (zero(T),), (LxT,))
    load_maxwellian!(ps, rng, (zero(T), zero(T), zero(T)), (vpar, vperp, vperp))
    set_density_weight!(ps, one(T), g)
    fill!(st.fields.B[1], one(T))                  # B₀ = x̂
    fill!(st.fields.B[2], zero(T))
    fill!(st.fields.B[3], zero(T))
    init!(st, ps)

    wB0 = T(0.5) * LxT                             # ½ B₀² Lx (background magnetic energy, 1-D)
    wperp() = T(0.5) * (sum(abs2, st.fields.B[2]) + sum(abs2, st.fields.B[3])) * dx
    wmax = wperp()
    wfin = wmax
    ok = 0
    for s = 1:Int(nsteps)
        step!(st, ps, dtT; NB = Int(NB))
        all(isfinite, st.fields.B[2]) && all(isfinite, st.fields.B[3]) || break
        wfin = wperp()
        wmax = max(wmax, wfin)
        ok += 1
    end
    return (;
        wperp_max = wmax,
        wperp_final = wfin,
        wB0 = wB0,
        ratio_max = wmax / wB0,
        anisotropy = vpar^2 - vperp^2,
        unstable_theory = (vpar^2 - vperp^2) > one(T),
        nsamples = ok,
    )
end
