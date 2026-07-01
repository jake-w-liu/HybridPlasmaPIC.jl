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

"""
    ion_cyclotron_growth(; vth_par, vth_perp, kwargs...)

Drive the **electromagnetic ion-cyclotron (EMIC) anisotropy instability** — the
`T_⊥ > T_∥` counterpart of the firehose (same parallel bi-Maxwellian run and
transverse-δB diagnostic, opposite anisotropy). Reuses [`firehose_growth`] for the
run and adds the EMIC threshold: with anisotropy `A = T_⊥/T_∥ − 1` and `β_∥ = 2T_∥`,
the plasma is EMIC-unstable when `A > 0.43 / β_∥^0.43` (Gary 1993 marginal-stability
fit). Returns `(; wperp_max, wperp_final, wB0, ratio_max, T_anisotropy, beta_par,
unstable_theory, nsamples)`.
"""
function ion_cyclotron_growth(; vth_par::Real, vth_perp::Real, kwargs...)
    r = firehose_growth(; vth_par = vth_par, vth_perp = vth_perp, kwargs...)
    T = Float64
    β_par = 2 * T(vth_par)^2                              # β_∥ = 2 T_∥ (n=1, B=1)
    A = (T(vth_perp) / T(vth_par))^2 - 1                  # T_⊥/T_∥ − 1
    A_c = 0.43 / β_par^T(0.43)                            # Gary 1993 EMIC threshold
    return (;
        r.wperp_max,
        r.wperp_final,
        r.wB0,
        r.ratio_max,
        T_anisotropy = A + 1,
        beta_par = β_par,
        unstable_theory = A > A_c,
        r.nsamples,
    )
end

"""
    weibel_growth(; u0, vth=0.1, N=(8,96), L=(4π,12π), nppc=50, c=3.0,
                    dt=0.05, nsteps=600, seed=1)

Drive the **Weibel / current-filamentation instability** in the full electromagnetic
PIC model (kinetic electrons + immobile neutralizing ion background `+n₀`). Two
counter-streaming cold electron beams `±u₀ x̂` carry the free energy — an effective
velocity-space anisotropy `⟨vₓ²⟩ = u₀² + vth² > vth² = ⟨v_y²⟩` — and the unstable
wavevector is `k ∥ ŷ` (⊥ the streaming), so the out-of-plane field `B_z(y)` grows from
particle shot noise. The full-PIC electron dynamics are essential: in the quasineutral
massless-electron hybrid model `B = 0` is an exact fixed point (E is the curl-free `∇pₑ`
term, so `−∇×E = 0`) and only the *bulk* ion moments couple to the field, so the beams'
`u_y` cancels and the ion-Weibel decays — verified. Here `B_z` grows exponentially, then
saturates as the beams isotropize.

    unstable ⇔ A = (u₀/vth)² > 1   (bimodal counter-streaming: streaming anisotropy
                                    exceeds the thermal spread; Weibel 1959, Fried 1959)

Returns `(; wBz_max, anisotropy, unstable_theory, nsamples)` where `wBz_max` is the peak
`B_z` magnetic energy `½∫B_z² dV` and `anisotropy = (u₀/vth)²`.
"""
function weibel_growth(;
    u0::Real,
    vth::Real = 0.1,
    N::NTuple{2,Integer} = (8, 96),
    L::NTuple{2,Real} = (4π, 12π),
    nppc::Integer = 50,
    c::Real = 3.0,
    dt::Real = 0.05,
    nsteps::Integer = 600,
    seed::Integer = 1,
)
    T = Float64
    all(n -> n >= 4, N) || throw(ArgumentError("both N components must be ≥ 4"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    u0T = _require_finite_nonnegative_real("u0", u0, T)      # 0 allowed (stable reference)
    vthT = _require_finite_positive_real("vth", vth, T)
    cT = _require_finite_positive_real("c", c, T)
    dtT = _require_finite_positive_real("dt", dt, T)
    LT = (
        _require_finite_positive_real("L[1]", L[1], T),
        _require_finite_positive_real("L[2]", L[2], T),
    )

    g = FourierGrid((Int(N[1]), Int(N[2])), LT)
    # spectral-leapfrog EM Courant: transverse wave ω=c·k, k_max≈π/dx ⇒ c·dt ≲ 0.64 dx
    cT * dtT <= T(0.6) * minimum(g.dx) || throw(
        ArgumentError(
            "c·dt=$(round(cT * dtT, sigdigits = 3)) exceeds the EM Courant limit " *
            "0.6·min(dx)=$(round(T(0.6) * minimum(g.dx), sigdigits = 3)); reduce dt or c",
        ),
    )

    Np = Int(nppc) * prod(Int.(N))
    es = EMPIC(g, Np; n0 = one(T), c = cT)
    e = ParticleSet{2,T}(Np; q = -one(T), m = one(T))
    rng = MersenneTwister(Int(seed))
    load_uniform!(e, rng, (zero(T), zero(T)), LT)
    load_maxwellian!(e, rng, (zero(T), zero(T), zero(T)), (vthT, vthT, vthT))
    half = Np ÷ 2                                           # two counter-streaming halves (±u₀ x̂)
    @inbounds for p = 1:Np
        e.v[1][p] += p <= half ? u0T : -u0T
    end
    set_density_weight!(e, one(T), g)
    init_empic!(es, e)

    dV = prod(g.dx)
    Bz = es.B[3]
    wBz() = T(0.5) * sum(abs2, Bz) * dV
    wmax = wBz()
    ok = 0
    for s = 1:Int(nsteps)
        step_empic!(es, e, dtT)
        all(isfinite, Bz) || break
        wmax = max(wmax, wBz())
        ok += 1
    end
    A = (u0T / vthT)^2
    return (; wBz_max = wmax, anisotropy = A, unstable_theory = A > one(T), nsamples = ok)
end
