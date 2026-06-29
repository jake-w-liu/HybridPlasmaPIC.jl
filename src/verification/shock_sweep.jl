# shock_sweep.jl — Phase 11 (1-D parts): collisionless perpendicular-shock
# research sweep over the upstream Alfvén Mach number M_A.
#
# Drives the verified reflecting-wall perpendicular hybrid shock (`PerpShock` /
# `step_shock!`, SHK-002): upstream ions drift at U0 = M_A·v_A (v_A=1, Ω_ci-units)
# into the piston wall at x=0; a shock forms and propagates back toward the
# inflow boundary at x=Lx. Each run measures
#   • downstream compression n2 and field Bz2 (averaged over the downstream slab),
#   • shock speed Vs from mass conservation U0/(n2−1),
#   • the frozen-in ratio (Bz2/B0)/n2 — should be ≈1 for flux freezing,
#   • the fluid Rankine–Hugoniot compression X_rh at the realized Mach
#     (U0+Vs)/v_A via `rankine_hugoniot` (an independent analytic oracle),
#   • the reflected-ion fraction via `classify_reflected`.
#
# The kinetic shock compresses LESS than the γ=5/3 fluid (hotter downstream with
# reflected ions); the EOS-independent checks — flux freezing and mass
# conservation Vs ≈ U0/(n2−1) — are the tight, physics-grade diagnostics.

"""
    run_perp_shock(; MA, N=512, Lx=120, Te=0.125, γe=5/3, vthi=0.35, η=0.02,
                     nppc=64, nsteps=900, seed=1, inject=false)
        -> (; n2, Bz2, Vs, X_rh, frozen_ratio, reflected_fraction,
              reflected_flux, M_real, xf)

Set up and run a 1-D perpendicular collisionless shock at upstream Alfvén Mach
number `MA` (upstream drift `U0 = MA·v_A`, with `v_A = 1` in Ω_ci-normalized
units), then return the downstream / shock-front diagnostics as a NamedTuple:

  • `n2`                — downstream compression ρ₂/ρ₁ (upstream n₁=1),
  • `Bz2`               — downstream magnetic field,
  • `Vs`                — shock-front speed from mass conservation (rest frame),
  • `X_rh`              — fluid Rankine–Hugoniot compression at the realized
                          Mach number `M_real = (U0+Vs)/v_A`,
  • `frozen_ratio`      — `(Bz2/B0)/n2` (=1 when the field is frozen to the flow),
  • `reflected_fraction`— whole-box weighted fraction of ions classified as reflected,
  • `reflected_flux`    — literature-style α = reflected flux at the front /
                          upstream flux `n₁·V₁` ([`reflected_flux_fraction`](@ref)),
  • `M_real`            — realized (shock-frame) Mach number,
  • `xf`                — final shock-front position.

The simulation is the verified reflecting-wall `PerpShock` model (SHK-002);
`seed` seeds the particle load so the kinetic noise / seed sensitivity can be
studied. `nsteps`/`nppc`/`N` are kept modest by default for fast research sweeps.
Set `inject=true` for sustained upstream injection (maintains `n₁=1` at the inflow
so the shock can run to a quasi-stationary state instead of draining the finite
reservoir) — needed to compare against sustained shocks (e.g. Leroy 1982).
"""
function run_perp_shock(;
    MA::Real,
    N::Integer = 512,
    Lx::Real = 120.0,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    nppc::Integer = 64,
    nsteps::Integer = 900,
    seed::Integer = 1,
    inject::Bool = false,
)
    N >= 3 || throw(ArgumentError("N must be at least 3"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    T = Float64
    MAT = _require_valid_positive_shock_ma(MA, T)
    LxT = _require_finite_positive_real("Lx", Lx, T)
    B0 = one(T)
    vA = one(T)
    U0 = MAT * vA                        # upstream drift speed (toward the wall)
    vth = _require_finite_nonnegative_real("vthi", vthi, T)

    # field state: inflow SAT strength τ ≈ inflow speed U0 (as in SHK-002)
    sh = PerpShock(N, LxT; Te, γe, η, τ = U0, B0 = B0)

    # load drifting ions: x uniform on [0,Lx], v drifting at −U0 toward the wall
    Np = Int(nppc) * Int(N)
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(Int(seed))
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    @inbounds for p = 1:Np
        xp[p] = LxT * rand(rng)
        vx[p] = -U0 + vth * randn(rng)
        vy[p] = vth * randn(rng)
        vz[p] = vth * randn(rng)
    end
    wp = shock_density_weight(one(T), LxT, Np)
    ps.weight .= wp                                      # so n₁ = 1
    init_shock!(sh, ps)

    # opt-in sustained upstream injection (maintains n₁=1 at the inflow so the
    # shock can reach a quasi-stationary state instead of depleting the reservoir).
    injector =
        inject ?
        ShockInjector(
            MersenneTwister(Int(seed) + 1);
            n0 = one(T),
            drift = U0,
            vthi = vth,
            σt = vth,
            weight = wp,
            first_id = Np + 1,
        ) : nothing

    # time-march; sample the front position as a fallback speed diagnostic. The
    # outermost half-amplitude crossing can be corrupted by kinetic Bz ripples in
    # short low-ppc sweeps, so it is no longer the primary speed measurement.
    dt = T(0.02)
    rec_every = max(1, nsteps ÷ 28)        # ~28 samples ⇒ a well-conditioned fit
    pos = T[]
    tt = T[]
    for st = 1:nsteps
        step_shock!(sh, ps, dt; NB = 2, injector = injector)
        if st % rec_every == 0
            push!(pos, _front_crossing(sh.Bz, sh.x, B0))
            push!(tt, st * dt)
        end
    end
    _require_all_finite("Bz", sh.Bz, "unstable run")

    # reported front position uses the steepest-gradient locator (shock_front)
    xf, _ = shock_front(sh.Bz, sh.x)

    # downstream slab: between the wall and the front (avoid both edges)
    dmask = (sh.x .> T(5)) .& (sh.x .< xf - T(5))
    if !any(dmask)                         # front too close to the wall: fall back
        dmask = (sh.x .> T(2)) .& (sh.x .< max(xf - T(2), T(3)))
    end
    n2 = _slab_mean(sh.n, dmask)
    Bz2 = _slab_mean(sh.Bz, dmask)

    # Shock speed from mass conservation, Vs = U0/(n2−1). This is the robust
    # EOS-independent measure used by the 3-D shock diagnostics too. The sampled
    # front track is retained only as a fallback for non-compressive/invalid runs:
    # kinetic Bz ripples can create distant half-level crossings and corrupt a
    # slope fit in short, low-ppc sweeps.
    Vs_mass = isfinite(n2) && n2 > one(T) ? U0 / (n2 - one(T)) : T(NaN)
    Vs = isfinite(Vs_mass) ? Vs_mass : _front_speed(pos, tt, LxT)

    # realized shock-frame Mach number and fluid RH compression at that Mach
    M_real = (U0 + Vs) / vA
    p1 = vth^2 + T(Te)                     # upstream pressure (ion + electron)
    rh = rankine_hugoniot(MHDState(one(T), U0 + Vs, zero(T), p1, zero(T), B0), T(γe))
    X_rh = rh.X

    frozen_ratio = (Bz2 / B0) / n2

    # reflected-ion fraction (weighted): shock front moving in +x at Vs, plasma
    # upstream at +x. classify_reflected uses the rest-frame front position and
    # the front speed in the same (rest) frame.
    refl = classify_reflected(ps, xf, Vs)
    wtot = zero(T)
    wrefl = zero(T)
    w = ps.weight
    @inbounds for p in eachindex(w)
        wtot += w[p]
        refl[p] && (wrefl += w[p])
    end
    reflected_fraction = wtot > 0 ? wrefl / wtot : zero(T)

    # literature-style reflected fraction α (flux at the front / upstream flux
    # n₁·V₁), the apples-to-apples quantity for Leroy 1982 / Hellinger 2002.
    reflected_flux =
        isfinite(xf) && isfinite(Vs) && isfinite(M_real) ?
        reflected_flux_fraction(ps, xf, Vs, U0 + Vs) : T(NaN)

    return (; n2, Bz2, Vs, X_rh, frozen_ratio, reflected_fraction, reflected_flux, M_real, xf)
end

# mean / population std / median of a Float64 vector (NaN-safe, no Statistics dep)
function _avg_std(v::AbstractVector{T}) where {T}
    isempty(v) && return (T(NaN), T(NaN))
    m = sum(v) / length(v)
    s = sqrt(sum(x -> (x - m)^2, v) / length(v))
    return m, s
end
function _median(v::AbstractVector{T}) where {T}
    isempty(v) && return T(NaN)
    s = sort(v)
    n = length(s)
    return isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1]) / 2
end

"""
    run_perp_shock_rh(; MA, β=1.0, N=512, Lx=200, γe=5/3, η=0.02, nppc=160,
                        nsteps=1500, seed=1, t_avg_start=8.0)
        -> (; compression, compression_std, reflected_flux, reflected_flux_std,
              overshoot, overshoot_std, X_rh, M_real, nsamples)

Sustained perpendicular shock initialized from the **two-state Rankine–Hugoniot
profile** (Leroy et al. 1982 setup): the downstream half is pre-loaded at the
fluid-RH compression and temperature, the upstream at `(n=1, B=B0, u=−U0)`, joined
by a thin `tanh` ramp, and upstream plasma is injected at the inflow so the shock
runs to a quasi-stationary state. Diagnostics (downstream median compression, the
front overshoot `B_max/B₂`, and the flux-based reflected fraction α) are
time-averaged once `t ≥ t_avg_start`. Use this — not the piston `run_perp_shock` —
to compare against sustained published perpendicular shocks.

`MA` is the shock-frame Alfvén Mach number; `β = β_e = β_i` sets the upstream
temperatures (`Te = β/2`, `vthi = √(β/2)`).
"""
function run_perp_shock_rh(;
    MA::Real,
    β::Real = 1.0,
    N::Integer = 512,
    Lx::Real = 200.0,
    γe::Real = 5 / 3,
    η::Real = 0.02,
    nppc::Integer = 160,
    nsteps::Integer = 1500,
    seed::Integer = 1,
    t_avg_start::Real = 8.0,
)
    N >= 8 || throw(ArgumentError("N must be ≥ 8"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    T = Float64
    MAT = _require_valid_positive_shock_ma(MA, T)
    LxT = _require_finite_positive_real("Lx", Lx, T)
    βT = _require_finite_nonnegative_real("β", β, T)
    B0 = one(T)
    Te = βT / 2                              # β_e = 2 Te (n=1, B=1)
    vthi = sqrt(βT / 2)                      # β_i = 2 vthi²
    p1 = vthi^2 + Te                         # upstream total pressure

    # fluid Rankine–Hugoniot downstream at shock-frame Mach MA (perpendicular: Bn=0)
    rh = rankine_hugoniot(MHDState(one(T), MAT, zero(T), p1, zero(T), B0), T(γe))
    n2 = rh.X
    n2 > one(T) || throw(ArgumentError("no compressive RH shock at MA=$MA, β=$β"))
    Bz2 = n2 * B0
    pe2 = Te * n2^T(γe)
    Ti2 = max((rh.down.p - pe2) / n2, zero(T))
    vth2 = sqrt(Ti2)                          # downstream ion thermal speed
    U0 = MAT * (n2 - one(T)) / n2             # so M_real = U0·n2/(n2−1) = MA
    Vs = U0 / (n2 - one(T))

    sh = PerpShock(N, LxT; Te, γe, η, τ = U0, B0)
    Np = Int(nppc) * Int(N)
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(Int(seed))
    xr = LxT / 3                              # initial shock position (downstream = x<xr)
    initial_ramp!(sh, ps, xr, T(2), one(T), n2, B0, Bz2; rng = rng)
    @inbounds for p = 1:nparticles(ps)
        if ps.x[1][p] < xr                    # downstream: at rest (wall frame), hot
            ps.v[1][p] = vth2 * randn(rng)
            ps.v[2][p] = vth2 * randn(rng)
            ps.v[3][p] = vth2 * randn(rng)
        else                                  # upstream: drift −U0, cool
            ps.v[1][p] = -U0 + vthi * randn(rng)
            ps.v[2][p] = vthi * randn(rng)
            ps.v[3][p] = vthi * randn(rng)
        end
    end
    init_shock!(sh, ps)

    wp = shock_density_weight(one(T), LxT, Np)
    injector = ShockInjector(
        MersenneTwister(Int(seed) + 1);
        n0 = one(T),
        drift = U0,
        vthi = vthi,
        σt = vthi,
        weight = wp,
        first_id = Np + 1,
    )

    dt = T(0.02)
    comp = T[]
    over = T[]
    alpha = T[]
    for st = 1:nsteps
        step_shock!(sh, ps, dt; NB = 2, injector = injector)
        if st * dt >= T(t_avg_start) && st % 20 == 0
            xf, _ = shock_front(sh.Bz, sh.x)
            isfinite(xf) || continue
            T(15) < xf < LxT - T(20) || continue       # shock well inside the box
            dmask = (sh.x .> T(8)) .& (sh.x .< xf - T(8))
            any(dmask) || continue
            nd = sh.n[dmask]
            Bzd = sh.Bz[dmask]
            ndm = _median(nd)
            ndm > one(T) || continue
            push!(comp, ndm)
            fmask = (sh.x .> xf - T(8)) .& (sh.x .< xf + T(2))
            push!(over, maximum(sh.Bz[fmask]) / _median(Bzd))
            Vsr = U0 / (ndm - one(T))
            push!(alpha, reflected_flux_fraction(ps, xf, Vsr, U0 + Vsr))
        end
    end
    cm, cs = _avg_std(comp)
    om, os = _avg_std(over)
    am, as = _avg_std(alpha)
    return (;
        compression = cm,
        compression_std = cs,
        reflected_flux = am,
        reflected_flux_std = as,
        overshoot = om,
        overshoot_std = os,
        X_rh = n2,
        M_real = U0 + Vs,
        nsamples = length(comp),
    )
end

"""
    run_perp_shock_leroy(; MA, β=1.0, N=512, Lx=200.0, γe=5/3, η=0.02,
                           nppc=160, nsteps=1500, seed=1, t_avg_start=8.0, window=8.0)

§11.3 Rankine–Hugoniot two-state shock in the wall-less, shock-REST frame
(Leroy et al. 1982 setup): upstream plasma flows IN at x=Lx at the shock-frame
speed `V1 = MA`, downstream flows OUT at x=0, the shock is held stationary by the
two-ended flux balance, and the downstream boundary is a thermal RESERVOIR (not a
specular wall) so self-consistent ion reflection / the energetic foot can develop.

Unlike [`run_perp_shock_rh`] (reflecting-wall, §11.2), this is the configuration
that produces Leroy's reflected fraction α: returns the same fields, with
`reflected_flux` the rest-frame α (flux of back-streaming ions in `[xf, xf+window]`
over the upstream flux n₁·V1) and `M_real = MA` (the inflow Mach, exact by
construction).
"""
function run_perp_shock_leroy(;
    MA::Real,
    β::Real = 1.0,
    N::Integer = 512,
    Lx::Real = 200.0,
    γe::Real = 5 / 3,
    η::Real = 0.02,
    nppc::Integer = 160,
    nsteps::Integer = 1500,
    seed::Integer = 1,
    t_avg_start::Real = 8.0,
    window::Real = 8.0,
)
    N >= 8 || throw(ArgumentError("N must be ≥ 8"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    T = Float64
    MAT = _require_valid_positive_shock_ma(MA, T)
    LxT = _require_finite_positive_real("Lx", Lx, T)
    βT = _require_finite_nonnegative_real("β", β, T)
    winT = _require_finite_positive_real("window", window, T)
    B0 = one(T)
    Te = βT / 2
    vthi = sqrt(βT / 2)
    p1 = vthi^2 + Te

    rh = rankine_hugoniot(MHDState(one(T), MAT, zero(T), p1, zero(T), B0), T(γe))
    n2 = rh.X
    n2 > one(T) || throw(ArgumentError("no compressive RH shock at MA=$MA, β=$β"))
    Bz2 = n2 * B0
    Ti2 = max((rh.down.p - Te * n2^T(γe)) / n2, zero(T))
    vth2 = sqrt(Ti2)
    V1 = MAT                                  # shock-frame upstream inflow speed
    V2 = V1 / n2                              # downstream outflow speed (n₁V1 = n₂V2)
    p_up = V2 / flux_per_density(V2, vth2)    # exiting-ion fraction recycled upstream

    sh = PerpShock(N, LxT; Te, γe, η, τ = V1, B0)
    Np = Int(nppc) * Int(N)
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(Int(seed))
    xr = LxT / 2                              # shock sits mid-box (room both sides)
    initial_ramp!(sh, ps, xr, T(2), one(T), n2, B0, Bz2; rng = rng)
    @inbounds for p = 1:nparticles(ps)
        if ps.x[1][p] < xr                    # downstream: flows out at −V2, hot
            ps.v[1][p] = -V2 + vth2 * randn(rng)
            ps.v[2][p] = vth2 * randn(rng)
            ps.v[3][p] = vth2 * randn(rng)
        else                                  # upstream: flows in at −V1, cool
            ps.v[1][p] = -V1 + vthi * randn(rng)
            ps.v[2][p] = vthi * randn(rng)
            ps.v[3][p] = vthi * randn(rng)
        end
    end
    init_shock!(sh, ps)

    bc = LeroyBoundary(
        MersenneTwister(Int(seed) + 1);
        V1 = V1,
        vthi = vthi,
        V2 = V2,
        vth2 = vth2,
        B_down = Bz2,
        p_up = p_up,
    )

    dt = T(0.02)
    comp = T[]
    over = T[]
    alpha = T[]
    for st = 1:nsteps
        step_leroy_shock!(sh, ps, dt; NB = 2, bc = bc)
        if st * dt >= T(t_avg_start) && st % 20 == 0
            xf, _ = shock_front(sh.Bz, sh.x)
            isfinite(xf) || continue
            T(12) < xf < LxT - T(12) || continue
            dmask = (sh.x .> T(8)) .& (sh.x .< xf - T(8))      # downstream = x < xf
            any(dmask) || continue
            ndm = _median(sh.n[dmask])
            ndm > one(T) || continue
            push!(comp, ndm)
            fmask = (sh.x .> xf - T(8)) .& (sh.x .< xf + T(2))
            push!(over, maximum(sh.Bz[fmask]) / _median(sh.Bz[dmask]))
            # rest frame: Vs=0, reflected ions back-stream (vx>0) into [xf, xf+window]
            push!(alpha, reflected_flux_fraction(ps, xf, zero(T), V1; window = winT, n1 = one(T)))
        end
    end
    cm, cs = _avg_std(comp)
    om, os = _avg_std(over)
    am, as = _avg_std(alpha)
    return (;
        compression = cm,
        compression_std = cs,
        reflected_flux = am,
        reflected_flux_std = as,
        overshoot = om,
        overshoot_std = os,
        X_rh = n2,
        M_real = MAT,
        nsamples = length(comp),
    )
end

# weighted-free arithmetic mean of `v` over the masked nodes (NaN if empty)
function _slab_mean(v::AbstractVector{T}, mask::AbstractVector{Bool}) where {T}
    s = zero(T)
    c = 0
    @inbounds for i in eachindex(v)
        if mask[i]
            s += v[i]
            c += 1
        end
    end
    return c > 0 ? s / c : T(NaN)
end

# outermost (largest-x) crossing of the half-amplitude level between the
# downstream field Bz[1] (wall side) and the upstream B0 (inflow side). Linearly
# interpolated; returns the wall node if no crossing exists yet (pre-formation).
function _front_crossing(Bz::AbstractVector{T}, x::AbstractVector{T}, B0::T) where {T}
    n = length(Bz)
    lvl = (Bz[1] + B0) / 2
    @inbounds for i = n-1:-1:1
        d0 = Bz[i] - lvl
        d1 = Bz[i+1] - lvl
        if d0 * d1 <= 0 && Bz[i] != Bz[i+1]
            f = (lvl - Bz[i]) / (Bz[i+1] - Bz[i])
            return x[i] + f * (x[i+1] - x[i])
        end
    end
    return x[1]
end

# least-squares slope dpos/dt over the back half of the sampled front track
# (the early samples include the wall start-up transient). Pre-formation samples
# pinned at/near the inflow boundary `Lx` are dropped first. Falls back to the
# full clean record / a finite difference for very short tracks.
function _front_speed(pos::AbstractVector{T}, tt::AbstractVector{T}, Lx::T) where {T}
    keep = [pos[i] > T(1) && pos[i] < Lx - T(1) for i in eachindex(pos)]
    p = pos[keep]
    t = tt[keep]
    m = length(p)
    m >= 2 || return T(NaN)
    i0 = max(1, m ÷ 2)                     # back half
    (m - i0 + 1) >= 2 || (i0 = m - 1)      # need ≥2 points
    xs = @view t[i0:m]
    ys = @view p[i0:m]
    n = length(xs)
    sx = zero(T)
    sy = zero(T)
    sxx = zero(T)
    sxy = zero(T)
    @inbounds for k = 1:n
        sx += xs[k]
        sy += ys[k]
        sxx += xs[k]^2
        sxy += xs[k] * ys[k]
    end
    denom = n * sxx - sx^2
    return denom != 0 ? (n * sxy - sx * sy) / denom : (ys[end] - ys[1]) / (xs[end] - xs[1])
end

"""
    perp_shock_sweep(MAs; kwargs...) -> Vector{NamedTuple}

Run [`run_perp_shock`](@ref) for each Alfvén Mach number in `MAs`, forwarding the
keyword arguments, and return the per-`MA` diagnostic NamedTuples (each augmented
with the input `MA`).
"""
function perp_shock_sweep(MAs; kwargs...)
    out = NamedTuple[]
    for MA in MAs
        r = run_perp_shock(; MA = MA, kwargs...)
        push!(out, (; MA = MA, r...))
    end
    return out
end
