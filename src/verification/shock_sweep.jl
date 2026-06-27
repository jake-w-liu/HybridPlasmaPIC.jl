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
                     nppc=64, nsteps=900, seed=1)
        -> (; n2, Bz2, Vs, X_rh, frozen_ratio, reflected_fraction, M_real, xf)

Set up and run a 1-D perpendicular collisionless shock at upstream Alfvén Mach
number `MA` (upstream drift `U0 = MA·v_A`, with `v_A = 1` in Ω_ci-normalized
units), then return the downstream / shock-front diagnostics as a NamedTuple:

  • `n2`                — downstream compression ρ₂/ρ₁ (upstream n₁=1),
  • `Bz2`               — downstream magnetic field,
  • `Vs`                — shock-front speed from mass conservation (rest frame),
  • `X_rh`              — fluid Rankine–Hugoniot compression at the realized
                          Mach number `M_real = (U0+Vs)/v_A`,
  • `frozen_ratio`      — `(Bz2/B0)/n2` (=1 when the field is frozen to the flow),
  • `reflected_fraction`— weighted fraction of ions classified as reflected,
  • `M_real`            — realized (shock-frame) Mach number,
  • `xf`                — final shock-front position.

The simulation is the verified reflecting-wall `PerpShock` model (SHK-002);
`seed` seeds the particle load so the kinetic noise / seed sensitivity can be
studied. `nsteps`/`nppc`/`N` are kept modest by default for fast research sweeps.
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
)
    N >= 3 || throw(ArgumentError("N must be at least 3"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    T = Float64
    LxT = T(Lx)
    B0 = one(T)
    vA = one(T)
    U0 = T(MA) * vA                      # upstream drift speed (toward the wall)
    vth = T(vthi)

    # field state: inflow SAT strength τ ≈ inflow speed U0 (as in SHK-002)
    sh = PerpShock(N, LxT; Te = T(Te), γe = T(γe), η = T(η), τ = U0, B0 = B0)

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
    ps.weight .= shock_density_weight(one(T), LxT, Np)   # so n₁ = 1
    init_shock!(sh, ps)

    # time-march; sample the front position as a fallback speed diagnostic. The
    # outermost half-amplitude crossing can be corrupted by kinetic Bz ripples in
    # short low-ppc sweeps, so it is no longer the primary speed measurement.
    dt = T(0.02)
    rec_every = max(1, nsteps ÷ 28)        # ~28 samples ⇒ a well-conditioned fit
    pos = T[]
    tt = T[]
    for st = 1:nsteps
        step_shock!(sh, ps, dt; NB = 2)
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

    return (; n2, Bz2, Vs, X_rh, frozen_ratio, reflected_fraction, M_real, xf)
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
