# shock_ramp.jl — §11.3 / Phase-10: initial-ramp generator for the 1-D
# perpendicular hybrid shock plus ramp-width / box-length sensitivity scans.
#
# A perpendicular shock is conveniently *initialized* (rather than self-formed
# from a piston) by laying down a smooth, finite-width tanh transition between an
# upstream state (n1, B1) and a downstream state (n2, B2). With the shock normal
# = x and the convention used here, the UPSTREAM half-space is x > x_ramp and the
# DOWNSTREAM half-space is x < x_ramp (the reflecting wall lives at x = 0, i.e. on
# the downstream side). The field profile is the canonical tanh ramp
#
#     Bz(x) = B1 + (B2 − B1)·½·(1 − tanh((x − x_ramp)/width))
#
# which → B1 as x → +∞ (upstream) and → B2 as x → −∞ (downstream). The ion
# density is initialized to follow the SAME tanh shape,
#
#     n(x)  = n1 + (n2 − n1)·½·(1 − tanh((x − x_ramp)/width)),
#
# so the field is frozen-in (Bz/n constant) at t = 0 when B1/n1 = B2/n2.
#
# Density is realized by WEIGHT MODULATION: particles are placed uniformly over
# [0, Lx] and each particle's weight is set to w_p = n(x_p)·Lx/Np. The CIC
# deposit (deposit_moments!) then reproduces n(x) in the mode-resolved mean to
# within particle-statistics noise (no accept/reject loop is needed and the
# total particle count is fixed, which keeps the SoA layout and id/tag arrays
# untouched). This matches `shock_density_weight` in the uniform limit
# (n(x) ≡ n0 ⇒ w_p = n0·Lx/Np).

"""
    initial_ramp!(sh::PerpShock, ps, x_ramp, width, n1, n2, B1, B2; rng=Random.GLOBAL_RNG)

Initialize a finite-width `tanh` shock ramp on the `PerpShock` field state `sh`
and (re)load the particle set `ps` so the ion density follows the same profile.

The magnetic field is set node-wise to the analytic ramp

    sh.Bz[i] = B1 + (B2 − B1)·½·(1 − tanh((sh.x[i] − x_ramp)/width))

(`→ B1` upstream `x ≫ x_ramp`, `→ B2` downstream `x ≪ x_ramp`). Particles are
placed uniformly on `[0, Lx]` (Lx = sh.x[end]) and weighted by

    w_p = n(x_p)·Lx/Np,   n(x) = n1 + (n2 − n1)·½·(1 − tanh((x − x_ramp)/width)),

so the CIC-deposited density reproduces `n(x)` in the mean (upstream endpoint
`n1`, downstream endpoint `n2`). `ps` must be a 1-D `ParticleSet` whose length is
the desired particle count `Np`; positions, weights, ids and tags are (re)set.
Velocities are left untouched (set them before/after as the experiment needs).
Finally `init_shock!` is called to deposit the moments and compute the carried E.
Returns `sh`.
"""
function initial_ramp!(
    sh::PerpShock{T},
    ps::ParticleSet{1,T},
    x_ramp::Real,
    width::Real,
    n1::Real,
    n2::Real,
    B1::Real,
    B2::Real;
    rng = Random.GLOBAL_RNG,
) where {T}
    xr = T(x_ramp)
    w = T(width)
    w > 0 || throw(ArgumentError("ramp width must be > 0"))
    N1 = T(n1)
    N2 = T(n2)
    Bup = T(B1)
    Bdn = T(B2)
    Lx = sh.x[end]

    # analytic tanh field on the SBP nodes
    @inbounds for i in eachindex(sh.Bz)
        sh.Bz[i] = Bup + (Bdn - Bup) * T(0.5) * (one(T) - tanh((sh.x[i] - xr) / w))
    end

    # uniform placement + tanh weight modulation so the deposited density is n(x)
    Np = nparticles(ps)
    Np > 0 || throw(ArgumentError("particle set must be non-empty"))
    xp = ps.x[1]
    wt = ps.weight
    cell = Lx / Np                                   # = shock_density_weight(1, Lx, Np)
    @inbounds for p = 1:Np
        xx = Lx * rand(rng, T)
        xp[p] = xx
        nx = N1 + (N2 - N1) * T(0.5) * (one(T) - tanh((xx - xr) / w))
        wt[p] = nx * cell
    end
    # refresh provenance arrays to a clean state for the (re)loaded set
    @inbounds for p = 1:Np
        ps.id[p] = UInt64(p)
        ps.tag[p] = zero(UInt32)
    end

    init_shock!(sh, ps)
    return sh
end

"""
    ramp_width_scan(; widths, N=256, Lx=120.0, x_ramp=Lx/2, n1=1.0, n2=3.0,
                      B1=1.0, B2=3.0, vthi=0.35, Te=0.125, γe=5/3, η=0.02,
                      nppc=64, nsteps=60, dt=0.02, NB=2, seed=1)
        -> Vector{NamedTuple}

For each initial ramp `width`, build a `PerpShock`, lay down the tanh ramp with
[`initial_ramp!`](@ref), march it a SHORT time, and report the measured shock
front width (`shock_front`). Returns one
`(; width0, xf, width_measured, n2_meas, Bz_jump)` per input width, where
`width0` is the input ramp width, `width_measured` the measured ramp width after
relaxation, `xf` the front position, `n2_meas` the downstream-slab mean density,
and `Bz_jump` the realized field jump `Bz_down − Bz_up`. Use it to study how the
measured front width tracks the imposed initial width.
"""
function ramp_width_scan(;
    widths,
    N::Integer = 256,
    Lx::Real = 120.0,
    x_ramp::Real = Lx / 2,
    n1::Real = 1.0,
    n2::Real = 3.0,
    B1::Real = 1.0,
    B2::Real = 3.0,
    vthi::Real = 0.35,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    η::Real = 0.02,
    nppc::Integer = 64,
    nsteps::Integer = 60,
    dt::Real = 0.02,
    NB::Integer = 2,
    seed::Integer = 1,
)
    T = Float64
    LxT = T(Lx)
    out = NamedTuple[]
    for wd in widths
        sh = PerpShock(N, LxT; Te = T(Te), γe = T(γe), η = T(η), τ = zero(T), B0 = T(B1))
        Np = Int(nppc) * Int(N)
        ps = ParticleSet{1,T}(Np)
        rng = MersenneTwister(Int(seed))
        initial_ramp!(sh, ps, T(x_ramp), T(wd), T(n1), T(n2), T(B1), T(B2); rng = rng)
        # thermal velocities (no net drift): a quiescent ramp relaxation test
        vth = T(vthi)
        for c = 1:3
            vc = ps.v[c]
            @inbounds for p in eachindex(vc)
                vc[p] = vth * randn(rng)
            end
        end
        init_shock!(sh, ps)
        for _ = 1:Int(nsteps)
            step_shock!(sh, ps, T(dt); NB = Int(NB))
        end
        _require_all_finite("Bz", sh.Bz, "unstable ramp run")
        xf, wmeas = shock_front(sh.Bz, sh.x)
        # downstream slab between wall and front
        dmask = (sh.x .> T(2)) .& (sh.x .< xf - T(2))
        any(dmask) || (dmask = sh.x .< xf)
        n2_meas = _ramp_slab_mean(sh.n, dmask)
        Bz_jump = sh.Bz[1] - sh.Bz[end]
        push!(out, (; width0 = T(wd), xf, width_measured = wmeas, n2_meas, Bz_jump))
    end
    return out
end

"""
    box_length_scan(; Lxs, MA=4.0, kwargs...) -> Vector{NamedTuple}

Run the verified reflecting-wall piston shock ([`run_perp_shock`](@ref)) at
several box lengths `Lxs` (keeping the node spacing comparable by scaling `N`
with `Lx`), and return the measured downstream compression per box. Returns one
`(; Lx, N, n2, Bz2, Vs)` per box. The compression should be INSENSITIVE to `Lx`
once the box is long enough that the shock and its upstream region fit without
the front reaching the inflow boundary during the run.
"""
function box_length_scan(;
    Lxs,
    MA::Real = 4.0,
    N0::Integer = 256,
    Lx0::Real = 120.0,
    nsteps::Integer = 500,
    nppc::Integer = 64,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    seed::Integer = 1,
)
    out = NamedTuple[]
    for Lx in Lxs
        # keep dx ≈ Lx0/N0 by scaling node count with the box length
        Nb = max(64, round(Int, Int(N0) * (Lx / Lx0)))
        r = run_perp_shock(;
            MA = MA,
            N = Nb,
            Lx = Lx,
            Te = Te,
            γe = γe,
            vthi = vthi,
            η = η,
            nppc = nppc,
            nsteps = Int(nsteps),
            seed = Int(seed),
        )
        push!(out, (; Lx = Float64(Lx), N = Nb, n2 = r.n2, Bz2 = r.Bz2, Vs = r.Vs))
    end
    return out
end

# arithmetic mean of `v` over masked nodes (NaN if none) — local helper
function _ramp_slab_mean(v::AbstractVector{T}, mask::AbstractVector{Bool}) where {T}
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
