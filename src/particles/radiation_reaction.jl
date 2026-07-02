# radiation_reaction.jl — Landau-Lifshitz radiation-reaction drag (Phase-4 forcing).
#
# A relativistic charged particle radiating in strong fields feels a back-reaction. This
# applies the dominant (γ²) term of the reduced Landau-Lifshitz force as a post-push
# velocity kick (uniform fields), analogous to how the collision operators act on velocity.

"""
    apply_radiation_reaction!(ps::ParticleSet{D,T}, E, B, dt; K, c=1.0) -> ps

One substep of the **Landau-Lifshitz radiation-reaction** drag (dominant `γ²` term) in
spatially uniform fields `E, B` (3-tuples, per unit `q/m`). For each particle of velocity
`v`, `γ = 1/√(1 − |v|²/c²)`:

    w = E + v×B,   χ² = |w|² − (v·E/c)²,   F_rad = −K γ² χ² v

The relativistic momentum `p = γ v` is advanced `p ← p + F_rad·dt` and `v` recovered via
`γ = √(1 + |p|²/c²)`, so `|v| < c` is preserved. `K` is the radiation-reaction coefficient
(∝ the classical electron radius). A gyrating ultrarelativistic particle in a magnetic
field cools synchrotron-like, `d(1/γ)/dt ≈ K B²` (so `1/γ ≈ 1/γ₀ + K B² t`).

The momentum decay is integrated exactly as `p ← p·exp(−K γ χ² dt)` (an exponential
integrator), so it is **unconditionally stable**: `|v|` always shrinks toward zero, never
reverses or spuriously heats, for any `K·dt`. `K ≥ 0`, `dt ≥ 0`, `c > 0`; `E, B` finite;
`K=0` or `dt=0` is a no-op. Returns `ps`.
"""
function apply_radiation_reaction!(
    ps::ParticleSet{D,T},
    E::NTuple{3,<:Real},
    B::NTuple{3,<:Real},
    dt::Real;
    K::Real,
    c::Real = 1.0,
) where {D,T}
    K >= 0 || throw(ArgumentError("radiation-reaction coefficient K must be ≥ 0"))
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    cT = _require_finite_positive_real("c", c, T)
    (all(isfinite, E) && all(isfinite, B)) ||
        throw(ArgumentError("apply_radiation_reaction!: E and B must be finite"))
    (K == 0 || dt == 0) && return ps
    Ex, Ey, Ez = T(E[1]), T(E[2]), T(E[3])
    Bx, By, Bz = T(B[1]), T(B[2]), T(B[3])
    KT = T(K)
    dtT = T(dt)
    c2 = cT * cT
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(vx)
        ux = vx[p]
        uy = vy[p]
        uz = vz[p]
        v2 = ux * ux + uy * uy + uz * uz
        v2 < c2 || continue                                # skip (numerically) superluminal
        γ = one(T) / sqrt(one(T) - v2 / c2)
        wx = Ex + (uy * Bz - uz * By)                      # w = E + v×B
        wy = Ey + (uz * Bx - ux * Bz)
        wz = Ez + (ux * By - uy * Bx)
        vE = ux * Ex + uy * Ey + uz * Ez
        χ2 = (wx * wx + wy * wy + wz * wz) - (vE * vE) / c2 # radiation invariant
        χ2 > 0 || continue
        # dp/dt = F_rad = −K γ² χ² v = −(K γ χ²) p, so p decays by exactly exp(−K γ χ² dt)
        # over the substep (an exponential integrator): |p| always shrinks, never reverses or
        # heats — unconditionally stable, unlike the 1st-order γ − K γ² χ² dt which overshoots.
        s = γ * exp(-KT * γ * χ2 * dtT)
        px = s * ux
        py = s * uy
        pz = s * uz
        γn = sqrt(one(T) + (px * px + py * py + pz * pz) / c2)
        vx[p] = px / γn
        vy[p] = py / γn
        vz[p] = pz / γn
    end
    return ps
end
