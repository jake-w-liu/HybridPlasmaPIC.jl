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

The relativistic momentum is `p = γ v`; after its drag substep, recover `v` via
`γ = √(1 + |p|²/c²)`. `K` is the radiation-reaction coefficient
(∝ the classical electron radius). A gyrating particle in a magnetic field cools
synchrotron-like, `d(1/γ)/dt = K v⁴ B²/c² ≈ K c² B²` ultrarelativistically (so
`1/γ ≈ 1/γ₀ + K c² B² t`; this reduces to the familiar `K B²` rate only when `c = 1`).

The exponential substep `p ← p·exp(−K γ χ² dt)` holds the initial drag rate fixed.
It approximates the nonlinear drag equation; it is exact only for that frozen
rate. Extreme finite scales can still lose numeric range in the invariant or
momentum arithmetic. `K ≥ 0`, `dt ≥ 0`, `c > 0`; these
parameters and `E, B` must remain finite in the particle precision. Inputs are
validated before a `K=0` or `dt=0` no-op. Returns `ps`.
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
    KT = _require_finite_nonnegative_real("radiation-reaction coefficient K", K, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "apply_radiation_reaction!")
    cT = _require_finite_positive_real("c", c, T)
    Ex, Ey, Ez = _validated_uniform_field("E", E, T)
    Bx, By, Bz = _validated_uniform_field("B", B, T)
    (iszero(KT) || iszero(dtT)) && return ps
    _require_subluminal_velocities(ps, cT, "apply_radiation_reaction!")
    c2 = cT * cT
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(vx)
        ux = vx[p]
        uy = vy[p]
        uz = vz[p]
        v2 = ux * ux + uy * uy + uz * uz
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
