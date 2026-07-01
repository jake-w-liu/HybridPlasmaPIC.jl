# collisions.jl — BGK collision operator (§3.3 / §21.5). Relaxes the particle
# velocity distribution toward the LOCAL drifting Maxwellian at collision rate ν,
# while conserving total momentum Σ w v and total kinetic energy Σ w |v|² of the
# whole set EXACTLY (to roundoff).
#
# Standard conservative BGK scatter: a fraction (1 − exp(−ν dt)) of the particles
# are resampled from a Maxwellian built from the set's own mean velocity u and
# temperature T; the scattered subset is then momentum- and energy-corrected
# (shift + rescale) so the set totals are unchanged. Because only the scattered
# subset is touched and that subset's own momentum/energy are restored to their
# pre-scatter values, the global totals are conserved independently of the
# untouched particles.

"""
    collide_bgk!(ps::ParticleSet{D,T}, ν, dt; rng=Random.default_rng()) -> ps

Apply one BGK collision substep of duration `dt` at collision frequency `ν`.

The velocity distribution relaxes toward the local drifting Maxwellian: a random
fraction `1 − exp(−ν·dt)` of the particles is resampled from an isotropic
Maxwellian with the set's weighted mean velocity `u` and scalar temperature
`T = Σ w |v−u|² / (3 Σ w)` (so `v_c = u_c + √T·𝒩(0,1)`). The scattered subset is
then corrected (a uniform velocity shift to restore its momentum, followed by an
isotropic rescaling of its velocity fluctuations to restore its energy) so that
the whole-set totals `Σ w v` and `Σ w |v|²` are conserved exactly to roundoff.

`ν ≥ 0` and `dt ≥ 0` are required. `ν·dt = 0` or fewer than two particles is a
no-op. Returns `ps`.
"""
function collide_bgk!(
    ps::ParticleSet{D,T},
    ν::Real,
    dt::Real;
    rng = Random.default_rng(),
) where {D,T}
    ν >= 0 || throw(ArgumentError("collision frequency ν must be ≥ 0"))
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    N = nparticles(ps)
    # Need at least two particles to have any fluctuation to relax/correct.
    (N < 2 || ν == 0 || dt == 0) && return ps

    vx, vy, vz = ps.v
    w = ps.weight

    # --- whole-set weighted mean velocity u and scalar temperature T ----------
    Wtot = zero(T)
    sx = zero(T)
    sy = zero(T)
    sz = zero(T)
    @inbounds for p = 1:N
        wp = w[p]
        Wtot += wp
        sx += wp * vx[p]
        sy += wp * vy[p]
        sz += wp * vz[p]
    end
    Wtot > 0 || return ps                       # all-zero weights: nothing to do
    ux = sx / Wtot
    uy = sy / Wtot
    uz = sz / Wtot

    s2 = zero(T)
    @inbounds for p = 1:N
        wp = w[p]
        dx = vx[p] - ux
        dy = vy[p] - uy
        dz = vz[p] - uz
        s2 += wp * (dx * dx + dy * dy + dz * dz)
    end
    # scalar temperature (per-component variance) → isotropic thermal speed
    Tscalar = s2 / (3 * Wtot)
    vth = sqrt(max(Tscalar, zero(T)))

    # --- select the scattered subset (Bernoulli with p = 1 − exp(−ν dt)) ------
    Pcoll = -expm1(-T(ν) * T(dt))                # = 1 − exp(−ν dt), accurate near 0
    sel = falses(N)
    nsel = 0
    @inbounds for p = 1:N
        if rand(rng, T) < Pcoll
            sel[p] = true
            nsel += 1
        end
    end
    # A subset of <2 particles carries no fluctuation we can rescale, so leave it
    # untouched — its momentum/energy are already "conserved" trivially.
    nsel < 2 && return ps

    # --- record the subset's pre-scatter weighted momentum and energy ---------
    Wsub = zero(T)
    Px = zero(T)
    Py = zero(T)
    Pz = zero(T)
    Eold = zero(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p]
        Wsub += wp
        Px += wp * vx[p]
        Py += wp * vy[p]
        Pz += wp * vz[p]
        Eold += wp * (vx[p]^2 + vy[p]^2 + vz[p]^2)
    end
    Wsub > 0 || return ps

    # --- resample the subset from the local drifting Maxwellian ---------------
    @inbounds for p = 1:N
        sel[p] || continue
        vx[p] = ux + vth * randn(rng, T)
        vy[p] = uy + vth * randn(rng, T)
        vz[p] = uz + vth * randn(rng, T)
    end

    # --- momentum correction: shift so subset momentum = pre-scatter value ----
    # Desired subset weighted mean (restores Σ_S w v exactly):
    ūx = Px / Wsub
    ūy = Py / Wsub
    ūz = Pz / Wsub
    # Current post-resample subset weighted mean:
    nsx = zero(T)
    nsy = zero(T)
    nsz = zero(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p]
        nsx += wp * vx[p]
        nsy += wp * vy[p]
        nsz += wp * vz[p]
    end
    n̄x = nsx / Wsub
    n̄y = nsy / Wsub
    n̄z = nsz / Wsub
    δx = ūx - n̄x
    δy = ūy - n̄y
    δz = ūz - n̄z
    @inbounds for p = 1:N
        sel[p] || continue
        vx[p] += δx
        vy[p] += δy
        vz[p] += δz
    end
    # Now the subset weighted mean is (ūx,ūy,ūz) → subset momentum = (Px,Py,Pz).

    # --- energy correction: rescale fluctuations about the (now-exact) mean ----
    # Peculiar (mean-removed) kinetic energy currently in the subset:
    K1 = zero(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p]
        dx = vx[p] - ūx
        dy = vy[p] - ūy
        dz = vz[p] - ūz
        K1 += wp * (dx * dx + dy * dy + dz * dz)
    end
    # Target peculiar energy = pre-scatter total energy minus mean-flow energy.
    # Eold = Wsub|ū|² + K_target  (exact identity for the pre-scatter subset),
    # and Wsub|ū|² is exactly the mean-flow energy of the corrected subset, so
    # K_target ≥ 0 by construction (Eold ≥ Wsub|ū|² by Cauchy–Schwarz).
    Emean = Wsub * (ūx^2 + ūy^2 + ūz^2)
    Ktarget = Eold - Emean
    if K1 > 0 && Ktarget > 0
        α = sqrt(Ktarget / K1)
        @inbounds for p = 1:N
            sel[p] || continue
            vx[p] = ūx + α * (vx[p] - ūx)
            vy[p] = ūy + α * (vy[p] - ūy)
            vz[p] = ūz + α * (vz[p] - ūz)
        end
    end
    # Rescaling about ū leaves the subset mean (hence momentum) unchanged and sets
    # subset energy to Wsub|ū|² + Ktarget = Eold. Both subset totals restored ⇒
    # whole-set Σ w v and Σ w |v|² conserved exactly (to roundoff).

    return ps
end

# ---------------------------------------------------------------- Coulomb (Takizuka-Abe)

"""
    collide_coulomb!(ps::ParticleSet{D,T}, gcoeff, dt; rng=Random.default_rng(),
                     u_floor=1e-3) -> ps

One **Takizuka-Abe (1977)** binary-Coulomb collision substep. The particles are
randomly paired; each pair's relative velocity `g = v_i − v_j` is rotated by a polar
scattering angle `Θ` (azimuth `Φ` uniform) with `δ = tan(Θ/2)` drawn from
`𝒩(0, ⟨δ²⟩)`, `⟨δ²⟩ = gcoeff·dt / max(|g|, u_floor)³`. The `|g|⁻³` velocity dependence
is the physical Coulomb law (slow particles scatter through larger angles); `gcoeff`
bundles the collisional prefactor `n·lnΛ·(qᵢqⱼ/…)²` in the code's normalized units, so
it sets the collisionality. Each pair is updated about its **weighted centre of mass**
`V = (wᵢvᵢ+wⱼvⱼ)/(wᵢ+wⱼ)`:

    v_i ← V + (wⱼ/(wᵢ+wⱼ)) g',   v_j ← V − (wᵢ/(wᵢ+wⱼ)) g',   g' = R(Θ,Φ) g

Because `|g'| = |g|` (a pure rotation) this conserves each pair's momentum **and**
energy exactly for arbitrary weights (to roundoff), hence the whole-set `Σwv` and
`Σw|v|²`. Unlike BGK ([`collide_bgk!`]) this reproduces true Coulomb velocity-space
diffusion and relaxes a temperature **anisotropy** toward isotropy at the physical,
speed-dependent rate.

Whole-set pairing (0-D velocity-space relaxation, same scope as `collide_bgk!`);
cell-local pairing is the spatially-resolved upgrade. `gcoeff ≥ 0`, `dt ≥ 0`;
fewer than two particles is a no-op. Returns `ps`.
"""
function collide_coulomb!(
    ps::ParticleSet{D,T},
    gcoeff::Real,
    dt::Real;
    rng = Random.default_rng(),
    u_floor::Real = 1e-3,
) where {D,T}
    gcoeff >= 0 || throw(ArgumentError("collision coefficient gcoeff must be ≥ 0"))
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    N = nparticles(ps)
    (N < 2 || gcoeff == 0 || dt == 0) && return ps

    vx, vy, vz = ps.v
    w = ps.weight
    gc = T(gcoeff)
    dtT = T(dt)
    twoπ = 2 * T(π)

    idx = randperm(rng, N)                     # random pairing (like collide_bgk!'s falses(N))
    npair = N ÷ 2
    @inbounds for k = 1:npair
        i = idx[2k-1]
        j = idx[2k]
        wi = w[i]
        wj = w[j]
        wsum = wi + wj
        wsum > 0 || continue

        gx = vx[i] - vx[j]
        gy = vy[i] - vy[j]
        gz = vz[i] - vz[j]
        gmag = sqrt(gx * gx + gy * gy + gz * gz)
        gmag > 0 || continue                   # identical velocities: no relative motion

        Vx = (wi * vx[i] + wj * vx[j]) / wsum
        Vy = (wi * vy[i] + wj * vy[j]) / wsum
        Vz = (wi * vz[i] + wj * vz[j]) / wsum

        var = gc * dtT / max(gmag, uf)^3       # Takizuka-Abe ⟨δ²⟩
        δ = sqrt(var) * randn(rng, T)
        δ2 = δ * δ
        cosθ = (one(T) - δ2) / (one(T) + δ2)   # δ = tan(Θ/2) ⇒ cosΘ,sinΘ
        sinθ = 2δ / (one(T) + δ2)
        φ = twoπ * rand(rng, T)
        cosφ = cos(φ)
        sinφ = sin(φ)

        gperp = sqrt(gx * gx + gy * gy)
        if gperp > 0                            # rotate g by (Θ,Φ) — Takizuka & Abe 1977 eq.
            Δgx =
                (gx / gperp) * gz * sinθ * cosφ - (gy / gperp) * gmag * sinθ * sinφ -
                gx * (one(T) - cosθ)
            Δgy =
                (gy / gperp) * gz * sinθ * cosφ + (gx / gperp) * gmag * sinθ * sinφ -
                gy * (one(T) - cosθ)
            Δgz = -gperp * sinθ * cosφ - gz * (one(T) - cosθ)
        else                                    # g along ±z: rotate the z-aligned vector directly
            Δgx = gmag * sinθ * cosφ
            Δgy = gmag * sinθ * sinφ
            Δgz = -gz * (one(T) - cosθ)          # gz = ±gmag ⇒ |g'| = gmag preserved
        end

        gpx = gx + Δgx
        gpy = gy + Δgy
        gpz = gz + Δgz
        fi = wj / wsum
        fj = wi / wsum
        vx[i] = Vx + fi * gpx
        vy[i] = Vy + fi * gpy
        vz[i] = Vz + fi * gpz
        vx[j] = Vx - fj * gpx
        vy[j] = Vy - fj * gpy
        vz[j] = Vz - fj * gpz
    end
    return ps
end

# ---------------------------------------------------------------- neutral MCC (elastic)

"""
    collide_neutral_mcc!(ps::ParticleSet{D,T}, dt; nσ, T_n, m_n=1.0,
                         u_n=(0.0,0.0,0.0), rng=Random.default_rng()) -> ps

One **Monte-Carlo-collision (MCC)** substep of elastic scattering off a background
**neutral gas** — a thermal reservoir at temperature `T_n`, mass `m_n`, bulk drift
`u_n`, and density×cross-section product `nσ = n_n·σ`. For each charged particle
(mass `m_p = ps.m`) a neutral partner velocity `v_n` is drawn from the neutral
Maxwellian; with probability `P = 1 − exp(−nσ·|g|·dt)` (relative speed `|g| = |v−v_n|`)
the pair scatters **elastically and isotropically in the centre-of-mass frame**:

    V = (m_p v + m_n v_n)/(m_p+m_n),   g' = |g| n̂  (n̂ isotropic),   v ← V + (m_n/(m_p+m_n)) g'

Each binary collision conserves the (particle+neutral) momentum and energy exactly,
but the neutral is a **reservoir** (freshly sampled and discarded each collision), so
the charged population relaxes toward the neutral distribution: its temperature toward
`T_n` and its drift toward `u_n` (full thermalization when `m_p = m_n`).

`nσ ≥ 0`, `T_n ≥ 0`, `dt ≥ 0`, `m_n > 0`. Elastic only (inelastic excitation and
ionization are the upgrade path). Returns `ps`.
"""
function collide_neutral_mcc!(
    ps::ParticleSet{D,T},
    dt::Real;
    nσ::Real,
    T_n::Real,
    m_n::Real = 1.0,
    u_n::NTuple{3,<:Real} = (0.0, 0.0, 0.0),
    rng = Random.default_rng(),
) where {D,T}
    nσ >= 0 || throw(ArgumentError("nσ (density×cross-section) must be ≥ 0"))
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    TnT = _require_finite_nonnegative_real("T_n", T_n, T)
    mnT = _require_finite_positive_real("m_n", m_n, T)
    N = nparticles(ps)
    (N == 0 || nσ == 0 || dt == 0) && return ps

    vx, vy, vz = ps.v
    mp = T(ps.m)
    nσT = T(nσ)
    dtT = T(dt)
    vthn = sqrt(TnT / mnT)                          # neutral thermal speed √(T_n/m_n)
    unx, uny, unz = T(u_n[1]), T(u_n[2]), T(u_n[3])
    invM = one(T) / (mp + mnT)
    μn = mnT * invM                                 # v ← V + μn g'
    twoπ = 2 * T(π)
    @inbounds for p = 1:N
        vnx = unx + vthn * randn(rng, T)            # sample a neutral partner
        vny = uny + vthn * randn(rng, T)
        vnz = unz + vthn * randn(rng, T)
        gx = vx[p] - vnx
        gy = vy[p] - vny
        gz = vz[p] - vnz
        gmag = sqrt(gx * gx + gy * gy + gz * gz)
        gmag > 0 || continue
        Pcoll = -expm1(-nσT * gmag * dtT)           # 1 − exp(−nσ|g|dt)
        rand(rng, T) < Pcoll || continue
        Vx = (mp * vx[p] + mnT * vnx) * invM        # centre-of-mass velocity
        Vy = (mp * vy[p] + mnT * vny) * invM
        Vz = (mp * vz[p] + mnT * vnz) * invM
        cosχ = 2 * rand(rng, T) - one(T)            # isotropic elastic scatter in CM
        sinχ = sqrt(max(zero(T), one(T) - cosχ * cosχ))
        φ = twoπ * rand(rng, T)
        vx[p] = Vx + μn * gmag * sinχ * cos(φ)
        vy[p] = Vy + μn * gmag * sinχ * sin(φ)
        vz[p] = Vz + μn * gmag * cosχ
    end
    return ps
end
