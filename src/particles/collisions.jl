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

@inline function _scaled_centered_component(rms_scale::T, value::T, center::T) where {T}
    difference = value - center
    isfinite(difference) && return rms_scale * difference
    return rms_scale * value - rms_scale * center
end

@inline function _scaled_ssq_add(scale::T, ssq::T, value::T) where {T}
    isfinite(value) || return T(Inf), one(T)
    iszero(value) && return scale, ssq
    magnitude = abs(value)
    if scale < magnitude
        ratio = scale / magnitude
        return magnitude, one(T) + ssq * ratio * ratio
    end
    ratio = magnitude / scale
    return scale, ssq + ratio * ratio
end

@inline function _accumulate_normalized_deviation(
    scale::T,
    ssq::T,
    wp::T,
    sqrt_total_weight::T,
    x::T,
    y::T,
    z::T,
    center_x::T,
    center_y::T,
    center_z::T,
) where {T}
    rms_scale = sqrt(wp) / sqrt_total_weight
    dx = _scaled_centered_component(rms_scale, x, center_x)
    dy = _scaled_centered_component(rms_scale, y, center_y)
    dz = _scaled_centered_component(rms_scale, z, center_z)
    scale, ssq = _scaled_ssq_add(scale, ssq, dx)
    scale, ssq = _scaled_ssq_add(scale, ssq, dy)
    return _scaled_ssq_add(scale, ssq, dz)
end

@inline function _scaled_ssq_energy(scale::T, ssq::T) where {T}
    iszero(scale) && return zero(T)
    norm = scale * sqrt(ssq)
    return norm * norm
end

@inline function _finite_muldiv(a::T, b::T, denominator::T) where {T}
    (iszero(a) || iszero(b)) && return zero(T)
    af, ae = frexp(a)
    bf, be = frexp(b)
    df, de = frexp(denominator)
    return ldexp((af * bf) / df, ae + be - de)
end

@inline function _convex_weighted_mean(
    current::T,
    value::T,
    old_weight::T,
    new_weight::T,
    total_weight::T,
) where {T}
    current == value && return current
    old_fraction = old_weight / total_weight
    new_fraction = new_weight / total_weight
    old_term =
        iszero(old_fraction) ? _finite_muldiv(old_weight, current, total_weight) :
        old_fraction * current
    new_term =
        iszero(new_fraction) ? _finite_muldiv(new_weight, value, total_weight) :
        new_fraction * value
    return old_term + new_term
end

"""
    collide_bgk!(ps::ParticleSet{D,T}, ν, dt; rng=Random.default_rng(),
                 work=nothing, selection_work=nothing) -> ps

Apply one BGK collision substep of duration `dt` at collision frequency `ν`.

The velocity distribution relaxes toward the local drifting Maxwellian: a random
fraction `1 − exp(−ν·dt)` of the particles is resampled from an isotropic
Maxwellian with the set's weighted mean velocity `u` and scalar temperature
`T = Σ w |v−u|² / (3 Σ w)` (so `v_c = u_c + √T·𝒩(0,1)`). The scattered subset is
then corrected (a uniform velocity shift to restore its momentum, followed by an
isotropic rescaling of its velocity fluctuations to restore its energy) so that
the whole-set totals `Σ w v` and `Σ w |v|²` are conserved exactly to roundoff.

Finite `ν ≥ 0`, finite `dt ≥ 0`, finite velocities, and finite non-negative
macro-particle weights are required. `ν·dt = 0`, fewer than two particles, or
an all-zero-weight set is a no-op. Every positive weight must remain nonzero
when normalized by the largest weight. The sampled update is transactional: an
invalid sampled state throws without changing `ps`, while an affine correction
that is not representable leaves the selected state unchanged. By default this
uses three temporary velocity vectors; repeated callers can pass a non-aliasing
3-tuple of particle-length vectors as `work` and a particle-length
`AbstractVector{Bool}` as `selection_work`. Returns `ps`.
"""
function collide_bgk!(
    ps::ParticleSet{D,T},
    ν::Real,
    dt::Real;
    rng = Random.default_rng(),
    work = nothing,
    selection_work = nothing,
) where {D,T}
    νT = _require_finite_nonnegative_real("collision frequency ν", ν, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "collide_bgk!")
    _require_finite_particle_velocities(ps, "collide_bgk!")
    N = nparticles(ps)
    vx, vy, vz = ps.v
    w = ps.weight

    # Normalize weights once so a uniform, physically irrelevant macro-weight
    # scale cannot overflow the conserved-moment arithmetic.
    wmax = zero(T)
    @inbounds for p = 1:N
        wp = w[p]
        (isfinite(wp) && wp >= zero(T)) ||
            throw(ArgumentError("collide_bgk!: particle $p weight must be finite and non-negative"))
        wmax = max(wmax, wp)
    end
    # Need at least two particles to have a fluctuation to relax/correct.
    (N < 2 || iszero(νT) || iszero(dtT) || iszero(wmax)) && return ps
    _require_nonaliasing_collision_scratch(selection_work, "selection_work", work, "collide_bgk!")

    # --- whole-set weighted mean velocity u and scalar temperature T ----------
    Wtot = zero(T)
    ux = zero(T)
    uy = zero(T)
    uz = zero(T)
    @inbounds for p = 1:N
        wp = w[p] / wmax
        (w[p] > zero(T) && iszero(wp)) && throw(
            ArgumentError(
                "collide_bgk!: particle $p positive weight is below the normalized numeric range relative to the largest weight",
            ),
        )
        iszero(wp) && continue
        Wnew = Wtot + wp
        ux = _convex_weighted_mean(ux, vx[p], Wtot, wp, Wnew)
        uy = _convex_weighted_mean(uy, vy[p], Wtot, wp, Wnew)
        uz = _convex_weighted_mean(uz, vz[p], Wtot, wp, Wnew)
        Wtot = Wnew
    end
    (isfinite(ux) && isfinite(uy) && isfinite(uz)) ||
        throw(ArgumentError("collide_bgk!: weighted mean velocity exceeds numeric range"))

    thermal_scale = zero(T)
    thermal_ssq = one(T)
    sqrt_Wtot = sqrt(Wtot)
    @inbounds for p = 1:N
        wp = w[p] / wmax
        iszero(wp) && continue
        thermal_scale, thermal_ssq = _accumulate_normalized_deviation(
            thermal_scale,
            thermal_ssq,
            wp,
            sqrt_Wtot,
            vx[p],
            vy[p],
            vz[p],
            ux,
            uy,
            uz,
        )
    end
    s2 = _scaled_ssq_energy(thermal_scale, thermal_ssq)
    isfinite(s2) || throw(ArgumentError("collide_bgk!: thermal energy exceeds numeric range"))
    # scalar temperature (per-component variance) → isotropic thermal speed
    Tscalar = s2 / 3
    vth = sqrt(Tscalar)

    # --- select the scattered subset (Bernoulli with p = 1 − exp(−ν dt)) ------
    Pcoll = -expm1(-νT * dtT)                    # = 1 − exp(−ν dt), accurate near 0
    sel = _collision_selection_work(ps, selection_work, "collide_bgk!")
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

    # --- record the subset's normalized weight, mean, and thermal energy ------
    Wsub = zero(T)
    ūx = zero(T)
    ūy = zero(T)
    ūz = zero(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p] / wmax
        iszero(wp) && continue
        Wnew = Wsub + wp
        ūx = _convex_weighted_mean(ūx, vx[p], Wsub, wp, Wnew)
        ūy = _convex_weighted_mean(ūy, vy[p], Wsub, wp, Wnew)
        ūz = _convex_weighted_mean(ūz, vz[p], Wsub, wp, Wnew)
        Wsub = Wnew
    end
    Wsub > 0 || return ps
    target_scale = zero(T)
    target_ssq = one(T)
    sqrt_Wsub = sqrt(Wsub)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p] / wmax
        iszero(wp) && continue
        target_scale, target_ssq = _accumulate_normalized_deviation(
            target_scale,
            target_ssq,
            wp,
            sqrt_Wsub,
            vx[p],
            vy[p],
            vz[p],
            ūx,
            ūy,
            ūz,
        )
    end
    Ktarget = _scaled_ssq_energy(target_scale, target_ssq)
    isfinite(Ktarget) ||
        throw(ArgumentError("collide_bgk!: selected thermal energy exceeds numeric range"))

    # --- resample the subset from the local drifting Maxwellian ---------------
    vnew = _collision_velocity_work(ps, work, "collide_bgk!")
    nvx, nvy, nvz = vnew
    @inbounds for p = 1:N
        sel[p] || continue
        sx = ux + vth * randn(rng, T)
        sy = uy + vth * randn(rng, T)
        sz = uz + vth * randn(rng, T)
        (isfinite(sx) && isfinite(sy) && isfinite(sz)) ||
            throw(ArgumentError("collide_bgk!: sampled velocity must be finite"))
        nvx[p] = sx
        nvy[p] = sy
        nvz[p] = sz
    end

    # --- affine correction: restore the subset's momentum and thermal energy ---
    # Current post-resample subset weighted mean:
    n̄x = zero(T)
    n̄y = zero(T)
    n̄z = zero(T)
    Wnewsub = zero(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p] / wmax
        iszero(wp) && continue
        Wnew = Wnewsub + wp
        n̄x = _convex_weighted_mean(n̄x, nvx[p], Wnewsub, wp, Wnew)
        n̄y = _convex_weighted_mean(n̄y, nvy[p], Wnewsub, wp, Wnew)
        n̄z = _convex_weighted_mean(n̄z, nvz[p], Wnewsub, wp, Wnew)
        Wnewsub = Wnew
    end
    # Measure the original sampled deviations. Applying a rounded shift first can
    # turn identical samples into a common nonzero residual that is then amplified.
    sampled_scale = zero(T)
    sampled_ssq = one(T)
    @inbounds for p = 1:N
        sel[p] || continue
        wp = w[p] / wmax
        iszero(wp) && continue
        sampled_scale, sampled_ssq = _accumulate_normalized_deviation(
            sampled_scale,
            sampled_ssq,
            wp,
            sqrt_Wsub,
            nvx[p],
            nvy[p],
            nvz[p],
            n̄x,
            n̄y,
            n̄z,
        )
    end
    K1 = _scaled_ssq_energy(sampled_scale, sampled_ssq)
    isfinite(K1) ||
        throw(ArgumentError("collide_bgk!: corrected thermal energy exceeds numeric range"))
    restore_original = false
    if K1 > 0
        α = sqrt(Ktarget) / sqrt(K1)
        restore_original = !isfinite(α)
        if !restore_original
            @inbounds for p = 1:N
                sel[p] || continue
                corrected_x = ūx + α * (nvx[p] - n̄x)
                corrected_y = ūy + α * (nvy[p] - n̄y)
                corrected_z = ūz + α * (nvz[p] - n̄z)
                if !(isfinite(corrected_x) && isfinite(corrected_y) && isfinite(corrected_z))
                    restore_original = true
                    break
                end
                nvx[p] = corrected_x
                nvy[p] = corrected_y
                nvz[p] = corrected_z
            end
        end
    elseif Ktarget > 0
        # A degenerate random source (or finite-precision collapse) can give every
        # resampled particle the same velocity, so there is no fluctuation to
        # rescale. Preserve the original selected state transactionally instead
        # of synthesizing deviations that may not be representable at this scale.
        restore_original = true
    else
        @inbounds for p = 1:N
            sel[p] || continue
            nvx[p] = ūx
            nvy[p] = ūy
            nvz[p] = ūz
        end
    end
    if restore_original
        @inbounds for p = 1:N
            sel[p] || continue
            nvx[p] = vx[p]
            nvy[p] = vy[p]
            nvz[p] = vz[p]
        end
    end

    @inbounds for p = 1:N
        sel[p] || continue
        (isfinite(nvx[p]) && isfinite(nvy[p]) && isfinite(nvz[p])) ||
            throw(ArgumentError("collide_bgk!: corrected velocity must be finite"))
    end
    for c = 1:3
        copyto!(ps.v[c], vnew[c])
    end
    return ps
end

# ---------------------------------------------------------------- Coulomb (Takizuka-Abe)

# One Takizuka-Abe pair scatter with variance ⟨δ²⟩ = gcdt·(max(wᵢ,wⱼ)/w̄) / max(|g|, uf)³.
# All particles of a set share the ONE physical mass ps.m, so the physical kick is the
# symmetric half-kick about the true midpoint (vᵢ+vⱼ)/2 regardless of macro-particle
# weight. Unequal weights use the CORRECTED rejection scheme (Higginson, Holod & Link,
# JCP 413 (2020) 109450, eqs. 31-32, fixing Nanbu & Yonemura 1998 / Sentoku & Kemp 2008):
# each partner accepts its update with probability w_other/max(wᵢ,wⱼ) AND the pair
# variance is scaled by max(wᵢ,wⱼ)/w̄ (w̄ = mean set weight), so every particle's expected
# per-step variance is E_partner[(w_other/wmax)·(wmax/w̄)]·⟨δ²⟩ = ⟨δ²⟩ — the physical
# relaxation rate independent of the weight distribution — while momentum and energy are
# conserved in expectation (exactly, per pair, for equal weights: both accept, factor 1).
# When ⟨δ²⟩ > 1 the Gaussian tan(Θ/2) model is invalid (the sampled Θ → π, g' ≈ −g is a
# statistical no-op, so relaxation would STALL as collisionality rises); the pair then
# draws an isotropic scattering angle instead (strong-collision limit, Nanbu 1997 /
# Pérez et al. 2012), so relaxation saturates at the collision rate.
@inline function _coulomb_pair_scatter!(
    vx,
    vy,
    vz,
    w,
    i::Int,
    j::Int,
    gcdt::T,
    log_gcdt::T,
    wbar::T,
    log_wbar::T,
    weight_scale::T,
    uf::T,
    twoπ::T,
    rng,
) where {T}
    @inbounds begin
        wi = w[i] / weight_scale
        wj = w[j] / weight_scale
        wmax = max(wi, wj)
        wmax > 0 || return nothing

        gx = vx[i] - vx[j]
        gy = vy[i] - vy[j]
        gz = vz[i] - vz[j]
        (isfinite(gx) && isfinite(gy) && isfinite(gz)) ||
            throw(ArgumentError("collide_coulomb!: relative velocity exceeds numeric range"))
        gmag = hypot(gx, gy, gz)
        isfinite(gmag) ||
            throw(ArgumentError("collide_coulomb!: relative speed exceeds numeric range"))
        gmag > 0 || return nothing             # identical velocities: no relative motion

        # Divide before adding so a representable same-mass midpoint does not
        # overflow when both velocities are large and have the same sign.
        Mx = vx[i] / 2 + vx[j] / 2
        My = vy[i] / 2 + vy[j] / 2
        Mz = vz[i] / 2 + vz[j] / 2

        # Takizuka-Abe ⟨δ²⟩, weight-corrected: × max(wᵢ,wⱼ)/w̄ (Higginson 2020 eq. 31).
        # The ordinary-range path avoids transcendental work in this hot loop.
        # Fall back to logs only when an intermediate product over/underflows.
        speed = max(gmag, uf)
        numerator = gcdt * (wmax / wbar)
        speed2 = speed * speed
        denominator = speed2 * speed
        use_fast =
            isfinite(numerator) &&
            numerator > zero(T) &&
            isfinite(denominator) &&
            denominator > zero(T)
        if use_fast
            var = numerator / denominator
            large_angle = !isfinite(var) || var > one(T)
        else
            logvar = log_gcdt + log(wmax) - log_wbar - 3 * log(speed)
            large_angle = logvar > zero(T)
            var = large_angle ? one(T) : exp(logvar)
        end
        if large_angle                           # large-angle regime → isotropic scatter
            cosθ = 2 * rand(rng, T) - one(T)
            sinθ = sqrt(max(zero(T), one(T) - cosθ * cosθ))
        else
            δ = sqrt(var) * randn(rng, T)
            isfinite(δ) ||
                throw(ArgumentError("collide_coulomb!: sampled scattering angle must be finite"))
            if abs(δ) <= one(T)
                δ2 = δ * δ
                denom = one(T) + δ2
                cosθ = (one(T) - δ2) / denom
                sinθ = 2δ / denom
            else
                invδ = one(T) / δ
                invδ2 = invδ * invδ
                denom = one(T) + invδ2
                cosθ = (invδ2 - one(T)) / denom
                sinθ = 2invδ / denom
            end
        end
        φ = twoπ * rand(rng, T)
        (isfinite(cosθ) && isfinite(sinθ) && isfinite(φ)) ||
            throw(ArgumentError("collide_coulomb!: sampled scattering direction must be finite"))
        sinφ, cosφ = sincos(φ)

        # Rotate the unit relative-velocity vector, renormalize against roundoff,
        # and only then restore its magnitude. Direct O(gmag) sums can overflow
        # transiently even though every component of the result is bounded.
        invgmag = inv(gmag)
        if isfinite(invgmag)
            ugx = gx * invgmag
            ugy = gy * invgmag
            ugz = gz * invgmag
        else
            ugx = gx / gmag
            ugy = gy / gmag
            ugz = gz / gmag
        end
        ugperp = hypot(ugx, ugy)
        if ugperp > 0                           # Takizuka & Abe 1977 rotation
            invugperp = inv(ugperp)
            if isfinite(invugperp)
                xperp = ugx * invugperp
                yperp = ugy * invugperp
            else
                xperp = ugx / ugperp
                yperp = ugy / ugperp
            end
            polar_term = ugz * sinθ * cosφ
            azimuthal_term = sinθ * sinφ
            ugpx = ugx * cosθ + xperp * polar_term - yperp * azimuthal_term
            ugpy = ugy * cosθ + yperp * polar_term + xperp * azimuthal_term
            ugpz = ugz * cosθ - ugperp * sinθ * cosφ
        else                                    # g along ±z
            ugpx = sinθ * cosφ
            ugpy = sinθ * sinφ
            ugpz = ugz * cosθ
        end
        ugpmag = hypot(ugpx, ugpy, ugpz)
        (isfinite(ugpmag) && ugpmag > zero(T)) ||
            throw(ArgumentError("collide_coulomb!: rotated relative direction must be finite"))
        invugpmag = inv(ugpmag)
        gpx = gmag * (ugpx * invugpmag)
        gpy = gmag * (ugpy * invugpmag)
        gpz = gmag * (ugpz * invugpmag)
        (isfinite(gpx) && isfinite(gpy) && isfinite(gpz)) ||
            throw(ArgumentError("collide_coulomb!: rotated relative velocity must be finite"))
        # Rejection step: the max-weight partner is probabilistically skipped (equal
        # weights short-circuit — no draw, both always accept).
        accept_i = wi == wj || rand(rng, T) < wj / wmax
        accept_j = wi == wj || rand(rng, T) < wi / wmax
        vix = accept_i ? Mx + gpx / 2 : vx[i]
        viy = accept_i ? My + gpy / 2 : vy[i]
        viz = accept_i ? Mz + gpz / 2 : vz[i]
        vjx = accept_j ? Mx - gpx / 2 : vx[j]
        vjy = accept_j ? My - gpy / 2 : vy[j]
        vjz = accept_j ? Mz - gpz / 2 : vz[j]
        (
            isfinite(vix) &&
            isfinite(viy) &&
            isfinite(viz) &&
            isfinite(vjx) &&
            isfinite(vjy) &&
            isfinite(vjz)
        ) || throw(ArgumentError("collide_coulomb!: scattered velocity must be finite"))
        vx[i] = vix
        vy[i] = viy
        vz[i] = viz
        vx[j] = vjx
        vy[j] = vjy
        vz[j] = vjz
    end
    return nothing
end

function _collision_permutation(rng, ps::ParticleSet, index_work, context::AbstractString)
    N = nparticles(ps)
    index_work === nothing && return randperm(rng, N)
    index_work isa AbstractVector{Int} ||
        throw(ArgumentError("$context: index_work must be an AbstractVector{Int}"))
    axes(index_work, 1) == Base.OneTo(N) || throw(
        DimensionMismatch(
            "$context: index_work axes $(axes(index_work, 1)) do not match particle axes $(Base.OneTo(N))",
        ),
    )
    _mightalias_particle_storage(index_work, ps) &&
        throw(ArgumentError("$context: index_work must not alias particle storage"))
    @inbounds for p = 1:N
        index_work[p] = p
    end
    shuffle!(rng, index_work)
    return index_work
end

"""
    collide_coulomb!(ps::ParticleSet{D,T}, gcoeff, dt; rng=Random.default_rng(),
                     u_floor=1e-3, work=nothing, index_work=nothing) -> ps

One **Takizuka-Abe (1977)** binary-Coulomb collision substep. The particles are
randomly paired; each pair's relative velocity `g = v_i − v_j` is rotated by a polar
scattering angle `Θ` (azimuth `Φ` uniform) with `δ = tan(Θ/2)` drawn from
`𝒩(0, ⟨δ²⟩)`, `⟨δ²⟩ = gcoeff·dt / max(|g|, u_floor)³`. The `|g|⁻³` velocity dependence
is the physical Coulomb law (slow particles scatter through larger angles); `gcoeff`
bundles the collisional prefactor `n·lnΛ·(qᵢqⱼ/…)²` in the code's normalized units, so
it sets the collisionality. When `⟨δ²⟩ > 1` the small-angle Gaussian model is invalid
(`Θ → π` degenerates into a velocity swap), so the pair draws an **isotropic**
scattering angle instead (`cosΘ` uniform on [−1,1] — the strong-collision limit, Nanbu
1997 / Pérez et al. 2012): relaxation saturates with collisionality instead of stalling.

All particles of a set share the one physical mass `ps.m`, so each partner is kicked
**symmetrically about the true midpoint** `M = (vᵢ+vⱼ)/2`:

    v_i ← M + g'/2,   v_j ← M − g'/2,   g' = R(Θ,Φ) g

For **unequal macro-particle weights** the update is applied via the corrected rejection
scheme (Higginson, Holod & Link, JCP 413 (2020) 109450, fixing Nanbu & Yonemura 1998 /
Sentoku & Kemp 2008): the pair's `⟨δ²⟩` is scaled by `max(wᵢ,wⱼ)/w̄` (`w̄` = mean set
weight), particle `i` accepts its update with probability `wⱼ/max(wᵢ,wⱼ)` and `j` with
`wᵢ/max(wᵢ,wⱼ)`. Every particle's expected per-step scattering variance is then the
nominal `⟨δ²⟩` — the physical relaxation rate independent of the weight distribution —
and `Σwv` and `Σw|v|²` are conserved **in expectation**. Equal weights accept both
updates always (variance factor ≡ 1), and `|g'| = |g|` (a pure rotation) then conserves
each pair's momentum **and** energy exactly (to roundoff), hence the whole-set totals.
Unlike BGK ([`collide_bgk!`]) this reproduces true Coulomb velocity-space diffusion and
relaxes a temperature **anisotropy** toward isotropy at the physical, speed-dependent
rate.

Whole-set pairing (0-D velocity-space relaxation, same scope as `collide_bgk!`);
cell-local pairing is the spatially-resolved upgrade. Odd `N` uses the Takizuka-Abe
triplet: the first three particles of the permutation form the pairs (1,2),(2,3),(3,1),
each scattered with **half** the variance (small-angle variances add linearly), so every
particle collides every step at the nominal rate. `gcoeff ≥ 0`, `dt ≥ 0`; fewer than
two particles or an all-zero-weight set is a no-op. Velocities and macro-particle
weights must be finite and non-negative, and every positive weight must remain
nonzero when normalized by the largest weight. Updates are staged transactionally
in three temporary velocity vectors; repeated callers can pass non-aliasing
particle-length vectors as `work`. A supplied particle-length
`AbstractVector{Int}` in `index_work` also reuses the pairing permutation storage.
Returns `ps`.
"""
function collide_coulomb!(
    ps::ParticleSet{D,T},
    gcoeff::Real,
    dt::Real;
    rng = Random.default_rng(),
    u_floor::Real = 1e-3,
    work = nothing,
    index_work = nothing,
) where {D,T}
    gc = _require_finite_nonnegative_real("collision coefficient gcoeff", gcoeff, T)
    dtT = _require_finite_nonnegative_real("dt", dt, T)
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    _require_finite_particle_velocities(ps, "collide_coulomb!")
    N = nparticles(ps)
    w = ps.weight
    weight_scale = zero(T)
    @inbounds for p = 1:N
        wp = w[p]
        (isfinite(wp) && wp >= zero(T)) || throw(
            ArgumentError("collide_coulomb!: particle $p weight must be finite and non-negative"),
        )
        weight_scale = max(weight_scale, wp)
    end
    (N < 2 || iszero(gc) || iszero(dtT) || iszero(weight_scale)) && return ps
    _require_nonaliasing_collision_scratch(index_work, "index_work", work, "collide_coulomb!")
    gcdt = gc * dtT
    log_gcdt = log(gc) + log(dtT)

    # Mean macro-particle weight w̄ normalizes the per-pair variance factor max(wᵢ,wⱼ)/w̄
    # (Higginson 2020): the scheme is then invariant under re-partitioning the same
    # physical plasma into different macro-weights at fixed gcoeff, and reduces exactly
    # to plain Takizuka-Abe for uniform weights (factor ≡ 1).
    Wtot = zero(T)
    @inbounds for p = 1:N
        normalized_weight = w[p] / weight_scale
        (w[p] > zero(T) && iszero(normalized_weight)) && throw(
            ArgumentError(
                "collide_coulomb!: particle $p positive weight is below the normalized numeric range relative to the largest weight",
            ),
        )
        Wtot += normalized_weight
    end
    isfinite(Wtot) ||
        throw(ArgumentError("collide_coulomb!: normalized total weight exceeds numeric range"))
    wbar = Wtot / T(N)
    log_wbar = log(wbar)

    vnew = _collision_velocity_work(ps, work, "collide_coulomb!")
    vx, vy, vz = vnew
    idx = _collision_permutation(rng, ps, index_work, "collide_coulomb!")
    twoπ = 2 * T(π)
    base = 0
    if isodd(N)
        # Takizuka-Abe odd-N prescription: the first three particles of the permutation
        # form the pairs (1,2),(2,3),(3,1), each scattered with HALF the variance (the
        # small-angle variances add linearly, so each triplet particle still accumulates
        # the full per-step ⟨δ²⟩), and every particle collides every step.
        i1 = idx[1]
        i2 = idx[2]
        i3 = idx[3]
        halfgcdt = gcdt / 2
        half_log_gcdt = log_gcdt - log(T(2))
        _coulomb_pair_scatter!(
            vx,
            vy,
            vz,
            w,
            i1,
            i2,
            halfgcdt,
            half_log_gcdt,
            wbar,
            log_wbar,
            weight_scale,
            uf,
            twoπ,
            rng,
        )
        _coulomb_pair_scatter!(
            vx,
            vy,
            vz,
            w,
            i2,
            i3,
            halfgcdt,
            half_log_gcdt,
            wbar,
            log_wbar,
            weight_scale,
            uf,
            twoπ,
            rng,
        )
        _coulomb_pair_scatter!(
            vx,
            vy,
            vz,
            w,
            i3,
            i1,
            halfgcdt,
            half_log_gcdt,
            wbar,
            log_wbar,
            weight_scale,
            uf,
            twoπ,
            rng,
        )
        base = 3
    end
    npair = (N - base) ÷ 2
    @inbounds for k = 1:npair
        _coulomb_pair_scatter!(
            vx,
            vy,
            vz,
            w,
            idx[base+2k-1],
            idx[base+2k],
            gcdt,
            log_gcdt,
            wbar,
            log_wbar,
            weight_scale,
            uf,
            twoπ,
            rng,
        )
    end
    for c = 1:3
        copyto!(ps.v[c], vnew[c])
    end
    return ps
end

# ---------------------------------------------------------------- neutral MCC (elastic)

function _require_finite_particle_velocities(ps::ParticleSet, context::AbstractString)
    @inbounds for p in eachindex(ps.weight)
        for c = 1:3
            isfinite(ps.v[c][p]) ||
                throw(ArgumentError("$context: particle $p velocity component $c must be finite"))
        end
    end
    return nothing
end

function _mightalias_particle_storage(candidate, ps::ParticleSet)
    for d = 1:length(ps.x)
        Base.mightalias(candidate, ps.x[d]) && return true
    end
    for c = 1:3
        Base.mightalias(candidate, ps.v[c]) && return true
    end
    Base.mightalias(candidate, ps.weight) && return true
    Base.mightalias(candidate, ps.id) && return true
    Base.mightalias(candidate, ps.tag) && return true
    return false
end

function _require_nonaliasing_collision_scratch(candidate, candidate_name, work, context)
    (candidate === nothing || work === nothing) && return nothing
    (work isa Tuple && length(work) == 3) || return nothing
    for c = 1:3
        Base.mightalias(candidate, work[c]) &&
            throw(ArgumentError("$context: $candidate_name must not alias work[$c]"))
    end
    return nothing
end

function _collision_selection_work(ps::ParticleSet, selection_work, context::AbstractString)
    N = nparticles(ps)
    if selection_work === nothing
        return falses(N)
    end
    selection_work isa AbstractVector{Bool} ||
        throw(ArgumentError("$context: selection_work must be an AbstractVector{Bool}"))
    axes(selection_work, 1) == Base.OneTo(N) || throw(
        DimensionMismatch(
            "$context: selection_work axes $(axes(selection_work, 1)) do not match particle axes $(Base.OneTo(N))",
        ),
    )
    _mightalias_particle_storage(selection_work, ps) &&
        throw(ArgumentError("$context: selection_work must not alias particle storage"))
    fill!(selection_work, false)
    return selection_work
end

function _collision_velocity_work(ps::ParticleSet{D,T}, work, context::AbstractString) where {D,T}
    if work === nothing
        return ntuple(c -> copy(ps.v[c]), 3)
    end
    (work isa Tuple && length(work) == 3) ||
        throw(ArgumentError("$context: work must be a 3-tuple of velocity vectors"))
    for c = 1:3
        work[c] isa AbstractVector{T} ||
            throw(ArgumentError("$context: work[$c] must be an AbstractVector{$T}"))
        axes(work[c]) == axes(ps.v[c]) || throw(
            DimensionMismatch(
                "$context: work[$c] axes $(axes(work[c])) do not match particle axes $(axes(ps.v[c]))",
            ),
        )
    end
    for c = 1:3
        _mightalias_particle_storage(work[c], ps) &&
            throw(ArgumentError("$context: work must not alias particle storage"))
        for d = (c+1):3
            Base.mightalias(work[c], work[d]) &&
                throw(ArgumentError("$context: work components must not alias"))
        end
    end
    for c = 1:3
        copyto!(work[c], ps.v[c])
    end
    return work
end

"""
    collide_neutral_mcc!(ps::ParticleSet{D,T}, dt; nσ, T_n, m_n=1.0,
                         u_n=(0.0,0.0,0.0), rng=Random.default_rng(),
                         work=nothing) -> ps

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
ionization are the upgrade path). The update is transactional: invalid sampled or
computed states throw without changing `ps`. By default this uses three temporary
velocity vectors; repeated callers can pass a non-aliasing 3-tuple of particle-length
vectors as `work` to reuse that storage. Returns `ps`.
"""
function collide_neutral_mcc!(
    ps::ParticleSet{D,T},
    dt::Real;
    nσ::Real,
    T_n::Real,
    m_n::Real = 1.0,
    u_n::NTuple{3,<:Real} = (0.0, 0.0, 0.0),
    rng = Random.default_rng(),
    work = nothing,
) where {D,T}
    nσT = _require_finite_nonnegative_real("nσ (density×cross-section)", nσ, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "collide_neutral_mcc!")
    TnT = _require_finite_nonnegative_real("T_n", T_n, T)
    mnT = _require_finite_positive_real("m_n", m_n, T)
    mp = _require_finite_positive_real("particle mass ps.m", ps.m, T)
    unx = _require_finite_real("u_n[1]", u_n[1], T)
    uny = _require_finite_real("u_n[2]", u_n[2], T)
    unz = _require_finite_real("u_n[3]", u_n[3], T)
    _require_finite_particle_velocities(ps, "collide_neutral_mcc!")
    N = nparticles(ps)
    (N == 0 || iszero(nσT) || iszero(dtT)) && return ps

    vx, vy, vz = ps.v
    vthn = sqrt(TnT / mnT)                          # neutral thermal speed √(T_n/m_n)
    isfinite(vthn) ||
        throw(ArgumentError("collide_neutral_mcc!: neutral thermal speed must be finite"))
    # Form the mass fractions from a ratio ≤ 1. Directly evaluating mp + mnT or
    # mp*v can overflow even when the centre-of-mass velocity is representable.
    if mp >= mnT
        ratio = mnT / mp
        denom = one(T) + ratio
        μp = one(T) / denom
        μn = ratio / denom
    else
        ratio = mp / mnT
        denom = one(T) + ratio
        μp = ratio / denom
        μn = one(T) / denom
    end
    vnew = _collision_velocity_work(ps, work, "collide_neutral_mcc!")
    nvx, nvy, nvz = vnew
    twoπ = 2 * T(π)
    @inbounds for p = 1:N
        vnx = unx + vthn * randn(rng, T)            # sample a neutral partner
        vny = uny + vthn * randn(rng, T)
        vnz = unz + vthn * randn(rng, T)
        (isfinite(vnx) && isfinite(vny) && isfinite(vnz)) ||
            throw(ArgumentError("collide_neutral_mcc!: sampled neutral velocity must be finite"))
        gx = vx[p] - vnx
        gy = vy[p] - vny
        gz = vz[p] - vnz
        (isfinite(gx) && isfinite(gy) && isfinite(gz)) ||
            throw(ArgumentError("collide_neutral_mcc!: relative velocity exceeds numeric range"))
        gmag = hypot(gx, gy, gz)
        isfinite(gmag) ||
            throw(ArgumentError("collide_neutral_mcc!: relative speed exceeds numeric range"))
        gmag > 0 || continue
        Pcoll = -expm1(-nσT * gmag * dtT)           # 1 − exp(−nσ|g|dt)
        rand(rng, T) < Pcoll || continue
        Vx = μp * vx[p] + μn * vnx                  # centre-of-mass velocity
        Vy = μp * vy[p] + μn * vny
        Vz = μp * vz[p] + μn * vnz
        cosχ = 2 * rand(rng, T) - one(T)            # isotropic elastic scatter in CM
        sinχ = sqrt(max(zero(T), one(T) - cosχ * cosχ))
        φ = twoπ * rand(rng, T)
        outx = Vx + μn * gmag * sinχ * cos(φ)
        outy = Vy + μn * gmag * sinχ * sin(φ)
        outz = Vz + μn * gmag * cosχ
        (isfinite(outx) && isfinite(outy) && isfinite(outz)) ||
            throw(ArgumentError("collide_neutral_mcc!: scattered velocity must be finite"))
        nvx[p] = outx
        nvy[p] = outy
        nvz[p] = outz
    end
    for c = 1:3
        copyto!(ps.v[c], vnew[c])
    end
    return ps
end

# ---------------------------------------------------------------- electron-impact ionization

function _reserve_ionization_ids(
    ids::AbstractVector{UInt64},
    nextid::Union{Nothing,Base.RefValue{UInt64}},
    count::Int,
    species::AbstractString,
)
    count > 0 || throw(ArgumentError("ionization id reservation requires a positive count"))
    livemax = isempty(ids) ? zero(UInt64) : maximum(ids)
    livemax < typemax(UInt64) ||
        throw(ArgumentError("ionize_mcc!: $species particle id space is exhausted"))
    livefirst = livemax + one(UInt64)
    requested = nextid === nothing ? livefirst : max(max(nextid[], one(UInt64)), livefirst)
    countU = UInt64(count)
    countU <= typemax(UInt64) - requested ||
        throw(ArgumentError("ionize_mcc!: $species particle id counter is exhausted"))
    return requested, requested + countU
end

"""
    ionize_mcc!(electrons::ParticleSet{D,T}, ions::ParticleSet{D,T}, dt;
                nσ_iz, E_iz, T_n=0.0, m_n=1.0, u_n=(0.0,0.0,0.0),
                rng=Random.default_rng()) -> nionized

One **electron-impact ionization** MCC substep: `e_fast + N → e_primary + e_secondary + N⁺`.
Each electron with kinetic energy `KE = ½ mₑ|v|² > E_iz` ionizes a background neutral with
probability `P = 1 − exp(−nσ_iz·|v|·dt)`. On ionization:

  * the incident (primary) electron is **cooled by exactly `E_iz`** — its speed is rescaled
    by `√((KE−E_iz)/KE)` (so its new `KE = KE−E_iz`), the energy going into unbinding;
  * a **secondary electron** and a **positive ion** are created at the incident position, each
    with a velocity drawn from the neutral Maxwellian (temperature `T_n`, mass `m_n`, drift
    `u_n`), inheriting the incident macro-particle weight.

`electrons` and `ions` grow by the number of ionizations (in place, via `append_particles!`),
which is returned. Net charge change per event is zero (`+1 e⁻`, `+1 ion` from a neutral). The
neutral reservoir supplies the secondary/ion energy and absorbs recoil momentum; with `T_n=0`
the secondaries are born at rest, so the weighted electron-population energy loses **exactly**
`E_iz Σ_{p∈ionized} w_p` (equal to `nionized·E_iz` only for unit macro-particle weights).

`nσ_iz ≥ 0`, `E_iz ≥ 0`, `T_n ≥ 0`, `dt ≥ 0`, `m_n > 0`. The electron species must
have finite negative charge and the ion species the exactly opposite positive charge, as required
by the neutral pair-creation model. Elastic-secondary reservoir model (a
differential-cross-section secondary spectrum is the upgrade). Returns the ionization count.

**Particle ids.** Pass a persistent per-species monotonic counter `e_nextid`/`i_nextid`
(`Ref{UInt64}`, as [`inject_face_1d!`](@ref) does) to give newborns **globally unique** ids
that are never reused — required when id-keyed provenance/restart is active and particles are
removed between calls. Without a counter the ids are only unique against the current **live**
set (`max(existing)+k`); that is safe for the immediate physics but a removed particle's id can
later be reissued.
"""
function ionize_mcc!(
    electrons::ParticleSet{D,T},
    ions::ParticleSet{D,T},
    dt::Real;
    nσ_iz::Real,
    E_iz::Real,
    T_n::Real = 0.0,
    m_n::Real = 1.0,
    u_n::NTuple{3,<:Real} = (0.0, 0.0, 0.0),
    e_nextid::Union{Nothing,Base.RefValue{UInt64}} = nothing,
    i_nextid::Union{Nothing,Base.RefValue{UInt64}} = nothing,
    rng = Random.default_rng(),
) where {D,T}
    nσT = _require_finite_nonnegative_real("nσ_iz", nσ_iz, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "ionize_mcc!")
    Eiz = _require_finite_nonnegative_real("E_iz", E_iz, T)
    TnT = _require_finite_nonnegative_real("T_n", T_n, T)
    mnT = _require_finite_positive_real("m_n", m_n, T)
    me = _require_finite_positive_real("electron mass", electrons.m, T)
    _require_finite_positive_real("ion mass", ions.m, T)
    qe = _require_finite_real("electron charge", electrons.q, T)
    qi = _require_finite_real("ion charge", ions.q, T)
    qe < zero(T) || throw(ArgumentError("ionize_mcc!: electron species charge must be negative"))
    qi > zero(T) || throw(ArgumentError("ionize_mcc!: ion species charge must be positive"))
    qi == -qe || throw(
        ArgumentError("ionize_mcc!: electron and ion species charges must be exactly opposite"),
    )
    unx = _require_finite_real("u_n[1]", u_n[1], T)
    uny = _require_finite_real("u_n[2]", u_n[2], T)
    unz = _require_finite_real("u_n[3]", u_n[3], T)
    _require_finite_particle_velocities(electrons, "ionize_mcc!")
    Ne = nparticles(electrons)
    (Ne == 0 || iszero(nσT) || iszero(dtT)) && return 0

    vthn = sqrt(TnT / mnT)
    isfinite(vthn) || throw(ArgumentError("ionize_mcc!: neutral thermal speed must be finite"))
    evx, evy, evz = electrons.v
    ex = electrons.x
    ew = electrons.weight
    ionization_speed = sqrt(T(2)) * sqrt(Eiz) / sqrt(me)

    born = Int[]                                   # incident electrons that ionized
    @inbounds for p = 1:Ne
        speed = hypot(evx[p], evy[p], evz[p])
        speed > ionization_speed || continue       # below the ionization threshold
        Pcoll = -expm1(-nσT * speed * dtT)
        rand(rng, T) < Pcoll || continue
        push!(born, p)
    end
    nb = length(born)
    nb == 0 && return 0

    # These fields are copied into newborns. Validate the complete selected batch
    # before reserving ids, cooling a primary, or appending either species.
    @inbounds for p in born
        wp = ew[p]
        (isfinite(wp) && wp >= zero(T)) || throw(
            ArgumentError(
                "ionize_mcc!: ionized electron $p weight must be finite and non-negative",
            ),
        )
        for d = 1:D
            isfinite(ex[d][p]) || throw(
                ArgumentError(
                    "ionize_mcc!: ionized electron $p position component $d must be finite",
                ),
            )
        end
    end

    # Reserve both species' complete id ranges before cooling a primary or appending a
    # newborn. This keeps id exhaustion an exception-safe input error instead of leaving
    # cooled primaries, duplicate ids, or a wrapped monotonic counter.
    first_e, next_e = _reserve_ionization_ids(electrons.id, e_nextid, nb, "electron")
    first_i, next_i = _reserve_ionization_ids(ions.id, i_nextid, nb, "ion")

    # Build newborn secondary electrons + ions off-state (batched: one append each).
    new_e = ParticleSet{D,T}(nb; q = electrons.q, m = electrons.m)
    new_i = ParticleSet{D,T}(nb; q = ions.q, m = ions.m)
    @inbounds for (k, p) in enumerate(born)
        for d = 1:D
            new_e.x[d][k] = ex[d][p]               # born at the incident position
            new_i.x[d][k] = ex[d][p]
        end
        ev1 = unx + vthn * randn(rng, T)           # secondary e⁻ from neutral Maxwellian
        ev2 = uny + vthn * randn(rng, T)
        ev3 = unz + vthn * randn(rng, T)
        iv1 = unx + vthn * randn(rng, T)           # ion from the neutral Maxwellian
        iv2 = uny + vthn * randn(rng, T)
        iv3 = unz + vthn * randn(rng, T)
        (
            isfinite(ev1) &&
            isfinite(ev2) &&
            isfinite(ev3) &&
            isfinite(iv1) &&
            isfinite(iv2) &&
            isfinite(iv3)
        ) || throw(ArgumentError("ionize_mcc!: sampled newborn velocities must be finite"))
        new_e.v[1][k] = ev1
        new_e.v[2][k] = ev2
        new_e.v[3][k] = ev3
        new_i.v[1][k] = iv1
        new_i.v[2][k] = iv2
        new_i.v[3][k] = iv3
        new_e.weight[k] = ew[p]                    # inherit the incident macro-particle weight
        new_i.weight[k] = ew[p]
        new_e.id[k] = first_e + UInt64(k - 1)
        new_i.id[k] = first_i + UInt64(k - 1)
    end
    @inbounds for p in born
        speed = hypot(evx[p], evy[p], evz[p])
        loss_ratio = ionization_speed / speed
        # KE'/KE = 1 - Eiz/KE = 1 - (v_threshold/|v|)^2.
        # This ratio neither squares a tiny speed nor forms Inf/Inf at high energy.
        scale = sqrt(max(zero(T), one(T) - loss_ratio * loss_ratio))
        evx[p] *= scale
        evy[p] *= scale
        evz[p] *= scale
    end
    append_particles!(electrons, new_e)
    append_particles!(ions, new_i)
    e_nextid === nothing || (e_nextid[] = next_e)
    i_nextid === nothing || (i_nextid[] = next_i)
    return nb
end
