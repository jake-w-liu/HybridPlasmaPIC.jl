# raycon_conversion.jl — mode-conversion analysis (port of the dispertok.m
# det-U path + the ray.m 'convert_list' algorithm + cgamma.m).
#
# The conversion machinery works with the DETERMINANT dispersion function
# U = D11·D22 − |D12|² of the 2×2 cold tensor (dispertok convention), expanded
# to second order around the detection point in the osculating plane spanned by
# the ray velocity ż and acceleration z̈ (Tracy, Kaufman & Jaun 2001/2007):
#
#   η² = |e₁ᴴ·DD·e₂|²/⟨∇α, J∇λ⟩   (coupling),  τ = e^{−πη²}  (transmission),
#   β = √(2πτ) / (η·Γ(−iη²))      (conversion phase/amplitude factor).
#
# Upstream parity notes (see the port-notes doc):
#  - the (r,z) chain here deliberately has NO field-curvature corrections and
#    the second-derivative chain omits ∂²(s,t)/∂(r,z)² terms — exactly as
#    dispertok.m; the exact-eigenvalue tracing lives in raycon_dispersion.jl;
#  - cld3x3 conversion is refused: upstream's 3×3 second derivatives carry a
#    known sign bug (dispertok.m "FIX THIS: sgn_fix NOT ADDED BELOW");
#  - a malformed hyperbola (|1+λ₄/λ₁| ≥ 0.1 or |λ₄/λ₃| ≤ 4) yields τ = 0 and
#    no ray split, as upstream.

"""
    cgamma(z) -> Complex

Complex gamma function (port of `cgamma.m`: Stirling series after shifting
`Re(z)` into [9, 10), conjugate symmetry for `Im(z) < 0`). Accurate to ~12
digits on the upstream test domain −10 < Re(z) < 10, |Im(z)| < 10.
"""
function cgamma(z::Number)
    zc = ComplexF64(z)
    flip = imag(zc) < 0
    flip && (zc = conj(zc))
    lng = zero(ComplexF64)
    nshift = floor(Int, 9 - real(zc))
    for _ = 1:max(nshift, 0)
        lng -= log(zc)
        zc += 1
    end
    zm1 = 1 / zc
    zm2 = zm1 * zm1
    zm3 = zm1 * zm2
    zm5 = zm3 * zm2
    zm7 = zm5 * zm2
    zm9 = zm7 * zm2
    zm11 = zm9 * zm2
    zm13 = zm11 * zm2
    zm15 = zm13 * zm2
    lng -= zm15 * (3617 / 122400)
    lng += zm13 / 156
    lng -= zm11 * (691 / 360360)
    lng += zm9 / 1188
    lng -= zm7 / 1680
    lng += zm5 / 1260
    lng -= zm3 / 360
    lng += zm1 / 12 + (zc - 0.5) * log(zc) - zc + 0.5 * log(2π)
    flip && (lng = conj(lng))
    return exp(lng)
end

# det-U (cld2x2) and its derivatives in dispertok's variables. Returns the
# tensor, U, V, first derivatives in (r,z,kr,kz,ω), the gDij phase-space
# gradients, and (need2nd) all second derivatives needed for the T2 Hessian.
function _detU_core(prob::RayconProblem, y::NTuple{4,Float64}; need2nd::Bool = false)
    prob.model === :cld2x2 || throw(
        ArgumentError(
            "mode-conversion analysis supports :cld2x2 only (upstream's 3×3 second " *
            "derivatives carry a known sign bug and are not ported)",
        ),
    )
    r, zc, kr, kz = y
    rho, theta = _poloidal_coords(prob.eq, r, zc)
    geo = magnetic_geometry(prob.eq, rho, theta)
    st = _stix_local(prob, geo)
    om = prob.omega
    coom = prob.cnst.c / om
    coomsq = coom^2
    g = geo

    kn = kr * g.ener + prob.kphi * g.enef + kz * g.enez
    kb = kr * g.eber + prob.kphi * g.ebef + kz * g.ebez
    kp = kr * g.eper + prob.kphi * g.epef + kz * g.epez
    Nn = coom * kn
    Nb = coom * kb
    Np = coom * kp
    Nn2 = Nn^2
    Nb2 = Nb^2
    Np2 = Np^2

    D11 = Nb2 + Np2 - st.S
    D12 = -Nn * Nb - im * st.D
    D22 = Nn2 + Np2 - st.S
    cD12 = conj(D12)
    U = real(D11 * D22 - D12 * cD12)
    V = 0.5 * (D11 + D22)

    # element derivatives in (s, t, kn, kb, kp, ω) — dispertok lines 141-190
    dD11ds = -st.dSds
    dD11dt = -st.dSdt
    dD12ds = -im * st.dDds
    dD12dt = -im * st.dDdt
    dD22ds = -st.dSds
    dD22dt = -st.dSdt
    dD11dkn = 0.0
    dD11dkb = 2 * Nb * coom
    dD11dkp = 2 * Np * coom
    dD12dkn = -Nb * coom
    dD12dkb = -Nn * coom
    dD12dkp = 0.0
    dD22dkn = 2 * Nn * coom
    dD22dkb = 0.0
    dD22dkp = 2 * Np * coom
    dD11dom = -2 / om * (Nb2 + Np2) - st.dSdom
    # upstream sign bug fixed as in raycon_dispersion.jl (∂(−NnNb)/∂ω = +2NnNb/ω)
    dD12dom = +2 / om * (Nn * Nb) - im * st.dDdom
    dD22dom = -2 / om * (Nn2 + Np2) - st.dSdom

    # U derivatives in dispertok variables
    _dU(d11, d22, d12) = real(d11 * D22 + D11 * d22 - 2 * real(D12 * conj(d12)))
    dUds = _dU(dD11ds, dD22ds, dD12ds)
    dUdt = _dU(dD11dt, dD22dt, dD12dt)
    dUdkn = _dU(dD11dkn, dD22dkn, dD12dkn)
    dUdkb = _dU(dD11dkb, dD22dkb, dD12dkb)
    dUdkp = _dU(dD11dkp, dD22dkp, dD12dkp)
    dUdom = _dU(dD11dom, dD22dom, dD12dom)

    # chain to (r, z, kr, kz) — NO curvature corrections (dispertok parity)
    dUdr = dUds * g.dsdr + dUdt * g.dtdr
    dUdz = dUds * g.dsdz + dUdt * g.dtdz
    dUdkr = dUdkn * g.ener + dUdkb * g.eber + dUdkp * g.eper
    dUdkz = dUdkn * g.enez + dUdkb * g.ebez + dUdkp * g.epez

    # gDij phase-space gradients (dispertok lines 230-249)
    dD11dr = dD11ds * g.dsdr + dD11dt * g.dtdr
    dD11dz = dD11ds * g.dsdz + dD11dt * g.dtdz
    dD12dr = dD12ds * g.dsdr + dD12dt * g.dtdr
    dD12dz = dD12ds * g.dsdz + dD12dt * g.dtdz
    dD22dr = dD22ds * g.dsdr + dD22dt * g.dtdr
    dD22dz = dD22ds * g.dsdz + dD22dt * g.dtdz
    dD11dkr = dD11dkn * g.ener + dD11dkb * g.eber + dD11dkp * g.eper
    dD11dkz = dD11dkn * g.enez + dD11dkb * g.ebez + dD11dkp * g.epez
    dD12dkr = dD12dkn * g.ener + dD12dkb * g.eber + dD12dkp * g.eper
    dD12dkz = dD12dkn * g.enez + dD12dkb * g.ebez + dD12dkp * g.epez
    dD22dkr = dD22dkn * g.ener + dD22dkb * g.eber + dD22dkp * g.eper
    dD22dkz = dD22dkn * g.enez + dD22dkb * g.ebez + dD22dkp * g.epez
    gD11 = ComplexF64[dD11dr, dD11dz, dD11dkr, dD11dkz]
    gD12 = ComplexF64[dD12dr, dD12dz, dD12dkr, dD12dkz]
    gD22 = ComplexF64[dD22dr, dD22dz, dD22dkr, dD22dkz]

    base = (;
        rho,
        theta,
        geo,
        st,
        DD = ComplexF64[D11 D12; cD12 D22],
        U,
        V,
        dUdr,
        dUdz,
        dUdkr,
        dUdkz,
        dUdom,
        gD11,
        gD12,
        gD22,
    )
    need2nd || return base

    # second derivatives of the elements (dispertok lines 204-212)
    dD11ds2 = -st.dSds2
    dD11dst = -st.dSdst
    dD11dt2 = -st.dSdt2
    dD12ds2 = -im * st.dDds2
    dD12dst = -im * st.dDdst
    dD12dt2 = -im * st.dDdt2
    dD22ds2 = -st.dSds2
    dD22dst = -st.dSdst
    dD22dt2 = -st.dSdt2

    # U second derivatives in dispertok variables (lines 377-391)
    _dU2(d11a, d22a, d12a, d11b, d22b, d12b, d11ab, d22ab, d12ab) = real(
        d11ab * D22 + d11a * d22b + d11b * d22a + D11 * d22ab -
        2 * real(d12b * conj(d12a) + D12 * conj(d12ab)),
    )
    dUds2 = _dU2(dD11ds, dD22ds, dD12ds, dD11ds, dD22ds, dD12ds, dD11ds2, dD22ds2, dD12ds2)
    dUdst = _dU2(dD11ds, dD22ds, dD12ds, dD11dt, dD22dt, dD12dt, dD11dst, dD22dst, dD12dst)
    dUdt2 = _dU2(dD11dt, dD22dt, dD12dt, dD11dt, dD22dt, dD12dt, dD11dt2, dD22dt2, dD12dt2)
    dUdkn2 = _dU2(dD11dkn, dD22dkn, dD12dkn, dD11dkn, dD22dkn, dD12dkn, 0.0, 2 * coomsq, 0.0)
    dUdknkb = _dU2(dD11dkn, dD22dkn, dD12dkn, dD11dkb, dD22dkb, dD12dkb, 0.0, 0.0, -coomsq)
    dUdknkp = _dU2(dD11dkn, dD22dkn, dD12dkn, dD11dkp, dD22dkp, dD12dkp, 0.0, 0.0, 0.0)
    dUdkb2 = _dU2(dD11dkb, dD22dkb, dD12dkb, dD11dkb, dD22dkb, dD12dkb, 2 * coomsq, 0.0, 0.0)
    dUdkbkp = _dU2(dD11dkb, dD22dkb, dD12dkb, dD11dkp, dD22dkp, dD12dkp, 0.0, 0.0, 0.0)
    dUdkp2 = _dU2(dD11dkp, dD22dkp, dD12dkp, dD11dkp, dD22dkp, dD12dkp, 2 * coomsq, 2 * coomsq, 0.0)
    _dUmix(d11a, d22a, d12a, d11b, d22b, d12b) =
        real(d11a * d22b + d11b * d22a - 2 * real(d12b * conj(d12a)))
    dUdskn = _dUmix(dD11ds, dD22ds, dD12ds, dD11dkn, dD22dkn, dD12dkn)
    dUdskb = _dUmix(dD11ds, dD22ds, dD12ds, dD11dkb, dD22dkb, dD12dkb)
    dUdskp = _dUmix(dD11ds, dD22ds, dD12ds, dD11dkp, dD22dkp, dD12dkp)
    dUdtkn = _dUmix(dD11dt, dD22dt, dD12dt, dD11dkn, dD22dkn, dD12dkn)
    dUdtkb = _dUmix(dD11dt, dD22dt, dD12dt, dD11dkb, dD22dkb, dD12dkb)
    dUdtkp = _dUmix(dD11dt, dD22dt, dD12dt, dD11dkp, dD22dkp, dD12dkp)

    # chain to (r, z, kr, kz) — upstream 2nd-order chain (lines 755-787),
    # which omits ∂²(s,t)/∂(r,z)² terms (parity-preserved)
    dUdrs = dUds2 * g.dsdr + dUdst * g.dtdr
    dUdzs = dUds2 * g.dsdz + dUdst * g.dtdz
    dUdrt = dUdst * g.dsdr + dUdt2 * g.dtdr
    dUdzt = dUdst * g.dsdz + dUdt2 * g.dtdz
    dUdskr = dUdskn * g.ener + dUdskb * g.eber + dUdskp * g.eper
    dUdskz = dUdskn * g.enez + dUdskb * g.ebez + dUdskp * g.epez
    dUdtkr = dUdtkn * g.ener + dUdtkb * g.eber + dUdtkp * g.eper
    dUdtkz = dUdtkn * g.enez + dUdtkb * g.ebez + dUdtkp * g.epez
    dUdkrkn = dUdkn2 * g.ener + dUdknkb * g.eber + dUdknkp * g.eper
    dUdkrkb = dUdknkb * g.ener + dUdkb2 * g.eber + dUdkbkp * g.eper
    dUdkrkp = dUdknkp * g.ener + dUdkbkp * g.eber + dUdkp2 * g.eper
    dUdkzkn = dUdkn2 * g.enez + dUdknkb * g.ebez + dUdknkp * g.epez
    dUdkzkb = dUdknkb * g.enez + dUdkb2 * g.ebez + dUdkbkp * g.epez
    dUdkzkp = dUdknkp * g.enez + dUdkbkp * g.ebez + dUdkp2 * g.epez
    dUdr2 = dUdrs * g.dsdr + dUdrt * g.dtdr
    dUdz2 = dUdzs * g.dsdz + dUdzt * g.dtdz
    dUdrz = dUdrs * g.dsdz + dUdrt * g.dtdz
    dUdrkr = dUdskr * g.dsdr + dUdtkr * g.dtdr
    dUdrkz = dUdskz * g.dsdr + dUdtkz * g.dtdr
    dUdzkr = dUdskr * g.dsdz + dUdtkr * g.dtdz
    dUdzkz = dUdskz * g.dsdz + dUdtkz * g.dtdz
    dUdkr2 = dUdkrkn * g.ener + dUdkrkb * g.eber + dUdkrkp * g.eper
    dUdkz2 = dUdkzkn * g.enez + dUdkzkb * g.ebez + dUdkzkp * g.epez
    dUdkrkz = dUdkrkn * g.enez + dUdkrkb * g.ebez + dUdkrkp * g.epez

    T2 = [
        dUdr2 dUdrz dUdrkr dUdrkz
        dUdrz dUdz2 dUdzkr dUdzkz
        dUdrkr dUdzkr dUdkr2 dUdkrkz
        dUdrkz dUdzkz dUdkrkz dUdkz2
    ]
    return (; base..., T2)
end

# Osculating-plane quantities at a point: unit velocity/acceleration directions
# (normalized by the symplectic area A), H Hessian, saddle displacement.
function _osculating(core, eqv::Vector{Float64}, epv::Vector{Float64})
    dD11dq = sum(eqv .* core.gD11)
    dD11dp = sum(epv .* core.gD11)
    dD12dq = sum(eqv .* core.gD12)
    dD12dp = sum(epv .* core.gD12)
    dD22dq = sum(eqv .* core.gD22)
    dD22dp = sum(epv .* core.gD22)
    D11 = core.DD[1, 1]
    D12 = core.DD[1, 2]
    D22 = core.DD[2, 2]
    H11 = real(2 * dD11dq * dD22dq - 2 * dD12dq * conj(dD12dq))
    H22 = real(2 * dD11dp * dD22dp - 2 * dD12dp * conj(dD12dp))
    H12 = real(dD11dq * dD22dp + dD11dp * dD22dq - 2 * real(dD12dq * conj(dD12dp)))
    dUdq = real(dD11dq * D22 + D11 * dD22dq - 2 * real(D12 * conj(dD12dq)))
    dUdp = real(dD11dp * D22 + D11 * dD22dp - 2 * real(D12 * conj(dD12dp)))
    detH = H11 * H22 - H12^2
    Hm11 = H22 / detH
    Hm12 = -H12 / detH
    Hm22 = H11 / detH
    eta2 = 0.5 / sqrt(abs(detH)) * abs(Hm11 * dUdq^2 + 2 * Hm12 * (dUdq * dUdp) + Hm22 * dUdp^2)
    qst = Hm11 * dUdq + Hm12 * dUdp
    pst = Hm12 * dUdq + Hm22 * dUdp
    zinzst = -(qst .* eqv .+ pst .* epv)          # displacement toward the saddle
    return (; eta2, zinzst)
end

# dominant-direction solve for a 2x2 complex matrix: upstream's power
# iteration (tol 1e-4, ≤100 iterations, start [1;1]) first, for parity. Its
# per-component relative-change metric never converges for some spectra (a
# real decaying component has constant relative change 1 − λ₂/λ₁), so fall
# back to the exact dominant eigenvector — verified to reproduce the iterated
# result to ~1e-6 where both work.
function _power_iterate(gd::Matrix{ComplexF64})
    new = ComplexF64[1, 1]
    old = ComplexF64[0, 0]
    for it = 1:100
        n1 = maximum(abs.((new .- old) ./ new))
        n2 = maximum(abs.((new .+ old) ./ new))
        (n1 <= 1e-4 || n2 <= 1e-4) && return new
        old = new
        new = gd * new
        new = new / norm(new)
    end
    F = eigen(gd)
    v = F.vectors[:, argmax(abs.(F.values))]
    return v / norm(v)
end

"""
    RayconConversion

Result of a mode-conversion analysis at a detection point: the `saddle` point
`z*`, the `converted` ray start (the incoming point, continuing on the followed
eigenvalue branch), the `transmitted` ray start (on the other branch, across
the avoided crossing), the coupling `eta2`, transmission `tau = e^{−πη²}`,
conversion factor `beta`, the a-priori estimate `eta2_estimate`, and validity
flags (`converged` saddle iteration, `hyperbola_ok` normal form).
"""
struct RayconConversion
    saddle::NTuple{4,Float64}
    incoming::NTuple{4,Float64}
    converted::NTuple{4,Float64}
    transmitted::NTuple{4,Float64}
    eta2::Float64
    tau::Float64
    beta::ComplexF64
    eta2_estimate::Float64
    converged::Bool
    hyperbola_ok::Bool
    transmitted_ok::Bool
end

"whether the analysis produced a usable ray split (all three validity gates)"
is_valid(c::RayconConversion) = c.converged && c.hyperbola_ok && c.transmitted_ok

"""
    analyze_conversion(prob, z0, zdot, zddot) -> RayconConversion

Full mode-conversion analysis at detection point `z0` with ray velocity `zdot`
and acceleration `zddot` (port of `ray.m 'convert_list'` + the dispertok
'Mon'/'Sdl'/'Trs'/'Cnv' operators): saddle-point Newton iteration (≤30 steps,
tol 1e-4), symplectic normal form `eig(J₄·∇∇U)` with the upstream hyperbola
quality gates, uncoupled polarizations by power iteration, coupling constant
`η`, `τ = e^{−πη²}`, `β = √(2πτ)/(η·Γ(−iη²))`, and the transmitted-ray launch
point found by a root solve along the incoming→saddle direction.
"""
function analyze_conversion(
    prob::RayconProblem,
    z0::AbstractVector{<:Real},
    zdot::AbstractVector{<:Real},
    zddot::AbstractVector{<:Real},
)
    (length(z0) == 4 && length(zdot) == 4 && length(zddot) == 4) ||
        throw(ArgumentError("z0, zdot, zddot must be 4-vectors (r, z, kr, kz)"))
    all(isfinite, z0) && all(isfinite, zdot) && all(isfinite, zddot) ||
        throw(ArgumentError("conversion inputs must be finite"))
    y0 = (Float64(z0[1]), Float64(z0[2]), Float64(z0[3]), Float64(z0[4]))

    # osculating directions, FIXED during the saddle iteration (upstream)
    A = zdot[1] * zddot[3] + zdot[2] * zddot[4] - zdot[3] * zddot[1] - zdot[4] * zddot[2]
    A != 0 || throw(ArgumentError("degenerate osculating plane (velocity ∥ acceleration)"))
    sqmA = 1 / sqrt(abs(A))
    eqv = sqmA .* Float64.(zdot[1:4])
    epv = sqmA .* Float64.(zddot[1:4])

    # the saddle search can wander outside the plasma (s ≥ 1 profile domain);
    # treat that like a non-converged saddle (upstream aborts the conversion)
    local eta2_estimate::Float64, saddle::NTuple{4,Float64}
    converged = false
    zst = collect(y0)
    try
        core0 = _detU_core(prob, y0; need2nd = true)
        osc0 = _osculating(core0, eqv, epv)
        eta2_estimate = osc0.eta2
        zst = collect(y0) .+ osc0.zinzst
        for _ = 1:30
            c = _detU_core(prob, (zst[1], zst[2], zst[3], zst[4]); need2nd = true)
            o = _osculating(c, eqv, epv)
            znew = zst .+ o.zinzst
            if maximum(abs.((znew .- zst) ./ znew)) <= 1e-4
                zst = znew
                converged = true
                break
            end
            zst = znew
        end
    catch e
        e isa DomainError || rethrow()
        eta2_estimate = NaN
        converged = false
    end
    saddle = (zst[1], zst[2], zst[3], zst[4])
    if !converged
        return RayconConversion(
            saddle,
            y0,
            y0,
            y0,
            NaN,
            NaN,
            complex(NaN, NaN),
            eta2_estimate,
            false,
            false,
            false,
        )
    end

    coeff = _conversion_coefficients(prob, saddle, y0)
    if !coeff.hyperbola_ok
        # upstream: tau = 0, no split
        return RayconConversion(
            saddle,
            y0,
            y0,
            y0,
            NaN,
            0.0,
            complex(0.0, 0.0),
            eta2_estimate,
            true,
            false,
            false,
        )
    end
    return RayconConversion(
        saddle,
        y0,
        y0,
        coeff.transmitted,
        coeff.eta2,
        coeff.tau,
        coeff.beta,
        eta2_estimate,
        true,
        true,
        coeff.transmitted_ok,
    )
end

# Normal form, coupling coefficients and transmitted-ray relaunch at a given
# saddle point (the dispertok 'Mch'/'Cnv'/'Trs' block). Split out so it can be
# exercised directly against MATLAB reference data.
function _conversion_coefficients(
    prob::RayconProblem,
    saddle::NTuple{4,Float64},
    y0::NTuple{4,Float64},
)
    cst = _detU_core(prob, saddle; need2nd = true)
    J4 = [0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; -1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0]
    F = eigen(J4 * cst.T2)
    ind = sortperm(real.(F.values))
    vals = real.(F.values[ind])
    opposite = vals[4] / vals[1]
    separate = abs(vals[4] / vals[3])
    hyperbola_ok = abs(1 + opposite) < 0.1 && separate > 4
    hyperbola_ok || return (;
        hyperbola_ok,
        eta2 = NaN,
        tau = 0.0,
        beta = complex(0.0, 0.0),
        transmitted = y0,
        transmitted_ok = false,
    )
    vp = F.vectors[:, ind[4]]
    vm = F.vectors[:, ind[1]]

    # uncoupled polarizations from the directional tensor gradients
    # (MATLAB reshape([gD11·v, gD12·v, cgD12·v, gD22·v], 2, 2) is column-major)
    _gdmat(v) = ComplexF64[
        sum(cst.gD11 .* v) sum(conj.(cst.gD12) .* v)
        sum(cst.gD12 .* v) sum(cst.gD22 .* v)
    ]
    pol1 = _power_iterate(_gdmat(vp))
    pol2 = _power_iterate(_gdmat(vm))

    # ∇(eᴴ·DD·e) for each uncoupled polarization
    _gdpol(e) =
        real.(
            cst.gD11 .* abs2(e[1]) .+ 2 .* real.(cst.gD12 .* (e[1] * conj(e[2]))) .+
            cst.gD22 .* abs2(e[2]),
        )
    gdalf = _gdpol(pol1)
    gdlam = _gdpol(pol2)
    braket = transpose(gdalf) * J4 * gdlam
    # upstream evaluates pol1'·reshape(DD,2,2)·pol2 with MATLAB's column-major
    # reshape, i.e. the TRANSPOSED tensor — consistent with the transposed gd
    # matrices the polarizations are iterated on (conjugate-polarization
    # convention throughout). Verified against the original code at the saddle.
    eta = (pol1' * transpose(cst.DD) * pol2) / sqrt(complex(braket))
    eta2 = real(eta * conj(eta))
    tau = exp(-π * eta2)
    beta = sqrt(2π * tau) / (eta * cgamma(-im * eta2))

    # transmitted ray: root of U along the incoming→saddle direction (≈ the
    # mirror point through the saddle, upstream guess factor 2)
    dirvec = collect(saddle) .- collect(y0)
    ftrs = fact -> begin
        pt = collect(y0) .+ fact .* dirvec
        core = _detU_core(prob, (pt[1], pt[2], pt[3], pt[4]))
        core.U
    end
    transmitted = y0
    transmitted_ok = true
    fact = 2.0
    try
        fact = _fzero_near(ftrs, 2.0)
    catch
        transmitted_ok = false
    end
    if transmitted_ok
        ptrans = collect(y0) .+ fact .* dirvec
        transmitted = (ptrans[1], ptrans[2], ptrans[3], ptrans[4])
    end
    return (; hyperbola_ok, eta2, tau, beta = ComplexF64(beta), transmitted, transmitted_ok)
end
