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
#  - cld3x3 conversion is the CORRECTED extension (upstream's 3×3 second
#    derivatives carry a known sign bug — dispertok.m "FIX THIS"): derivatives
#    come from exact determinant-sampling identities and the coupling is
#    evaluated in the Tracy–Kaufman 2-D near-null subspace of the tensor;
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

# real determinant of a small Hermitian matrix (for the conversion layer all
# sampled matrices A ± G stay Hermitian, so the determinant is exactly real)
@inline _hdet(A::Matrix{ComplexF64}) =
    size(A, 1) == 2 ? real(A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) :
    real(
        A[1, 1] * (A[2, 2] * A[3, 3] - A[2, 3] * A[3, 2]) -
        A[1, 2] * (A[2, 1] * A[3, 3] - A[2, 3] * A[3, 1]) +
        A[1, 3] * (A[2, 1] * A[3, 2] - A[2, 2] * A[3, 1]),
    )

# adjugate of a Hermitian 3×3 (transpose of the cofactor matrix)
function _adj3(A::Matrix{ComplexF64})
    return ComplexF64[
        (A[2, 2]*A[3, 3]-A[2, 3]*A[3, 2]) -(A[1, 2] * A[3, 3] - A[1, 3] * A[3, 2]) (A[1, 2]*A[2, 3]-A[1, 3]*A[2, 2])
        -(A[2, 1] * A[3, 3] - A[2, 3] * A[3, 1]) (A[1, 1]*A[3, 3]-A[1, 3]*A[3, 1]) -(A[1, 1] * A[2, 3] - A[1, 3] * A[2, 1])
        (A[2, 1]*A[3, 2]-A[2, 2]*A[3, 1]) -(A[1, 1] * A[3, 2] - A[1, 2] * A[3, 1]) (A[1, 1]*A[2, 2]-A[1, 2]*A[2, 1])
    ]
end

# det-U and its derivatives in dispertok's variables (:cld2x2 pinned against
# the original MATLAB; :cld3x3 is the corrected extension — upstream's 3×3
# second derivatives carry a known sign bug, so they are DERIVED here rather
# than transcribed, via exact determinant-sampling identities). Returns the
# tensor, U, first derivatives in (r,z,kr,kz,ω), the per-direction element-
# derivative matrices G[α] = ∂DD/∂z_α, and (need2nd) the T2 Hessian.
function _detU_core(prob::RayconProblem, y::NTuple{4,Float64}; need2nd::Bool = false)
    prob.model === :cld3x3 && return _detU_core3(prob, y; need2nd)
    prob.model === :cld2x2 || throw(
        ArgumentError("mode-conversion analysis requires :cld2x2 or :cld3x3 (got $(prob.model))"),
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
    G = ntuple(α -> ComplexF64[gD11[α] gD12[α]; conj(gD12[α]) gD22[α]], 4)

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
        G,
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

# cld3x3 det-U core: the CORRECTED 3×3 extension of the dispertok conversion
# layer (upstream's own 3×3 second derivatives are sign-broken; here U and its
# derivatives are derived from the real form of the Hermitian determinant
#
#   U = D11·D22·D33 + 2·aR·b·c − D11·c² − D22·b² − D33·(aR² + Dst²),
#
# with D12 = aR − i·Dst, D13 = b, D23 = c all built from real quantities, so
# every derivative is a compact product rule instead of upstream's ~200-line
# expansions). Variable structure and chains follow the validated 2×2 path
# (elements as functions of (s,t,kn,kb,kp,ω), fixed basis, no ∂²(s,t) terms).
function _detU_core3(prob::RayconProblem, y::NTuple{4,Float64}; need2nd::Bool = false)
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
    D22 = Nn2 + Np2 - st.S
    D33 = Nn2 + Nb2 - st.P
    aR = -Nn * Nb                 # Re(D12)
    Dst = st.D                    # −Im(D12)
    b = -Nn * Np                  # D13 (real)
    c = -Nb * Np                  # D23 (real)
    U = D11 * D22 * D33 + 2 * aR * b * c - D11 * c^2 - D22 * b^2 - D33 * (aR^2 + Dst^2)
    DD = ComplexF64[
        D11 (aR-im*Dst) b
        (aR+im*Dst) D22 c
        b c D33
    ]

    # ∂U/∂(element) partials of the real determinant form
    pD11 = D22 * D33 - c^2
    pD22 = D11 * D33 - b^2
    pD33 = D11 * D22 - aR^2 - Dst^2
    paR = 2 * (b * c - D33 * aR)
    pb = 2 * (aR * c - D22 * b)
    pc = 2 * (aR * b - D11 * c)
    pDst = -2 * D33 * Dst

    # element derivatives per variable x ∈ (s, t, kn, kb, kp, ω)
    twoc = 2 * coom
    # dU/dx = pD11·dD11 + pD22·dD22 + pD33·dD33 + paR·daR + pb·db + pc·dc + pDst·dDst
    dUds = pD11 * (-st.dSds) + pD22 * (-st.dSds) + pD33 * (-st.dPds) + pDst * st.dDds
    dUdt = pD11 * (-st.dSdt) + pD22 * (-st.dSdt) + pD33 * (-st.dPdt) + pDst * st.dDdt
    dUdkn = pD22 * (twoc * Nn) + pD33 * (twoc * Nn) + paR * (-coom * Nb) + pb * (-coom * Np)
    dUdkb = pD11 * (twoc * Nb) + pD33 * (twoc * Nb) + paR * (-coom * Nn) + pc * (-coom * Np)
    dUdkp = pD11 * (twoc * Np) + pD22 * (twoc * Np) + pb * (-coom * Nn) + pc * (-coom * Nb)
    dUdom =
        pD11 * (-2 / om * (Nb2 + Np2) - st.dSdom) +
        pD22 * (-2 / om * (Nn2 + Np2) - st.dSdom) +
        pD33 * (-2 / om * (Nn2 + Nb2) - st.dPdom) +
        paR * (-2 * aR / om) +
        pb * (-2 * b / om) +
        pc * (-2 * c / om) +
        pDst * st.dDdom

    # chain to (r, z, kr, kz) — dispertok conventions (no curvature terms)
    dUdr = dUds * g.dsdr + dUdt * g.dtdr
    dUdz = dUds * g.dsdz + dUdt * g.dtdz
    dUdkr = dUdkn * g.ener + dUdkb * g.eber + dUdkp * g.eper
    dUdkz = dUdkn * g.enez + dUdkb * g.ebez + dUdkp * g.epez

    # per-direction element-derivative matrices G[α] = ∂DD/∂z_α
    sα = (g.dsdr, g.dsdz, 0.0, 0.0)
    tα = (g.dtdr, g.dtdz, 0.0, 0.0)
    enα = (0.0, 0.0, g.ener, g.enez)
    ebα = (0.0, 0.0, g.eber, g.ebez)
    epα = (0.0, 0.0, g.eper, g.epez)
    G = ntuple(4) do α
        dD11 = -st.dSds * sα[α] - st.dSdt * tα[α] + twoc * Nb * ebα[α] + twoc * Np * epα[α]
        dD22 = -st.dSds * sα[α] - st.dSdt * tα[α] + twoc * Nn * enα[α] + twoc * Np * epα[α]
        dD33 = -st.dPds * sα[α] - st.dPdt * tα[α] + twoc * Nn * enα[α] + twoc * Nb * ebα[α]
        daR = -coom * (Nb * enα[α] + Nn * ebα[α])
        dDst = st.dDds * sα[α] + st.dDdt * tα[α]
        db = -coom * (Np * enα[α] + Nn * epα[α])
        dc = -coom * (Np * ebα[α] + Nb * epα[α])
        ComplexF64[
            dD11 (daR-im*dDst) db
            (daR+im*dDst) dD22 dc
            db dc dD33
        ]
    end

    base = (; rho, theta, geo, st, DD, U, V = 0.5 * (D11 + D22), dUdr, dUdz, dUdkr, dUdkz, dUdom, G)
    need2nd || return base

    # ----- T2 Hessian over (r, z, kr, kz) -----
    # T2 = bilinear part (products of element FIRST derivatives — computed by
    # EXACT polynomial extraction from determinant samples, see _detU_bilinear)
    # + adjugate part tr(adj(DD)·∂²DD) carrying the element second derivatives.
    adjA = _adj3(DD)
    # element second derivatives: (s,t) block via S/D/P, k-block constants
    function _hd_pos(du2, duv, dv2, i, j)      # d²(·)/dz_i dz_j from (s,t) 2nd derivs
        return du2 * sα[i] * sα[j] + duv * (sα[i] * tα[j] + tα[i] * sα[j]) + dv2 * tα[i] * tα[j]
    end
    T2 = zeros(4, 4)
    for i = 1:4, j = i:4
        # bilinear part
        bil = _detU_bilinear(DD, U, G[i], G[j], i == j)
        # adjugate part: HD = d²DD/dz_i dz_j
        if i <= 2 && j <= 2
            hd11 = -_hd_pos(st.dSds2, st.dSdst, st.dSdt2, i, j)
            hd33 = -_hd_pos(st.dPds2, st.dPdst, st.dPdt2, i, j)
            hdDst = _hd_pos(st.dDds2, st.dDdst, st.dDdt2, i, j)
            HD = ComplexF64[
                hd11 (-im*hdDst) 0.0
                (im*hdDst) hd11 0.0
                0.0 0.0 hd33
            ]
            T2[i, j] = bil + real(sum(adjA[a, b] * HD[b, a] for a = 1:3, b = 1:3))
        elseif i >= 3 && j >= 3
            hd11 = 2 * coomsq * (ebα[i] * ebα[j] + epα[i] * epα[j])
            hd22 = 2 * coomsq * (enα[i] * enα[j] + epα[i] * epα[j])
            hd33 = 2 * coomsq * (enα[i] * enα[j] + ebα[i] * ebα[j])
            hdaR = -coomsq * (enα[i] * ebα[j] + ebα[i] * enα[j])
            hdb = -coomsq * (enα[i] * epα[j] + epα[i] * enα[j])
            hdc = -coomsq * (ebα[i] * epα[j] + epα[i] * ebα[j])
            HD = ComplexF64[
                hd11 hdaR hdb
                hdaR hd22 hdc
                hdb hdc hd33
            ]
            T2[i, j] = bil + real(sum(adjA[a, b] * HD[b, a] for a = 1:3, b = 1:3))
        else
            # mixed position-k: element second derivatives vanish in the
            # dispertok variable structure (2×2 parity)
            T2[i, j] = bil
        end
        T2[j, i] = T2[i, j]
    end
    return (; base..., T2)
end

# EXACT second-order coefficient extraction for det(A + q·Gi + p·Gj): the
# determinant is a (bi)cubic polynomial, so symmetric sampling recovers the
# quadratic/mixed coefficients exactly (only floating-point rounding, tamed by
# scaling the increments to the matrix norm). For 2×2 this reduces exactly to
# upstream's H = products-of-first-derivatives formulas.
function _detU_bilinear(
    A::Matrix{ComplexF64},
    f0::Float64,
    Gi::Matrix{ComplexF64},
    Gj::Matrix{ComplexF64},
    diag::Bool,
)
    normA = max(norm(A), 1e-300)
    if diag
        λ = normA / max(norm(Gi), 1e-300)
        f⁺ = _hdet(A .+ λ .* Gi)
        f⁻ = _hdet(A .- λ .* Gi)
        return (f⁺ + f⁻ - 2 * f0) / λ^2
    else
        λi = normA / max(norm(Gi), 1e-300)
        λj = normA / max(norm(Gj), 1e-300)
        f⁺⁺ = _hdet(A .+ λi .* Gi .+ λj .* Gj)
        f⁺⁻ = _hdet(A .+ λi .* Gi .- λj .* Gj)
        f⁻⁺ = _hdet(A .- λi .* Gi .+ λj .* Gj)
        f⁻⁻ = _hdet(A .- λi .* Gi .- λj .* Gj)
        return (f⁺⁺ - f⁺⁻ - f⁻⁺ + f⁻⁻) / (4 * λi * λj)
    end
end

# first directional derivative of det(A + q·G) at q = 0, exact for 2×2 and 3×3
# (the cubic coefficient det(G) is subtracted for 3×3)
function _detU_directional(A::Matrix{ComplexF64}, G::Matrix{ComplexF64})
    normA = max(norm(A), 1e-300)
    λ = normA / max(norm(G), 1e-300)
    f⁺ = _hdet(A .+ λ .* G)
    f⁻ = _hdet(A .- λ .* G)
    c3 = size(A, 1) == 3 ? _hdet(Matrix(λ .* G)) : 0.0
    return (f⁺ - f⁻) / (2 * λ) - c3 / λ
end

# Osculating-plane quantities at a point: unit velocity/acceleration directions
# (normalized by the symplectic area A), H Hessian of det-U along (q, p) with
# elements linearized (the Tracy-Kaufman quadratic-form structure — for 2×2
# identical to upstream's product formulas by the polynomial identities above),
# and the saddle displacement.
function _osculating(core, eqv::Vector{Float64}, epv::Vector{Float64})
    A = core.DD
    Gq = eqv[1] .* core.G[1] .+ eqv[2] .* core.G[2] .+ eqv[3] .* core.G[3] .+ eqv[4] .* core.G[4]
    Gp = epv[1] .* core.G[1] .+ epv[2] .* core.G[2] .+ epv[3] .* core.G[3] .+ epv[4] .* core.G[4]
    dUdq = _detU_directional(A, Gq)
    dUdp = _detU_directional(A, Gp)
    H11 = _detU_bilinear(A, core.U, Gq, Gq, true)
    H22 = _detU_bilinear(A, core.U, Gp, Gp, true)
    H12 = _detU_bilinear(A, core.U, Gq, Gp, false)
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
    n = size(gd, 1)
    new = ones(ComplexF64, n)
    old = zeros(ComplexF64, n)
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
    gdalf::Vector{Float64}      # ∇(uncoupled α-Hamiltonian) at the saddle
    gdlam::Vector{Float64}      # ∇(uncoupled λ-Hamiltonian) at the saddle
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
            zeros(4),
            zeros(4),
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
            zeros(4),
            zeros(4),
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
        coeff.gdalf,
        coeff.gdlam,
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
        gdalf = zeros(4),
        gdlam = zeros(4),
    )
    vp = F.vectors[:, ind[4]]
    vm = F.vectors[:, ind[1]]

    # Effective 2×2 coupling problem. For :cld2x2 this is the tensor itself
    # (MATLAB-pinned path). For :cld3x3 the conversion couples the TWO
    # near-null branches only; the far-off-shell third branch (|λ₃| ≫ 0, the
    # electrostatic polarization at ICRF frequencies) would dominate the
    # polarization iteration if kept, so the problem is reduced to the
    # near-null eigen-subspace at the saddle (the Tracy–Kaufman 2-D crossing
    # subspace): D̃ = Vᴴ·DD·V, G̃[α] = Vᴴ·G[α]·V with V the two smallest-|λ|
    # eigenvectors. When the third row/column decouples exactly this reduces
    # to the 2×2 objects identically.
    local Deff::Matrix{ComplexF64}
    local Geff::NTuple{4,Matrix{ComplexF64}}
    if size(cst.DD, 1) == 2
        Deff = cst.DD
        Geff = cst.G
    else
        FE = eigen(Hermitian(cst.DD))
        order = sortperm(abs.(FE.values))
        Vsub = FE.vectors[:, order[1:2]]
        Deff = Vsub' * cst.DD * Vsub
        Geff = ntuple(α -> Matrix{ComplexF64}(Vsub' * cst.G[α] * Vsub), 4)
    end

    # uncoupled polarizations from the directional tensor gradients. The
    # MATLAB column-major reshape of the row layout [gD11·v gD12·v cgD12·v
    # gD22·v] is exactly transpose(Σ_α v_α G̃_α) — verified against the
    # original code for 2×2 at the reference saddle.
    _gdmat(v) = transpose(v[1] .* Geff[1] .+ v[2] .* Geff[2] .+ v[3] .* Geff[3] .+ v[4] .* Geff[4])
    pol1 = _power_iterate(Matrix(_gdmat(vp)))
    pol2 = _power_iterate(Matrix(_gdmat(vm)))

    # gradient of the uncoupled eikonal Hamiltonians in the conjugate-
    # polarization convention: gdalf_α = Re(eᵀ·G̃_α·ē) (2×2-verified form)
    _gdpol(e) = Float64[real(transpose(e) * (Geff[α] * conj(e))) for α = 1:4]
    gdalf = _gdpol(pol1)
    gdlam = _gdpol(pol2)
    braket = transpose(gdalf) * J4 * gdlam
    # upstream evaluates pol1'·reshape(DD,2,2)·pol2 with MATLAB's column-major
    # reshape, i.e. the TRANSPOSED tensor — consistent with the transposed gd
    # matrices the polarizations are iterated on (conjugate-polarization
    # convention throughout). Verified against the original code at the saddle.
    eta = (pol1' * transpose(Deff) * pol2) / sqrt(complex(braket))
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
    return (;
        hyperbola_ok,
        eta2,
        tau,
        beta = ComplexF64(beta),
        transmitted,
        transmitted_ok,
        gdalf,
        gdlam,
    )
end
