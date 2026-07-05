# raycon_dispersion.jl — cold-plasma dispersion via the near-zero eigenvalue of
# the Hermitian dispersion tensor (port of disp_eig.m; Stix-element section
# shared with dispertok.m).
#
# The traced Hamiltonian is U(z), the eigenvalue of the local dispersion tensor
# DD(z) closest to zero; the ray obeys dz/dσ = J·∇U with J = [0 I; −I 0] over
# (r, z | kr, kz) and σ the ray parameter (the physical time direction is
# sign(∂U/∂ω), `dUdomega`). Spatial derivatives of U include the
# field-curvature corrections from the (e_n, e_b, e_p) basis rotation.

@inline function _poloidal_coords(eq::SolovevEquilibrium, r::Float64, z::Float64)
    rho = sqrt((r - eq.r0)^2 + z^2)
    theta = atan(eq.r0 - r, -z) + π / 2
    return rho, theta
end

# ---------------------------------------------------------------- Stix layer

# Local plasma parameters + Stix elements with 1st AND 2nd derivatives.
# The 2nd-derivative block is a handful of scalar ops per species, so it is
# always computed: a single concrete NamedTuple return type keeps every caller
# type-stable (a runtime-flag-dependent return type is a JET/type-instability
# trap).
function _stix_local(prob::RayconProblem, geo)
    cnst = prob.cnst
    om = prob.omega
    om2 = om^2
    prof = plasma_profiles(prob, geo.sflx)
    ns = length(prob.amass)
    omp2 = Vector{Float64}(undef, ns)
    omc = Vector{Float64}(undef, ns)
    for k = 1:ns
        q = prob.acharge[k] * cnst.e
        m = prob.amass[k] * cnst.mp
        omp2[k] = q^2 * prof.n[k] / (m * cnst.eps0)
        omc[k] = geo.b * q / m                     # signed cyclotron frequency
    end
    omc2 = omc .^ 2
    caoc2 = 1 / sum(omp2 ./ omc2)                  # (c_A / c)²
    dLNomcds = geo.dbds / geo.b
    dLNomcdt = geo.dbdt / geo.b

    omc2Mom2 = omc2 .- om2
    Si = omp2 ./ omc2Mom2
    S = 1 + sum(Si)
    Di = (omc ./ om) .* Si
    D = sum(Di)
    Pi = omp2 ./ om2
    P = 1 - sum(Pi)

    dLNSids = prof.dLNnds .- 2 .* omc2 ./ omc2Mom2 .* dLNomcds
    dLNSidt = -2 .* omc2 ./ omc2Mom2 .* dLNomcdt
    dLNSidom = 2 * om ./ omc2Mom2
    dSds = sum(Si .* dLNSids)
    dSdt = sum(Si .* dLNSidt)
    dSdom = sum(Si .* dLNSidom)
    dLNDids = dLNSids .+ dLNomcds
    dLNDidt = dLNSidt .+ dLNomcdt
    dLNDidom = (3 * om2 .- omc2) ./ (om .* omc2Mom2)
    dDds = sum(Di .* dLNDids)
    dDdt = sum(Di .* dLNDidt)
    dDdom = sum(Di .* dLNDidom)
    dPds = -sum(Pi .* prof.dLNnds)
    dPdt = 0.0
    dPdom = sum(Pi .* (2 / om))

    dLNomcds2 = geo.dbds2 / geo.b - (geo.dbds / geo.b)^2
    dLNomcdst = geo.dbdst / geo.b - geo.dbds * geo.dbdt / geo.b^2
    dLNomcdt2 = geo.dbdt2 / geo.b - (geo.dbdt / geo.b)^2
    omOM = om2 ./ omc2Mom2
    ocOM = omc2 ./ omc2Mom2
    dLNSids2 = 2 .* ocOM .* (2 .* omOM .* dLNomcds .* dLNomcds .- dLNomcds2) .+ prof.dLNnds2
    dLNSidst = 2 .* ocOM .* (2 .* omOM .* dLNomcds .* dLNomcdt .- dLNomcdst)
    dLNSidt2 = 2 .* ocOM .* (2 .* omOM .* dLNomcdt .* dLNomcdt .- dLNomcdt2)
    dLNDids2 = dLNomcds2 .+ dLNSids2
    dLNDidst = dLNomcdst .+ dLNSidst
    dLNDidt2 = dLNomcdt2 .+ dLNSidt2
    dSds2 = sum(Si .* (dLNSids2 .+ dLNSids .^ 2))
    dSdst = sum(Si .* (dLNSidst .+ dLNSids .* dLNSidt))
    dSdt2 = sum(Si .* (dLNSidt2 .+ dLNSidt .^ 2))
    dDds2 = sum(Di .* (dLNDids2 .+ dLNDids .^ 2))
    dDdst = sum(Di .* (dLNDidst .+ dLNDids .* dLNDidt))
    dDdt2 = sum(Di .* (dLNDidt2 .+ dLNDidt .^ 2))
    dPds2 = -sum(Pi .* (prof.dLNnds2 .+ prof.dLNnds .^ 2))
    dPdst = 0.0
    dPdt2 = 0.0
    return (;
        omp2,
        omc,
        omc2,
        caoc2,
        S,
        D,
        P,
        dSds,
        dSdt,
        dSdom,
        dDds,
        dDdt,
        dDdom,
        dPds,
        dPdt,
        dPdom,
        T = prof.T,
        n = prof.n,
        dSds2,
        dSdst,
        dSdt2,
        dDds2,
        dDdst,
        dDdt2,
        dPds2,
        dPdst,
        dPdt2,
    )
end

"""
    stix_elements(prob, r, z) -> (; S, D, P, omp2, omc, caoc2)

Cold-plasma Stix elements and local species frequencies at poloidal position
`(r, z)` [m] for the antenna frequency of `prob`.
"""
function stix_elements(prob::RayconProblem, r::Real, z::Real)
    rho, theta = _poloidal_coords(prob.eq, Float64(r), Float64(z))
    geo = magnetic_geometry(prob.eq, rho, theta)
    st = _stix_local(prob, geo)
    return (; st.S, st.D, st.P, st.omp2, st.omc, st.caoc2)
end

"""
    cyclotron_frequencies(prob, r, z) -> Vector{Float64}

Port of `dispertok(...,'Frq')`: for a 3-species plasma returns
`[f_ii, f_c2, f_c3]` — the ion–ion hybrid frequency and the two ion cyclotron
frequencies in Hz; otherwise the ion cyclotron frequencies (species 2:end).
"""
function cyclotron_frequencies(prob::RayconProblem, r::Real, z::Real)
    rho, theta = _poloidal_coords(prob.eq, Float64(r), Float64(z))
    geo = magnetic_geometry(prob.eq, rho, theta)
    st = _stix_local(prob, geo)
    ns = length(prob.amass)
    if ns == 3
        omii =
            sqrt(st.omp2[2] * st.omc[3]^2 + st.omp2[3] * st.omc[2]^2) /
            sqrt(st.omp2[2] + st.omp2[3])
        return [omii, st.omc[2], st.omc[3]] ./ (2π)
    else
        return st.omc[2:ns] ./ (2π)
    end
end

# ---------------------------------------------------------------- tensor core

# Everything disp_eig computes at one phase-space point. `need1st` adds the
# exact eigenvalue gradient dU_vec = [dU/dω, dU/dr, dU/dz, dU/dkr, dU/dkz];
# `need2nd` adds the 4×4 Hessian dU_mat over (r, z, kr, kz).
function _disp_core(
    prob::RayconProblem,
    r::Float64,
    zc::Float64,
    kr::Float64,
    kz::Float64;
    need1st::Bool = false,
    need2nd::Bool = false,
)
    prob.model in (:cld2x2, :cld3x3) || throw(
        ArgumentError("eigenvalue dispersion supports :cld2x2 and :cld3x3 (got $(prob.model))"),
    )
    rho, theta = _poloidal_coords(prob.eq, r, zc)
    geo = magnetic_geometry(prob.eq, rho, theta)
    st = _stix_local(prob, geo)
    om = prob.omega
    coom = prob.cnst.c / om
    coomsq = coom^2
    kf = prob.kphi

    kn = kr * geo.ener + kf * geo.enef + kz * geo.enez
    kb = kr * geo.eber + kf * geo.ebef + kz * geo.ebez
    kp = kr * geo.eper + kf * geo.epef + kz * geo.epez
    Nn = coom * kn
    Nb = coom * kb
    Np = coom * kp
    Nn2 = Nn^2
    Nb2 = Nb^2
    Np2 = Np^2
    Nr = coom * kr
    Nf = coom * kf
    Nz = coom * kz

    D11 = Nb2 + Np2 - st.S
    D12 = -Nn * Nb - im * st.D
    D22 = Nn2 + Np2 - st.S
    cD12 = conj(D12)

    is3x3 = prob.model === :cld3x3
    # third-row/column elements: defined unconditionally (zero for 2×2) so
    # every later `is3x3` block sees definitely-assigned, concretely-typed locals
    D13 = zero(ComplexF64)
    D23 = zero(ComplexF64)
    D33 = 0.0
    if is3x3
        D13 = -Nn * Np + 0im
        D23 = -Nb * Np + 0im
        D33 = Nn2 + Nb2 - st.P
        DD = ComplexF64[D11 D12 D13; cD12 D22 D23; conj(D13) conj(D23) D33]
    else
        DD = ComplexF64[D11 D12; cD12 D22]
    end
    cD13 = conj(D13)
    cD23 = conj(D23)
    F = eigen(Hermitian(DD))
    ind = argmin(abs.(F.values))
    U = F.values[ind]
    pol = F.vectors[:, ind]
    mon2 = if is3x3
        abs((D22 * D33 - D23 * cD23) + (D11 * D33 - D13 * cD13) + (D11 * D22 - D12 * cD12))
    else
        abs(D11 + D22)
    end

    base = (; rho, theta, geo, st, kn, kb, kp, Nn, Nb, Np, DD, U, pol, mon2)
    need1st || return base

    # ----- first derivatives of the tensor elements -----
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

    g = geo
    dD11dkr = dD11dkn * g.ener + dD11dkb * g.eber + dD11dkp * g.eper
    dD11dkz = dD11dkn * g.enez + dD11dkb * g.ebez + dD11dkp * g.epez
    dD12dkr = dD12dkn * g.ener + dD12dkb * g.eber + dD12dkp * g.eper
    dD12dkz = dD12dkn * g.enez + dD12dkb * g.ebez + dD12dkp * g.epez
    dD22dkr = dD22dkn * g.ener + dD22dkb * g.eber + dD22dkp * g.eper
    dD22dkz = dD22dkn * g.enez + dD22dkb * g.ebez + dD22dkp * g.epez

    dD11dom = -2 / om * (Nb2 + Np2) - st.dSdom
    # UPSTREAM BUG FIX: D12 = −Nn·Nb − iD and N ∝ 1/ω give ∂(−NnNb)/∂ω =
    # +2NnNb/ω; disp_eig.m/dispertok.m write −2NnNb/ω (sign copied from the
    # diagonal entries). Verified against finite differences of U(ω).
    dD12dom = +2 / om * (Nn * Nb) - im * st.dDdom
    dD22dom = -2 / om * (Nn2 + Np2) - st.dSdom

    # curvature corrections: derivatives of the projected refraction indices
    dNndr = Nr * g.denerdr + Nf * g.denefdr + Nz * g.denezdr
    dNndz = Nr * g.denerdz + Nf * g.denefdz + Nz * g.denezdz
    dNbdr = Nr * g.deberdr + Nf * g.debefdr + Nz * g.debezdr
    dNbdz = Nr * g.deberdz + Nf * g.debefdz + Nz * g.debezdz
    dNpdr = Nr * g.deperdr + Nf * g.depefdr + Nz * g.depezdr
    dNpdz = Nr * g.deperdz + Nf * g.depefdz + Nz * g.depezdz

    dD11dr = dD11ds * g.dsdr + dD11dt * g.dtdr + 2 * (Nb * dNbdr + Np * dNpdr)
    dD11dz = dD11ds * g.dsdz + dD11dt * g.dtdz + 2 * (Nb * dNbdz + Np * dNpdz)
    dD12dr = dD12ds * g.dsdr + dD12dt * g.dtdr - Nn * dNbdr - Nb * dNndr
    dD12dz = dD12ds * g.dsdz + dD12dt * g.dtdz - Nn * dNbdz - Nb * dNndz
    dD22dr = dD22ds * g.dsdr + dD22dt * g.dtdr + 2 * (Nn * dNndr + Np * dNpdr)
    dD22dz = dD22ds * g.dsdz + dD22dt * g.dtdz + 2 * (Nn * dNndz + Np * dNpdz)

    dD11_vec = ComplexF64[dD11dom, dD11dr, dD11dz, dD11dkr, dD11dkz]
    dD12_vec = ComplexF64[dD12dom, dD12dr, dD12dz, dD12dkr, dD12dkz]
    dD22_vec = ComplexF64[dD22dom, dD22dr, dD22dz, dD22dkr, dD22dkz]

    local dU_vec::Vector{Float64}
    # 3×3-only quantities, definitely assigned for type stability (zero/unused
    # in the 2×2 branch)
    dD13_vec = zeros(ComplexF64, 5)
    dD23_vec = zeros(ComplexF64, 5)
    dD33_vec = zeros(ComplexF64, 5)
    U1 = 0.0
    U2 = 0.0
    U3 = 0.0
    subsum3 = 0.0
    if is3x3
        dD13ds = 0.0
        dD13dt = 0.0
        dD23ds = 0.0
        dD23dt = 0.0
        dD33ds = -st.dPds
        dD33dt = -st.dPdt
        dD13dkn = -Np * coom
        dD13dkb = 0.0
        dD13dkp = -Nn * coom
        dD23dkn = 0.0
        dD23dkb = -Np * coom
        dD23dkp = -Nb * coom
        dD33dkn = 2 * Nn * coom
        dD33dkb = 2 * Nb * coom
        dD33dkp = 0.0
        dD13dkr = dD13dkn * g.ener + dD13dkb * g.eber + dD13dkp * g.eper
        dD13dkz = dD13dkn * g.enez + dD13dkb * g.ebez + dD13dkp * g.epez
        dD23dkr = dD23dkn * g.ener + dD23dkb * g.eber + dD23dkp * g.eper
        dD23dkz = dD23dkn * g.enez + dD23dkb * g.ebez + dD23dkp * g.epez
        dD33dkr = dD33dkn * g.ener + dD33dkb * g.eber + dD33dkp * g.eper
        dD33dkz = dD33dkn * g.enez + dD33dkb * g.ebez + dD33dkp * g.epez
        # same upstream sign bug as dD12dom (elements are −Nn·Np, −Nb·Np)
        dD13dom = +2 / om * (Nn * Np)
        dD23dom = +2 / om * (Nb * Np)
        dD33dom = -2 / om * (Nn2 + Nb2) - st.dPdom
        dD13dr = dD13ds * g.dsdr + dD13dt * g.dtdr - dNndr * Np - Nn * dNpdr
        dD13dz = dD13ds * g.dsdz + dD13dt * g.dtdz - dNndz * Np - Nn * dNpdz
        dD23dr = dD23ds * g.dsdr + dD23dt * g.dtdr - dNbdr * Np - Nb * dNpdr
        dD23dz = dD23ds * g.dsdz + dD23dt * g.dtdz - dNbdz * Np - Nb * dNpdz
        dD33dr = dD33ds * g.dsdr + dD33dt * g.dtdr + 2 * (Nn * dNndr + Nb * dNbdr)
        dD33dz = dD33ds * g.dsdz + dD33dt * g.dtdz + 2 * (Nn * dNndz + Nb * dNbdz)
        dD13_vec = ComplexF64[dD13dom, dD13dr, dD13dz, dD13dkr, dD13dkz]
        dD23_vec = ComplexF64[dD23dom, dD23dr, dD23dz, dD23dkr, dD23dkz]
        dD33_vec = ComplexF64[dD33dom, dD33dr, dD33dz, dD33dkr, dD33dkz]

        U1 = real((D22 - U) * (D33 - U) - D23 * cD23)
        U2 = real((D11 - U) * (D33 - U) - D13 * cD13)
        U3 = real((D11 - U) * (D22 - U) - D12 * cD12)
        subsum3 = U1 + U2 + U3
        dU_vec =
            real.(
                (
                    dD11_vec .* U1 .- (D11 - U) .* 2 .* real.(cD23 .* dD23_vec) .+ dD22_vec .* U2 .-
                    (D22 - U) .* 2 .* real.(cD13 .* dD13_vec) .+ dD33_vec .* U3 .-
                    (D33 - U) .* 2 .* real.(cD12 .* dD12_vec)
                ) ./ subsum3 .+
                2 .*
                real.(
                    D12 .* D23 .* conj.(dD13_vec) .+ D12 .* dD23_vec .* cD13 .+
                    dD12_vec .* D23 .* cD13,
                ) ./ subsum3,
            )
    else
        subsum = real(D11 + D22 - 2 * U)
        dU_vec =
            real.(
                (dD11_vec .* (D22 - U) .+ dD22_vec .* (D11 - U)) ./ subsum .-
                2 .* real.(cD12 .* dD12_vec) ./ subsum,
            )
    end

    # toroidal-wavenumber derivative ∂U/∂kφ (needed by the eikonal-phase
    # transport; kφ is constant along rays but the phase advances with it)
    dD11dkf = dD11dkn * g.enef + dD11dkb * g.ebef + dD11dkp * g.epef
    dD12dkf = dD12dkn * g.enef + dD12dkb * g.ebef + dD12dkp * g.epef
    dD22dkf = dD22dkn * g.enef + dD22dkb * g.ebef + dD22dkp * g.epef
    local dUdkf::Float64
    if is3x3
        dD13dkf = (-Np * coom) * g.enef + 0.0 * g.ebef + (-Nn * coom) * g.epef
        dD23dkf = 0.0 * g.enef + (-Np * coom) * g.ebef + (-Nb * coom) * g.epef
        dD33dkf = (2 * Nn * coom) * g.enef + (2 * Nb * coom) * g.ebef + 0.0 * g.epef
        U1v = real((D22 - U) * (D33 - U) - D23 * cD23)
        U2v = real((D11 - U) * (D33 - U) - D13 * cD13)
        U3v = real((D11 - U) * (D22 - U) - D12 * cD12)
        ss3 = U1v + U2v + U3v
        dUdkf = real(
            (
                dD11dkf * U1v - (D11 - U) * 2 * real(cD23 * dD23dkf) + dD22dkf * U2v -
                (D22 - U) * 2 * real(cD13 * dD13dkf) + dD33dkf * U3v -
                (D33 - U) * 2 * real(cD12 * dD12dkf)
            ) / ss3 +
            2 * real(D12 * D23 * conj(dD13dkf) + D12 * dD23dkf * cD13 + dD12dkf * D23 * cD13) / ss3,
        )
    else
        ss2 = real(D11 + D22 - 2 * U)
        dUdkf =
            real((dD11dkf * (D22 - U) + dD22dkf * (D11 - U)) / ss2 - 2 * real(cD12 * dD12dkf) / ss2)
    end

    lvl1 = (;
        base...,
        dD11_vec,
        dD12_vec,
        dD22_vec,
        dU_vec,
        dUdkf,
        dNndr,
        dNndz,
        dNbdr,
        dNbdz,
        dNpdr,
        dNpdz,
    )
    need2nd || return lvl1

    # ----- second derivatives (dU_mat over (r, z, kr, kz)) -----
    dSdr2 =
        g.dsdr2 * st.dSds +
        g.dtdr2 * st.dSdt +
        g.dsdr * (g.dsdr * st.dSds2 + g.dtdr * st.dSdst) +
        g.dtdr * (g.dsdr * st.dSdst + g.dtdr * st.dSdt2)
    dSdrz =
        g.dsdrz * st.dSds +
        g.dtdrz * st.dSdt +
        g.dsdz * (g.dsdr * st.dSds2 + g.dtdr * st.dSdst) +
        g.dtdz * (g.dsdr * st.dSdst + g.dtdr * st.dSdt2)
    dSdz2 =
        g.dsdz2 * st.dSds +
        g.dtdz2 * st.dSdt +
        g.dsdz * (g.dsdz * st.dSds2 + g.dtdz * st.dSdst) +
        g.dtdz * (g.dsdz * st.dSdst + g.dtdz * st.dSdt2)
    dDdr2 =
        g.dsdr2 * st.dDds +
        g.dtdr2 * st.dDdt +
        g.dsdr * (g.dsdr * st.dDds2 + g.dtdr * st.dDdst) +
        g.dtdr * (g.dsdr * st.dDdst + g.dtdr * st.dDdt2)
    dDdrz =
        g.dsdrz * st.dDds +
        g.dtdrz * st.dDdt +
        g.dsdz * (g.dsdr * st.dDds2 + g.dtdr * st.dDdst) +
        g.dtdz * (g.dsdr * st.dDdst + g.dtdr * st.dDdt2)
    dDdz2 =
        g.dsdz2 * st.dDds +
        g.dtdz2 * st.dDdt +
        g.dsdz * (g.dsdz * st.dDds2 + g.dtdz * st.dDdst) +
        g.dtdz * (g.dsdz * st.dDdst + g.dtdz * st.dDdt2)

    dNndr2 = Nr * g.denerdr2 + Nf * g.denefdr2 + Nz * g.denezdr2
    dNndrz = Nr * g.denerdrz + Nf * g.denefdrz + Nz * g.denezdrz
    dNndz2 = Nr * g.denerdz2 + Nf * g.denefdz2 + Nz * g.denezdz2
    dNbdr2 = Nr * g.deberdr2 + Nf * g.debefdr2 + Nz * g.debezdr2
    dNbdrz = Nr * g.deberdrz + Nf * g.debefdrz + Nz * g.debezdrz
    dNbdz2 = Nr * g.deberdz2 + Nf * g.debefdz2 + Nz * g.debezdz2
    dNpdr2 = Nr * g.deperdr2 + Nf * g.depefdr2 + Nz * g.depezdr2
    dNpdrz = Nr * g.deperdrz + Nf * g.depefdrz + Nz * g.depezdrz
    dNpdz2 = Nr * g.deperdz2 + Nf * g.depefdz2 + Nz * g.depezdz2

    dD11dr2 = 2 * (dNbdr * dNbdr + Nb * dNbdr2 + dNpdr * dNpdr + Np * dNpdr2) - dSdr2
    dD11drz = 2 * (dNbdr * dNbdz + Nb * dNbdrz + dNpdr * dNpdz + Np * dNpdrz) - dSdrz
    dD11dz2 = 2 * (dNbdz * dNbdz + Nb * dNbdz2 + dNpdz * dNpdz + Np * dNpdz2) - dSdz2
    dD12dr2 = -(Nn * dNbdr2 + 2 * dNndr * dNbdr + Nb * dNndr2) - im * dDdr2
    dD12drz = -(Nn * dNbdrz + dNndr * dNbdz + dNndz * dNbdr + Nb * dNndrz) - im * dDdrz
    dD12dz2 = -(Nn * dNbdz2 + 2 * dNndz * dNbdz + Nb * dNndz2) - im * dDdz2
    dD22dr2 = 2 * (dNndr * dNndr + Nn * dNndr2 + dNpdr * dNpdr + Np * dNpdr2) - dSdr2
    dD22drz = 2 * (dNndr * dNndz + Nn * dNndrz + dNpdr * dNpdz + Np * dNpdrz) - dSdrz
    dD22dz2 = 2 * (dNndz * dNndz + Nn * dNndz2 + dNpdz * dNpdz + Np * dNpdz2) - dSdz2

    dD11drkr = 2 * coom * (g.eber * dNbdr + Nb * g.deberdr + g.eper * dNpdr + Np * g.deperdr)
    dD11drkz = 2 * coom * (g.ebez * dNbdr + Nb * g.debezdr + g.epez * dNpdr + Np * g.depezdr)
    dD11dzkr = 2 * coom * (g.eber * dNbdz + Nb * g.deberdz + g.eper * dNpdz + Np * g.deperdz)
    dD11dzkz = 2 * coom * (g.ebez * dNbdz + Nb * g.debezdz + g.epez * dNpdz + Np * g.depezdz)
    dD12drkr = -coom * (g.denerdr * Nb + dNndr * g.eber + g.ener * dNbdr + Nn * g.deberdr)
    dD12drkz = -coom * (g.denezdr * Nb + dNndr * g.ebez + g.enez * dNbdr + Nn * g.debezdr)
    dD12dzkr = -coom * (g.denerdz * Nb + dNndz * g.eber + g.ener * dNbdz + Nn * g.deberdz)
    dD12dzkz = -coom * (g.denezdz * Nb + dNndz * g.ebez + g.enez * dNbdz + Nn * g.debezdz)
    dD22drkr = 2 * coom * (g.ener * dNndr + Nn * g.denerdr + g.eper * dNpdr + Np * g.deperdr)
    dD22drkz = 2 * coom * (g.enez * dNndr + Nn * g.denezdr + g.epez * dNpdr + Np * g.depezdr)
    dD22dzkr = 2 * coom * (g.ener * dNndz + Nn * g.denerdz + g.eper * dNpdz + Np * g.deperdz)
    dD22dzkz = 2 * coom * (g.enez * dNndz + Nn * g.denezdz + g.epez * dNpdz + Np * g.depezdz)

    coomsq2 = 2 * coomsq
    dD11dkr2 = coomsq2 * (g.eber^2 + g.eper^2)
    dD11dkz2 = coomsq2 * (g.ebez^2 + g.epez^2)
    dD11dkrkz = coomsq2 * (g.eber * g.ebez + g.eper * g.epez)
    # UPSTREAM BUG FIX: disp_eig.m swaps the dkz² and dkr·dkz expressions for
    # the off-diagonal elements (their "dD12dkz2" holds the mixed derivative
    # and vice versa). Verified against finite differences of the Hessian.
    dD12dkr2 = -coomsq * (g.ener * g.eber + g.ener * g.eber)
    dD12dkz2 = -coomsq * (g.enez * g.ebez + g.enez * g.ebez)
    dD12dkrkz = -coomsq * (g.ener * g.ebez + g.enez * g.eber)
    dD22dkr2 = coomsq2 * (g.ener^2 + g.eper^2)
    dD22dkz2 = coomsq2 * (g.enez^2 + g.epez^2)
    dD22dkrkz = coomsq2 * (g.ener * g.enez + g.eper * g.epez)

    _sym4(a11, a12, a13, a14, a22, a23, a24, a33, a34, a44) = ComplexF64[
        a11 a12 a13 a14
        a12 a22 a23 a24
        a13 a23 a33 a34
        a14 a24 a34 a44
    ]
    dD11_mat = _sym4(
        dD11dr2,
        dD11drz,
        dD11drkr,
        dD11drkz,
        dD11dz2,
        dD11dzkr,
        dD11dzkz,
        dD11dkr2,
        dD11dkrkz,
        dD11dkz2,
    )
    dD12_mat = _sym4(
        dD12dr2,
        dD12drz,
        dD12drkr,
        dD12drkz,
        dD12dz2,
        dD12dzkr,
        dD12dzkz,
        dD12dkr2,
        dD12dkrkz,
        dD12dkz2,
    )
    dD22_mat = _sym4(
        dD22dr2,
        dD22drz,
        dD22drkr,
        dD22drkz,
        dD22dz2,
        dD22dzkr,
        dD22dzkz,
        dD22dkr2,
        dD22dkrkz,
        dD22dkz2,
    )

    v11 = dD11_vec[2:5]
    v12 = dD12_vec[2:5]
    v22 = dD22_vec[2:5]
    vU = ComplexF64.(dU_vec[2:5])

    local dU_mat::Matrix{Float64}
    if is3x3
        dPdr2 =
            g.dsdr2 * st.dPds +
            g.dtdr2 * st.dPdt +
            g.dsdr * (g.dsdr * st.dPds2 + g.dtdr * st.dPdst) +
            g.dtdr * (g.dsdr * st.dPdst + g.dtdr * st.dPdt2)
        dPdrz =
            g.dsdrz * st.dPds +
            g.dtdrz * st.dPdt +
            g.dsdz * (g.dsdr * st.dPds2 + g.dtdr * st.dPdst) +
            g.dtdz * (g.dsdr * st.dPdst + g.dtdr * st.dPdt2)
        dPdz2 =
            g.dsdz2 * st.dPds +
            g.dtdz2 * st.dPdt +
            g.dsdz * (g.dsdz * st.dPds2 + g.dtdz * st.dPdst) +
            g.dtdz * (g.dsdz * st.dPdst + g.dtdz * st.dPdt2)
        # D13/D23/D33/cD13/cD23 are the definitely-assigned locals from above
        dD13dr2 = -(Nn * dNpdr2 + 2 * dNndr * dNpdr + Np * dNndr2)
        dD13drz = -(Nn * dNpdrz + dNndr * dNpdz + dNndz * dNpdr + Np * dNndrz)
        dD13dz2 = -(Nn * dNpdz2 + 2 * dNndz * dNpdz + Np * dNndz2)
        dD23dr2 = -(Nb * dNpdr2 + 2 * dNbdr * dNpdr + Np * dNbdr2)
        dD23drz = -(Nb * dNpdrz + dNbdr * dNpdz + dNbdz * dNpdr + Np * dNbdrz)
        dD23dz2 = -(Nb * dNpdz2 + 2 * dNbdz * dNpdz + Np * dNbdz2)
        dD33dr2 = 2 * (dNndr * dNndr + Nn * dNndr2 + dNbdr * dNbdr + Nb * dNbdr2) - dPdr2
        dD33drz = 2 * (dNndr * dNndz + Nn * dNndrz + dNbdr * dNbdz + Nb * dNbdrz) - dPdrz
        dD33dz2 = 2 * (dNndz * dNndz + Nn * dNndz2 + dNbdz * dNbdz + Nb * dNbdz2) - dPdz2
        dD13drkr = -coom * (g.denerdr * Np + dNndr * g.eper + g.ener * dNpdr + Nn * g.deperdr)
        dD13drkz = -coom * (g.denezdr * Np + dNndr * g.epez + g.enez * dNpdr + Nn * g.depezdr)
        dD13dzkr = -coom * (g.denerdz * Np + dNndz * g.eper + g.ener * dNpdz + Nn * g.deperdz)
        dD13dzkz = -coom * (g.denezdz * Np + dNndz * g.epez + g.enez * dNpdz + Nn * g.depezdz)
        dD23drkr = -coom * (g.deperdr * Nb + dNpdr * g.eber + g.eper * dNbdr + Np * g.deberdr)
        dD23drkz = -coom * (g.depezdr * Nb + dNpdr * g.ebez + g.epez * dNbdr + Np * g.debezdr)
        dD23dzkr = -coom * (g.deperdz * Nb + dNpdz * g.eber + g.eper * dNbdz + Np * g.deberdz)
        dD23dzkz = -coom * (g.depezdz * Nb + dNpdz * g.ebez + g.epez * dNbdz + Np * g.debezdz)
        dD33drkr = 2 * coom * (g.ener * dNndr + Nn * g.denerdr + g.eber * dNbdr + Nb * g.deberdr)
        dD33drkz = 2 * coom * (g.enez * dNndr + Nn * g.denezdr + g.ebez * dNbdr + Nb * g.debezdr)
        dD33dzkr = 2 * coom * (g.ener * dNndz + Nn * g.denerdz + g.eber * dNbdz + Nb * g.deberdz)
        dD33dzkz = 2 * coom * (g.enez * dNndz + Nn * g.denezdz + g.ebez * dNbdz + Nb * g.debezdz)
        # same upstream dkz²/dkr·dkz swap as dD12 (see above)
        dD13dkr2 = -coomsq * (g.ener * g.eper + g.ener * g.eper)
        dD13dkz2 = -coomsq * (g.enez * g.epez + g.enez * g.epez)
        dD13dkrkz = -coomsq * (g.ener * g.epez + g.enez * g.eper)
        dD23dkr2 = -coomsq * (g.eper * g.eber + g.eper * g.eber)
        dD23dkz2 = -coomsq * (g.epez * g.ebez + g.epez * g.ebez)
        dD23dkrkz = -coomsq * (g.eper * g.ebez + g.epez * g.eber)
        dD33dkr2 = coomsq2 * (g.ener^2 + g.eber^2)
        dD33dkz2 = coomsq2 * (g.enez^2 + g.ebez^2)
        dD33dkrkz = coomsq2 * (g.ener * g.enez + g.eber * g.ebez)
        dD13_mat = _sym4(
            dD13dr2,
            dD13drz,
            dD13drkr,
            dD13drkz,
            dD13dz2,
            dD13dzkr,
            dD13dzkz,
            dD13dkr2,
            dD13dkrkz,
            dD13dkz2,
        )
        dD23_mat = _sym4(
            dD23dr2,
            dD23drz,
            dD23drkr,
            dD23drkz,
            dD23dz2,
            dD23dzkr,
            dD23dzkz,
            dD23dkr2,
            dD23dkrkz,
            dD23dkz2,
        )
        dD33_mat = _sym4(
            dD33dr2,
            dD33drz,
            dD33drkr,
            dD33drkz,
            dD33dz2,
            dD33dzkr,
            dD33dzkz,
            dD33dkr2,
            dD33dkrkz,
            dD33dkz2,
        )
        v13 = dD13_vec[2:5]
        v23 = dD23_vec[2:5]
        v33 = dD33_vec[2:5]
        dU1_vec = (v22 .- vU) .* (D33 - U) .+ (v33 .- vU) .* (D22 - U) .- 2 .* real.(cD23 .* v23)
        dU2_vec = (v11 .- vU) .* (D33 - U) .+ (v33 .- vU) .* (D11 - U) .- 2 .* real.(cD13 .* v13)
        dU3_vec = (v22 .- vU) .* (D11 - U) .+ (v11 .- vU) .* (D22 - U) .- 2 .* real.(cD12 .* v12)
        dB_mat =
            2 .*
            real.(
                D23 .* (v12 * transpose(conj.(v13))) .+ D12 .* (v23 * transpose(conj.(v13))) .+
                D12 * D23 .* conj.(dD13_mat) .+ (v12 * transpose(v23)) .* cD13 .+
                D12 .* dD23_mat .* cD13 .+ D12 .* (conj.(v13) * transpose(v23)) .+
                dD12_mat .* (D23 * cD13) .+ (v23 * transpose(v12)) .* cD13 .+
                (conj.(v13) * transpose(v12)) .* D23,
            )
        M = dD11_mat .* U1 .+ dD22_mat .* U2 .+ dD33_mat .* U3
        M =
            M .+ dU1_vec * transpose(v11 .- vU) .+ dU2_vec * transpose(v22 .- vU) .+
            dU3_vec * transpose(v33 .- vU)
        M =
            M .- 2 .* ((v11 .- vU) * transpose(real.(cD23 .* v23))) .-
            2 .* ((v22 .- vU) * transpose(real.(cD13 .* v13))) .-
            2 .* ((v33 .- vU) * transpose(real.(cD12 .* v12)))
        M =
            M .- (D11 - U) .* 2 .* real.(cD23 .* dD23_mat .+ conj.(v23) * transpose(v23)) .-
            (D22 - U) .* 2 .* real.(cD13 .* dD13_mat .+ conj.(v13) * transpose(v13)) .-
            (D33 - U) .* 2 .* real.(cD12 .* dD12_mat .+ conj.(v12) * transpose(v12))
        M = M .+ dB_mat
        dU_mat = real.(M ./ subsum3)
        maxasym = maximum(abs.(dU_mat .- transpose(dU_mat)))
        maxasym <= 1e-10 * max(1.0, maximum(abs.(dU_mat))) ||
            error("raycon dispersion Hessian lost symmetry (asymmetry $maxasym)")
    else
        subsum = real(D11 + D22 - 2 * U)
        dU_mat =
            real.(
                (
                    dD11_mat .* (D22 - U) .+ dD22_mat .* (D11 - U) .+
                    (v11 .- vU) * transpose(v22 .- vU) .+ (v22 .- vU) * transpose(v11 .- vU) .-
                    2 .* real.(cD12 .* dD12_mat .+ conj.(v12) * transpose(v12))
                ) ./ subsum,
            )
    end
    return (; lvl1..., dU_mat)
end

@inline function _validated_state4(y)
    length(y) == 4 || throw(ArgumentError("phase-space state must be (r, z, kr, kz)"))
    all(isfinite, y) || throw(ArgumentError("phase-space state must be finite"))
    return Float64(y[1]), Float64(y[2]), Float64(y[3]), Float64(y[4])
end

"""
    dispersion_U(prob, y) -> Float64

Dispersion function `U` (eigenvalue of the local dispersion tensor closest to
zero; `disp_eig(...,'Dsp')`) at phase-space point `y = (r, z, kr, kz)`.
`U = 0` on the dispersion manifold.
"""
function dispersion_U(prob::RayconProblem, y)
    r, zc, kr, kz = _validated_state4(y)
    return Float64(_disp_core(prob, r, zc, kr, kz).U)
end

"""
    conversion_monitors(prob, y) -> (; mon1, mon2)

Caustic (`mon1`, identically 0 for 4-dim tracing — upstream parity) and mode-
conversion (`mon2`) monitors: `mon2 = |tr DD|` (2x2) or the sum of principal
2×2 minors (3x3). A local minimum of `mon2` along a ray signals an avoided
eigenvalue crossing, i.e. a mode-conversion candidate (`disp_eig(...,'Mon')`).
"""
function conversion_monitors(prob::RayconProblem, y)
    r, zc, kr, kz = _validated_state4(y)
    core = _disp_core(prob, r, zc, kr, kz)
    return (; mon1 = 0.0, mon2 = core.mon2)
end

"""
    polarization(prob, y) -> Vector{ComplexF64}

Polarization eigenvector (of the near-zero eigenvalue) in the local Stix basis
`(e_n, e_b, e_p)` (`disp_eig(...,'Pol')`). Defined up to a complex phase.
"""
function polarization(prob::RayconProblem, y)
    r, zc, kr, kz = _validated_state4(y)
    return _disp_core(prob, r, zc, kr, kz).pol
end

"""
    dUdomega(prob, y) -> Float64

`∂U/∂ω` (`dispertok(...,'Sgn')` via the eigenvalue form): its sign relates the
ray parameter σ to physical time, `dt/dσ = ∂U/∂ω`.
"""
function dUdomega(prob::RayconProblem, y)
    r, zc, kr, kz = _validated_state4(y)
    return _disp_core(prob, r, zc, kr, kz; need1st = true).dU_vec[1]
end

"""
    trajectory_rhs(prob, y) -> Vector{Float64}

Ray-equation right-hand side `dz/dσ = J·∇U` (`disp_eig(...,'Trj')`). For a
4-vector state returns the 4 phase-space rates; for a 20-vector state
`[z; vec(S)]` also evolves the 4×4 symplectic tangent map `dS/dσ = J·(∇∇U)·S`
used by conversion-curvature diagnostics.
"""
function trajectory_rhs(prob::RayconProblem, y::AbstractVector{<:Real})
    if length(y) == 4
        core = _disp_core(
            prob,
            Float64(y[1]),
            Float64(y[2]),
            Float64(y[3]),
            Float64(y[4]);
            need1st = true,
        )
        dU = core.dU_vec
        return [dU[4], dU[5], -dU[2], -dU[3]]
    elseif length(y) == 20
        core = _disp_core(
            prob,
            Float64(y[1]),
            Float64(y[2]),
            Float64(y[3]),
            Float64(y[4]);
            need1st = true,
            need2nd = true,
        )
        dU = core.dU_vec
        S = reshape(Float64.(y[5:20]), 4, 4)
        J = [0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; -1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0]
        dS = J * core.dU_mat * S
        return vcat([dU[4], dU[5], -dU[2], -dU[3]], vec(dS))
    else
        throw(ArgumentError("state must have 4 or 20 components (got $(length(y)))"))
    end
end

"""
    msw_dispersion(prob, y) -> Float64

Fast magnetosonic 1×1 approximation `U = ½(c²/c_A² − N²)` (`'Msw'`), used to
estimate the launch wavenumber before polishing on the full dispersion surface.
"""
function msw_dispersion(prob::RayconProblem, y)
    r, zc, kr, kz = _validated_state4(y)
    rho, theta = _poloidal_coords(prob.eq, r, zc)
    geo = magnetic_geometry(prob.eq, rho, theta)
    st = _stix_local(prob, geo)
    coom = prob.cnst.c / prob.omega
    kn = kr * geo.ener + prob.kphi * geo.enef + kz * geo.enez
    kb = kr * geo.eber + prob.kphi * geo.ebef + kz * geo.ebez
    kp = kr * geo.eper + prob.kphi * geo.epef + kz * geo.epez
    N2 = (coom * kn)^2 + (coom * kb)^2 + (coom * kp)^2
    return 0.5 * (1 / st.caoc2 - N2)
end

"""
    adjust_to_dispersion(prob, y; m=0.0) -> Vector{Float64}

Adjust a launch point onto the dispersion surface (`adjust_disp_m.m`): keeping
the poloidal wavenumber contribution `kθ = m/ρ` fixed, solve for the radial
component `kρ` such that `U = 0`, starting from
`kρ₀ = −√(kr² + kz² − (m/ρ)²)` (inward launch). Returns `[r, z, kr, kz]`.
"""
function adjust_to_dispersion(prob::RayconProblem, y; m::Real = 0.0)
    r, zc, kr, kz = _validated_state4(y)
    rho, theta = _poloidal_coords(prob.eq, r, zc)
    mrho = Float64(m) / rho
    arg = kr^2 + kz^2 - mrho^2
    arg >= 0 || throw(ArgumentError("|k|² < (m/ρ)²: no radial wavenumber at this launch"))
    krho0 = -sqrt(arg)
    st = sin(theta)
    ct = cos(theta)
    f = krho -> dispersion_U(prob, (r, zc, krho * ct - mrho * st, -(krho * st + mrho * ct)))
    krho = _fzero_near(f, krho0)
    return [r, zc, krho * ct - mrho * st, -(krho * st + mrho * ct)]
end
