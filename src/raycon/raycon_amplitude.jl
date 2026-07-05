# raycon_amplitude.jl — WKB amplitude transport along rays: focusing tensor,
# log-amplitude, eikonal phase, Maslov/caustic k-space switching, and
# collisionless damping with per-species power deposition.
#
# This completes the layer upstream RAYCON ships DISABLED (the 'Amp' RHS call
# is commented out in trajectory.m), so unlike the ray/conversion layers it is
# verified against constructed oracles rather than the MATLAB original:
#
#   * the focusing tensor obeys the WKB Riccati equation driven by the SAME
#     phase-space Hessian dU_mat as the 20-dim symplectic tangent map, and the
#     exact relation  W(σ) = (S_kx + S_kk·W₀)·(S_xx + S_xk·W₀)⁻¹,
#     lnE²(σ) = lnE²₀ − ln|det(S_xx + S_xk·W₀)|  is enforced in the tests
#     (valid in the full inhomogeneous field, through caustics);
#   * the damping layer implements upstream's draft formulas (fundamental and
#     second-harmonic cyclotron, Landau/TTMP) with the time-direction sign
#     fixed via ∂U/∂ω; energy bookkeeping (deposition = amplitude decrement)
#     holds identically at the equation level and is tested.
#
# Transport equations in the ray parameter σ (U dimensionless, ẋ = ∂U/∂k,
# k̇ = −∂U/∂x; blocks of dU_mat: A = U_xx, B = U_xk, C = U_kk over (r,z|kr,kz)):
#
#   x-space:  Ẇ  = −(A + B·W + W·Bᵀ + W·C·W),     d lnE²/dσ = −(tr B + tr(C·W))
#   k-space:  Ẇ̃  = +(C + Bᵀ·W̃ + W̃·B + W̃·A·W̃),   d lnE²/dσ = +(tr B + tr(A·W̃))
#   phase:    dΘ/dσ = k·∂U/∂k + kφ·∂U/∂kφ   (x-space)
#             dΘ̃/dσ = x·∂U/∂x + kφ·∂U/∂kφ   (k-space; upstream drops the kφ
#                                            term there — fixed, see notes)
#
# Maslov transform at caustics (port of ray.m 'caustic_list'): W ← W⁻¹ with
# amplitude factor 2π/√|det W| (x→k) resp. 1/(2π√|det W̃|) (k→x), phase index
# (π/4)Σsign(eig W) and ∓x·k, so a round trip is amplitude-neutral.

"""
    AmplitudeTrace

Result of [`integrate_ray_amplitude`](@ref): ray parameter `sigma`, phase-space
trajectory `y` (4×N), focusing tensor components `W` (3×N: W11, W12, W22 — the
k-space form W̃ while `inkspace` is true), log-amplitude-squared `lnE2`,
eikonal phase `phase`, per-species deposited energy `dep` (nspec×N, units of
the launch E²), `inkspace` flags, number of Maslov transforms `nmaslov`, and
the stop `status` (as [`integrate_ray`](@ref), plus `:caustic_failure`).
"""
struct AmplitudeTrace
    sigma::Vector{Float64}
    y::Matrix{Float64}
    W::Matrix{Float64}
    lnE2::Vector{Float64}
    phase::Vector{Float64}
    dep::Matrix{Float64}
    inkspace::Vector{Bool}
    nmaslov::Int
    status::Symbol
end

# per-species damping rates γ_s/ω (dimensionless; upstream dispertok 'Amp'
# lines 864-887: fundamental + 2nd-harmonic cyclotron, Landau/TTMP; the
# anomalous-viscosity channel carries an explicit 0 coefficient upstream and
# is omitted). Cold species (T = 0) do not damp.
function _damping_rates(prob::RayconProblem, core, kn::Float64, kb::Float64, kp::Float64)
    st = core.st
    cnst = prob.cnst
    om = prob.omega
    ns = length(prob.amass)
    Q = zeros(ns)
    kp == 0.0 && return Q                     # k∥ = 0: no resonant interaction
    # 3-component polarization in the Stix basis: eigen pol extended by the
    # electrostatic component −(D13·e1 + D23·e2)/D33 for the 2×2 model
    pol = core.pol
    local En::ComplexF64, Eb::ComplexF64, Ep::ComplexF64
    if length(pol) == 3
        En, Eb, Ep = pol[1], pol[2], pol[3]
    else
        coom = cnst.c / om
        Nn = coom * kn
        Nb = coom * kb
        Np = coom * kp
        D13 = -Nn * Np
        D23 = -Nb * Np
        D33 = Nn^2 + Nb^2 - st.P
        En, Eb = pol[1], pol[2]
        Ep = -(D13 * En + D23 * Eb) / D33
    end
    iomBp = im * (kn * Eb - kb * En)
    TTMP2 = abs2(iomBp)
    for s = 1:ns
        Ts = st.T[s]
        Ts > 0 || continue
        vth2 = 2000 * Ts * cnst.e / (prob.amass[s] * cnst.mp)
        vth = sqrt(vth2)
        omc = st.omc[s]
        kvt = abs(kp) * vth
        argp0 = om / kvt
        argm1 = (om - omc) / kvt
        argm2 = (om - 2 * omc) / kvt
        lrmr2 = vth2 / omc^2
        f1 = 2 * om * omc / (kp * vth2)
        land2 = abs2(f1 * Ep - iomBp)
        Em1 = En - im * Eb
        Em2 = (im * kn + kb) * Em1
        Qresn =
            exp(-argm1^2) * abs2(Em1) +
            lrmr2 * exp(-argm2^2) * abs2(Em2) +
            lrmr2 * exp(-argp0^2) * (TTMP2 + land2)
        Q[s] = sqrt(π) / (om * abs(kp)) * st.omp2[s] / vth * Qresn
    end
    return Q
end

# amplitude RHS: extended state (x(4), W(3), lnE2, phase, dep(nspec))
function _amp_rhs(prob::RayconProblem, u::Vector{Float64}, inkspace::Bool, damping::Bool)
    core = _disp_core(prob, u[1], u[2], u[3], u[4]; need1st = true, need2nd = true)
    dU = core.dU_vec
    H = core.dU_mat
    A = H[1:2, 1:2]
    B = H[1:2, 3:4]
    C = H[3:4, 3:4]
    W = [u[5] u[6]; u[6] u[7]]
    du = zeros(length(u))
    du[1] = dU[4]
    du[2] = dU[5]
    du[3] = -dU[2]
    du[4] = -dU[3]
    if inkspace
        dW = C .+ transpose(B) * W .+ W * B .+ W * A * W
        dlnE2 = tr(B) + tr(A * W)
        dphase = u[1] * dU[2] + u[2] * dU[3] + prob.kphi * core.dUdkf
    else
        dW = -(A .+ B * W .+ W * transpose(B) .+ W * C * W)
        dlnE2 = -(tr(B) + tr(C * W))
        dphase = u[3] * dU[4] + u[4] * dU[5] + prob.kphi * core.dUdkf
    end
    du[5] = dW[1, 1]
    du[6] = 0.5 * (dW[1, 2] + dW[2, 1])
    du[7] = dW[2, 2]
    du[8] = dlnE2
    du[9] = dphase
    if damping
        g = core.geo
        kn = u[3] * g.ener + prob.kphi * g.enef + u[4] * g.enez
        kb = u[3] * g.eber + prob.kphi * g.ebef + u[4] * g.ebez
        kp = u[3] * g.eper + prob.kphi * g.epef + u[4] * g.epez
        Q = _damping_rates(prob, core, kn, kb, kp)
        # decay along the traced (group-velocity) direction: dx/dσ = ∂U/∂k =
        # −v_g·∂U/∂ω, so σ-forward follows the energy flow for either sign of
        # ∂U/∂ω and the temporal rate dlnE²/dt = −2γ (γ_s = ω·Q_s) maps to
        # −2γ·|∂U/∂ω| per unit σ
        absdt = abs(dU[1])
        du[8] -= 2 * prob.omega * sum(Q) * absdt
        # per-species deposition rate = the corresponding amplitude decrement
        E2 = exp(u[8])
        for s = 1:length(Q)
            du[9+s] = 2 * prob.omega * Q[s] * E2 * absdt
        end
    end
    return du
end

# Maslov/Legendre transform at a caustic (port of ray.m 'caustic_list'):
# invert the focusing tensor, update amplitude and phase, toggle the space.
function _maslov!(u::Vector{Float64}, inkspace::Bool)
    W = [u[5] u[6]; u[6] u[7]]
    dW = det(W)
    abs(dW) > 0 || return (inkspace, false)
    Wi = inv(W)
    u[5] = Wi[1, 1]
    u[6] = 0.5 * (Wi[1, 2] + Wi[2, 1])
    u[7] = Wi[2, 2]
    idx = π / 4 * sum(sign, eigvals(Symmetric(W)))
    if inkspace
        f = 1 / (2π * sqrt(abs(dW)))          # k → x
        sgn = +1.0
    else
        f = 2π / sqrt(abs(dW))                # x → k
        sgn = -1.0
    end
    u[8] += log(f^2)
    u[9] += idx + sgn * (u[1] * u[3] + u[2] * u[4])
    return (!inkspace, true)
end

# Converted-ray focusing-tensor matching (port of dispertok.m 'Cnv' Amp block,
# lines 987-998): from the uncoupled-Hamiltonian gradients gdalf (incoming α
# branch, split as za = ∂/∂x, zb = ∂/∂k) and gdlam (outgoing λ branch,
# zc/zd), the outgoing Hessian is Slam = Slam0 + λ·Slam1 with λ fixed by
# matching the incoming curvature Salf. Degenerate geometry (zd ≈ 0 or a
# vanishing matching denominator) falls back to the incoming tensor.
function _slam_matching(
    Salf::AbstractMatrix{<:Real},
    gdalf::Vector{Float64},
    gdlam::Vector{Float64},
)
    za = gdalf[1:2]
    zb = gdalf[3:4]
    zc = gdlam[1:2]
    zd = gdlam[3:4]
    J2 = [0.0 1.0; -1.0 0.0]
    zdsq = dot(zd, zd)
    zdsq > 0 || return Matrix{Float64}(Salf)
    Slam0 = (dot(zc, zd) / zdsq^2) .* (zd * zd') .- (zc * zd' .+ zd * zc') ./ zdsq
    J2zd = J2 * zd
    Slam1 = J2zd * J2zd'
    zvc = J2 * (zc .+ Salf * zd)
    den = dot(zb, Slam1 * zvc)
    den != 0 || return Matrix{Float64}(Salf)
    λm = -dot(za .+ Slam0 * zb, zvc) / den
    return Slam0 .+ λm .* Slam1
end

"""
    antenna_focusing(prob, y) -> Matrix

Initial focusing tensor `W₀` for a wavefront aligned with the flux surfaces at
the launch point (port of `dispertok(...,'Ant')`): solves the 3×3 system for
`(∂²Θ/∂R², ∂²Θ/∂R∂Z, ∂²Θ/∂Z²)` from the flux-surface slope and the local ray
direction. Ill-defined where the surface tangent is vertical (`∂s/∂Z = 0`,
e.g. exactly on the midplane) — launch slightly off `theta = 0`.
"""
function antenna_focusing(prob::RayconProblem, y::AbstractVector{<:Real})
    length(y) == 4 || throw(ArgumentError("state must be (r, z, kr, kz)"))
    core =
        _disp_core(prob, Float64(y[1]), Float64(y[2]), Float64(y[3]), Float64(y[4]); need1st = true)
    g = core.geo
    abs(g.dsdz) > 1e-12 * abs(g.dsdr) || throw(
        ArgumentError(
            "flux-surface-aligned wavefront is ill-defined here " *
            "(∂s/∂Z ≈ 0); launch slightly off the midplane",
        ),
    )
    dzdr = -g.dsdr / g.dsdz
    dU = core.dU_vec
    M = [1.0 2*dzdr dzdr^2; dU[4] dU[5] 0.0; 0.0 dU[4] dU[5]]
    rhs = [0.0, -dU[2], -dU[3]]
    w = M \ rhs
    return [w[1] w[2]; w[2] w[3]]
end

"""
    integrate_ray_amplitude(prob, y0, W0, sigma0, sigma_end;
                            lnE20=0.0, phase0=0.0, damping=true,
                            wmax_x=3000.0, wmax_k=0.01, max_maslov=8,
                            rtol=1e-6, atol=1e-7,
                            initial_step=1e-7*(sigma_end-sigma0),
                            max_steps=200_000) -> AmplitudeTrace

Integrate one ray WITH amplitude transport: focusing tensor `W0` (2×2
symmetric, e.g. from [`antenna_focusing`](@ref) or `zeros(2,2)` for a locally
plane wavefront), log-amplitude² and eikonal phase, per-species collisionless
damping/deposition, and automatic Maslov switching to the k-space
representation when `‖W‖ > wmax_x` (back when `‖W̃‖ > wmax_k`), which carries
rays through caustics with the correct π/2 phase jumps.
"""
function integrate_ray_amplitude(
    prob::RayconProblem,
    y0::AbstractVector{<:Real},
    W0::AbstractMatrix{<:Real},
    sigma0::Real,
    sigma_end::Real;
    lnE20::Real = 0.0,
    phase0::Real = 0.0,
    damping::Bool = true,
    wmax_x::Real = 3000.0,
    wmax_k::Real = 0.01,
    max_maslov::Integer = 8,
    rtol::Real = 1e-6,
    atol::Real = 1e-7,
    initial_step::Real = 1e-7 * (Float64(sigma_end) - Float64(sigma0)),
    max_steps::Integer = 200_000,
)
    length(y0) == 4 || throw(ArgumentError("state must be (r, z, kr, kz)"))
    (size(W0) == (2, 2) && isapprox(W0[1, 2], W0[2, 1]; atol = 1e-10 * (1 + abs(W0[1, 2])))) ||
        throw(ArgumentError("W0 must be a symmetric 2×2 focusing tensor"))
    all(isfinite, y0) && all(isfinite, W0) || throw(ArgumentError("initial state must be finite"))
    σ0 = Float64(sigma0)
    σe = Float64(sigma_end)
    σe > σ0 || throw(ArgumentError("sigma_end must exceed sigma0"))
    (isfinite(rtol) && rtol > 0 && isfinite(atol) && atol > 0) ||
        throw(ArgumentError("rtol and atol must be positive"))
    (isfinite(wmax_x) && wmax_x > 0 && isfinite(wmax_k) && wmax_k > 0) ||
        throw(ArgumentError("caustic thresholds must be positive"))
    h = Float64(initial_step)
    (isfinite(h) && h > 0) || throw(ArgumentError("initial_step must be positive"))
    ns = length(prob.amass)
    nu = 9 + ns

    u = zeros(nu)
    u[1:4] .= Float64.(y0)
    u[5] = W0[1, 1]
    u[6] = 0.5 * (W0[1, 2] + W0[2, 1])
    u[7] = W0[2, 2]
    u[8] = Float64(lnE20)
    u[9] = Float64(phase0)
    inkspace = false
    nmaslov = 0

    σ = σ0
    sigmas = [σ]
    ys = reshape(copy(u[1:4]), 4, 1)
    Ws = reshape(copy(u[5:7]), 3, 1)
    lnE2s = [u[8]]
    phases = [u[9]]
    deps = reshape(copy(u[10:end]), ns, 1)
    inks = [inkspace]
    status = :max_steps
    hmin = 16 * eps(max(abs(σ0), abs(σe)))

    local k1
    try
        k1 = _amp_rhs(prob, u, inkspace, damping)
    catch e
        e isa DomainError &&
            return AmplitudeTrace(sigmas, ys, Ws, lnE2s, phases, deps, inks, 0, :left_domain)
        rethrow()
    end

    ks = Vector{Vector{Float64}}(undef, 7)
    for _ = 1:max_steps
        h = min(h, σe - σ)
        if h <= hmin
            status = σe - σ <= hmin ? :end_of_span : :step_underflow
            break
        end
        ks[1] = k1
        unew = u
        try
            for i = 1:6
                acc = copy(u)
                a = _DP_A[i]
                for j = 1:i
                    acc .+= (h * a[j]) .* ks[j]
                end
                ks[i+1] = _amp_rhs(prob, acc, inkspace, damping)
                i == 6 && (unew = acc)
            end
        catch e
            if e isa DomainError
                status = :left_domain
                break
            end
            rethrow()
        end
        errv = zeros(nu)
        for j = 1:7
            _DP_E[j] == 0 && continue
            errv .+= (h * _DP_E[j]) .* ks[j]
        end
        sc = atol .+ rtol .* max.(abs.(u), abs.(unew))
        err = sqrt(sum((errv ./ sc) .^ 2) / nu)
        if err > 1.0
            h *= max(0.2, 0.9 * err^(-0.2))
            continue
        end
        σ += h
        u = unew
        k1 = ks[7]
        h *= min(5.0, max(0.2, err > 0 ? 0.9 * err^(-0.2) : 5.0))

        # caustic handling: switch representation when the focusing tensor
        # leaves its healthy range (ray.m 'caustic_which'/'caustic_list';
        # upstream thresholds: ‖W‖ > 3000 in x-space, ‖W̃‖ > 0.01 in k-space —
        # a hysteresis pair, since the transform maps each far inside the
        # other's healthy range)
        wnorm = sqrt(u[5]^2 + 2 * u[6]^2 + u[7]^2)
        if (!inkspace && wnorm > wmax_x) || (inkspace && wnorm > wmax_k)
            if nmaslov >= max_maslov
                status = :caustic_failure
                break
            end
            inkspace, ok = _maslov!(u, inkspace)
            if !ok
                status = :caustic_failure
                break
            end
            nmaslov += 1
            k1 = _amp_rhs(prob, u, inkspace, damping)  # FSAL invalid after transform
        end

        push!(sigmas, σ)
        ys = hcat(ys, u[1:4])
        Ws = hcat(Ws, u[5:7])
        push!(lnE2s, u[8])
        push!(phases, u[9])
        deps = hcat(deps, u[10:end])
        push!(inks, inkspace)
        if σ >= σe
            status = :end_of_span
            break
        end
    end
    return AmplitudeTrace(sigmas, ys, Ws, lnE2s, phases, deps, inks, nmaslov, status)
end
