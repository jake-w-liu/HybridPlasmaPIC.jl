# oracles/hybrid_dispersion_oracle.jl — HYB-006 independent dispersion oracle.
#
# A standalone, high-precision linearized warm Hall-MHD eigenvalue solver for the
# hybrid (massless-electron, quasineutral) model. It is DELIBERATELY independent of
# the production spatial operator: it uses only LinearAlgebra on a small dense
# matrix, so it can cross-check the PIC dispersion without sharing any code path.
#
# Physics. Linearize the hybrid/Hall-MHD equations about a uniform state
# (ρ0, u0=0, B0, p0) with perturbations ∝ exp(i(k·x − ωt)) (∂t→−iω, ∇→ik):
#
#   momentum   ρ0 ∂t δu = δJ×B0 − ∇δp           (δJ = ∇×δB = ik×δB)
#   induction  ∂t δB = ∇×(δu×B0) − ∇×(δJ×B0/n0)  (+ ∇×(∇δp_e/n0) ≡ 0 for uniform n0)
#   pressure   ∂t δp = −ρ0 c_s² ∇·δu             (adiabatic, combined sound speed)
#
# In the hybrid model the ELECTRON pressure enters the ion momentum through the
# generalized-Ohm E (the Boris push feels E = −u_i×B + (J×B − ∇p_e)/n), so the
# linear ion-momentum pressure gradient is the TOTAL pressure with the combined
# sound speed  c_s² = (γ_i p_i0 + γ_e p_e0)/ρ0  (HYB-005). For uniform n0 the
# electron-pressure term in induction is ∇×(∇p_e/n0) = (1/n0)∇×∇p_e = 0, so only
# c_s and the Hall term survive — exactly as the production model behaves.
#
# Writing each equation as ω·(var) = (linear map of vars) gives a standard complex
# eigenproblem  M q = ω q,  q = (δu_x,δu_y,δu_z, δB_x,δB_y,δB_z, δp). Its eigenvalues
# are the dispersion ω(k) of ALL branches (Alfvén, fast/slow magnetosonic, whistler,
# ion-cyclotron) at arbitrary angle; eigenvectors give polarization; dω/dk gives the
# group velocity. Verified (in test_dispersion_oracle.jl) to reproduce every known
# closed-form branch.
#
# Normalization: Ω_ci = 1, d_i = 1, v_A = 1 ⇒ |B0| = v_A = 1 for ρ0 = 1. The Hall
# term carries the implicit 1/(Ω_ci d_i) = 1 in these units.

module HybridDispersionOracle

using LinearAlgebra

export hybrid_dispersion_matrix,
    dispersion_frequencies, dispersion_eigen, fast_slow_speeds, group_velocity

@inline _cross(a, b) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

"""
    hybrid_dispersion_matrix(k, B0; cs=0.0, n0=1.0, rho0=1.0) -> 7×7 ComplexF64

Linear operator `M` of the warm Hall-MHD eigenproblem `M q = ω q` for wavevector
`k` (3-tuple) and background field `B0` (3-tuple), with combined sound speed `cs`,
background density `n0`, mass density `rho0`. State order
`q = (δu_x,δu_y,δu_z, δB_x,δB_y,δB_z, δp)`.
"""
function hybrid_dispersion_matrix(
    k::NTuple{3,<:Real},
    B0::NTuple{3,<:Real};
    cs::Real = 0.0,
    n0::Real = 1.0,
    rho0::Real = 1.0,
)
    all(isfinite, k) || throw(ArgumentError("k must be finite"))
    all(isfinite, B0) || throw(ArgumentError("B0 must be finite"))
    (isfinite(cs) && cs >= 0) || throw(ArgumentError("cs must be finite and ≥ 0"))
    (isfinite(n0) && n0 > 0) || throw(ArgumentError("n0 must be finite and > 0"))
    (isfinite(rho0) && rho0 > 0) || throw(ArgumentError("rho0 must be finite and > 0"))
    kk = (float(k[1]), float(k[2]), float(k[3]))
    bb = (float(B0[1]), float(B0[2]), float(B0[3]))
    cs2 = float(cs)^2
    M = zeros(ComplexF64, 7, 7)
    @inbounds for j = 1:7
        du = (j == 1, j == 2, j == 3) .* 1.0
        dB = (j == 4, j == 5, j == 6) .* 1.0
        dp = j == 7 ? 1.0 : 0.0
        kxB = _cross(kk, dB)
        kxB_x_B0 = _cross(kxB, bb)
        # ω δu = [k δp − (k×δB)×B0]/ρ0
        M[1, j] = (kk[1] * dp - kxB_x_B0[1]) / rho0
        M[2, j] = (kk[2] * dp - kxB_x_B0[2]) / rho0
        M[3, j] = (kk[3] * dp - kxB_x_B0[3]) / rho0
        # ω δB = −k×(δu×B0) + i k×((k×δB)×B0)/n0
        uxB0 = _cross(du, bb)
        kx_uxB0 = _cross(kk, uxB0)
        kx_hall = _cross(kk, kxB_x_B0)
        M[4, j] = -kx_uxB0[1] + im * kx_hall[1] / n0
        M[5, j] = -kx_uxB0[2] + im * kx_hall[2] / n0
        M[6, j] = -kx_uxB0[3] + im * kx_hall[3] / n0
        # ω δp = ρ0 c_s² (k·δu)
        M[7, j] = rho0 * cs2 * _dot3(kk, du)
    end
    return M
end

"""
    dispersion_frequencies(k, B0; cs=0.0, n0=1.0, rho0=1.0; atol=1e-8) -> Vector{Float64}

Sorted POSITIVE real-frequency branches ω(k) (the physical, undamped roots).
"""
function dispersion_frequencies(
    k::NTuple{3,<:Real},
    B0::NTuple{3,<:Real};
    cs::Real = 0.0,
    n0::Real = 1.0,
    rho0::Real = 1.0,
    atol::Real = 1e-8,
)
    λ = eigvals(hybrid_dispersion_matrix(k, B0; cs, n0, rho0))
    return sort([
        real(z) for z in λ if real(z) > atol && abs(imag(z)) <= atol * max(1, abs(real(z)))
    ])
end

"""
    dispersion_eigen(k, B0; cs=0.0, n0=1.0, rho0=1.0) -> (freqs, modes)

All eigenfrequencies `freqs` (complex) and eigenvectors `modes` (columns) for the
linear operator, for polarization analysis. `modes[:,i]` is the
`(δu_x,δu_y,δu_z, δB_x,δB_y,δB_z, δp)` eigenvector of `freqs[i]`.
"""
function dispersion_eigen(
    k::NTuple{3,<:Real},
    B0::NTuple{3,<:Real};
    cs::Real = 0.0,
    n0::Real = 1.0,
    rho0::Real = 1.0,
)
    F = eigen(hybrid_dispersion_matrix(k, B0; cs, n0, rho0))
    return F.values, F.vectors
end

"""
    fast_slow_speeds(vA, cs, θ) -> (c_fast, c_slow)

Closed-form oblique MHD fast/slow magnetosonic phase speeds (HYB-005):
`c_{f,s}² = ½[(v_A²+c_s²) ± √((v_A²+c_s²)² − 4 v_A² c_s² cos²θ)]`, with `θ` the
angle between `k` and `B0`.
"""
function fast_slow_speeds(vA::Real, cs::Real, θ::Real)
    (isfinite(vA) && vA >= 0) || throw(ArgumentError("vA must be finite and ≥ 0"))
    (isfinite(cs) && cs >= 0) || throw(ArgumentError("cs must be finite and ≥ 0"))
    isfinite(θ) || throw(ArgumentError("θ must be finite"))
    s = vA^2 + cs^2
    disc = s^2 - 4 * vA^2 * cs^2 * cos(θ)^2
    disc = max(disc, 0.0)                     # guard tiny negative round-off
    root = sqrt(disc)
    cf = sqrt(0.5 * (s + root))
    cslow = sqrt(max(0.5 * (s - root), 0.0))
    return cf, cslow
end

"""
    group_velocity(branch, khat, kmag, B0; cs=0.0, n0=1.0, rho0=1.0, dk=1e-4) -> Float64

Numerical group velocity `dω/dk` of the `branch`-th positive branch (1 = lowest),
along the unit wavevector `khat`, evaluated at `|k| = kmag` by central difference.
"""
function group_velocity(
    branch::Integer,
    khat::NTuple{3,<:Real},
    kmag::Real,
    B0::NTuple{3,<:Real};
    cs::Real = 0.0,
    n0::Real = 1.0,
    rho0::Real = 1.0,
    dk::Real = 1e-4,
)
    nrm = sqrt(_dot3(khat, khat))
    nrm > 0 || throw(ArgumentError("khat must be nonzero"))
    u = (khat[1] / nrm, khat[2] / nrm, khat[3] / nrm)
    kp = kmag + dk
    km = kmag - dk
    fp = dispersion_frequencies((u[1] * kp, u[2] * kp, u[3] * kp), B0; cs, n0, rho0)
    fm = dispersion_frequencies((u[1] * km, u[2] * km, u[3] * km), B0; cs, n0, rho0)
    (branch >= 1 && branch <= length(fp) && branch <= length(fm)) ||
        throw(ArgumentError("branch $branch out of range for this k"))
    return (fp[branch] - fm[branch]) / (2 * dk)
end

end # module
