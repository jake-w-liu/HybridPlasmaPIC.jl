# completeness.jl — closure/consistency diagnostics that exercise the full
# coupling chain rather than a single operator:
#
#   * particle_work!     — the energy the field does on the particle population,
#     ΔW_p += q·(v·E_gathered)·dt, the per-particle work increment that closes
#     the field↔particle energy budget (gather is the transpose of deposit).
#   * mixed_divcurl_residual — the discrete div(curl A) on the mixed
#     SBP-x × Fourier-y shock grid. Analytically zero; discretely bounded by the
#     SBP-x truncation error because ∂x (SBP) and ∂y (Fourier) do NOT commute
#     exactly (Fourier-y is spectrally exact, so the residual is entirely the
#     mixed second derivative ∂x∂y − ∂y∂x acting through the one-sided SBP
#     boundary closure). It is small and CONVERGES under x-refinement at the
#     SBP-(2,1) interior rate.

"""
    particle_work!(work, ps, Egrid::NTuple{3}, g, shape, dt) -> work

Accumulate, per particle, the work the electric field does over one step:

    work[p] += q · (v_p · E(x_p)) · dt

`E(x_p)` is the grid field `Egrid = (Ex, Ey, Ez)` gathered to the particle
position with `shape` (the transpose of [`deposit_scalar!`](@ref)), `v_p` is the
particle's full 3-velocity, `q = ps.q` is the species charge, and `dt` the step.
`work` is a length-`nparticles(ps)` vector that is *added to* (not zeroed), so it
can be summed across a run. Returns `work`.

For a spatially uniform field and constant velocity this reduces exactly to
`q (v·E) dt` per particle, and the population total to `q (v·E) dt · Σ_p` over the
sampled particles (a roundoff-level identity, see the test).
"""
function particle_work!(
    work::AbstractVector{T},
    ps::ParticleSet{D,T},
    Egrid::NTuple{3,<:Array{T,D}},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    dt::Real,
) where {D,T}
    Np = nparticles(ps)
    length(work) == Np ||
        throw(DimensionMismatch("work length $(length(work)) ≠ particle count $Np"))
    dtT = _require_finite_real("dt", dt, T)
    q = ps.q
    # gather the three field components to each particle (transpose of deposit)
    Ep = ntuple(_ -> Vector{T}(undef, Np), 3)
    for c = 1:3
        gather_scalar!(Ep[c], Egrid[c], ps, g, shape)
    end
    vx, vy, vz = ps.v
    ex, ey, ez = Ep
    @inbounds for p = 1:Np
        vdotE = vx[p] * ex[p] + vy[p] * ey[p] + vz[p] * ez[p]
        work[p] += q * vdotE * dtT
    end
    return work
end

"""
    mixed_divcurl_residual(sbp::SBP1D, nx, ny, Ly) -> r

Measure how well the mixed shock operator (non-periodic SBP along the shock
normal x with `nx` nodes over `[0, sbp.L]`, periodic Fourier along the
transverse y with `ny` nodes over `[0, Ly)`) preserves the solenoidal identity
`∇·(∇×A) = 0`. Builds a smooth manufactured 3-vector potential `A` (each
component an x-profile times a band-limited y-profile), takes `B = ∇×A` in
**closed form** (the manufactured `A` has an exact analytic curl), samples that
exact solenoidal `B` on the grid, applies the **mixed** divergence (`∂x` via
[`sbp_deriv_x!`](@ref), `∂y` via [`fourier_deriv_y!`](@ref)), and returns the
H-norm-normalized residual

    r = ‖∇·B‖_H / ‖B‖_H ,   B = ∇×A  (analytically ∇·B ≡ 0).

Expected magnitude and convergence (verified in `test_completeness.jl`):

  * `r` is **small but NOT machine zero** — it is the SBP-x truncation of `∂x Bx`
    (Fourier-y resolves `∂y By` to roundoff for a band-limited y-profile, so the
    residual is entirely the non-periodic SBP boundary/interior truncation).
  * `r` **CONVERGES as `nx` grows** at the diagonal-norm SBP-(2,1) H-norm rate
    (≈ 3/2; the one-sided 1st-order boundary closure caps the global 2nd-order
    interior). For `Lx = 1`, `ny = 16` it runs ~3.5e-3 (nx=33) → ~1.6e-4
    (nx=257).

This is the honest measure of the identity on the mixed grid: if `B` is instead
recomputed from `A` with the *same* mixed operators, `∇·(∇×A)` collapses to
machine zero because SBP-x and Fourier-y act on independent tensor axes and
therefore commute exactly — that cancellation hides the SBP truncation, so the
analytic curl is used to expose it. `sbp.dx` and `sbp.n` fix the x extent
`Lx = sbp.dx·(sbp.n−1)`; the operator is rebuilt at the requested `nx`.
"""
function mixed_divcurl_residual(sbp::SBP1D{T}, nx::Integer, ny::Integer, Ly::Real) where {T}
    nx >= 3 || throw(ArgumentError("need nx ≥ 3"))
    ny >= 2 || throw(ArgumentError("need ny ≥ 2"))
    Lx = sbp.dx * (sbp.n - 1)               # physical x extent of the supplied operator
    s = SBP1D(Int(nx), T(Lx))               # SBP operator at the requested x resolution
    LyT = _require_finite_positive_real("Ly", Ly, T)
    x = collect(range(zero(T), T(Lx); length = Int(nx)))
    y = T[(j - 1) * LyT / ny for j = 1:Int(ny)]

    # Smooth, non-periodic-in-x manufactured potential A and its exact curl.
    # Each component is an x-profile (which the SBP operator differentiates only
    # approximately) times a band-limited y-profile (spectrally exact in y).
    az(xi) = T(0.7) * sin(T(1.5) * xi) + T(0.2) * xi^2
    daz(xi) = T(1.05) * cos(T(1.5) * xi) + T(0.4) * xi      # az'(x)
    ky = T(2)                               # integer-mode-resolved on the y grid
    fy(yj) = cos(ky * yj) + T(0.4) * sin(T(2) * ky * yj)
    dfy(yj) = -ky * sin(ky * yj) + T(0.8) * ky * cos(T(2) * ky * yj)   # fy'(y)

    # B = ∇×A with ∂z = 0, only A_z = az(x)·fy(y) contributing to the in-plane B:
    #   Bx =  ∂y A_z =  az(x)·fy'(y)
    #   By = −∂x A_z = −az'(x)·fy(y)
    # ⇒ ∇·B = ∂x Bx + ∂y By = az'(x)fy'(y) − az'(x)fy'(y) ≡ 0 analytically.
    Bx = T[az(x[i]) * dfy(y[j]) for i = 1:Int(nx), j = 1:Int(ny)]
    By = T[-daz(x[i]) * fy(y[j]) for i = 1:Int(nx), j = 1:Int(ny)]

    dBx_dx = similar(Bx)
    dBy_dy = similar(By)
    sbp_deriv_x!(dBx_dx, Bx, s)             # SBP-x: truncation error survives
    ywork = FourierDerivYWorkspace(By, LyT)
    fourier_deriv_y!(dBy_dy, By, ywork)     # Fourier-y: roundoff-exact
    divB = dBx_dx .+ dBy_dy

    # H-norm over x (diagonal SBP quadrature), uniform sum over periodic y.
    H = s.H
    num = zero(T)
    den = zero(T)
    @inbounds for j = 1:Int(ny), i = 1:Int(nx)
        num += H[i] * divB[i, j]^2
        den += H[i] * (Bx[i, j]^2 + By[i, j]^2)
    end
    den > 0 || throw(ArgumentError("degenerate manufactured field"))
    return sqrt(num) / sqrt(den)
end
