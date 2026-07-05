# semiimplicit.jl — a conservative semi-implicit (Crank–Nicolson) field
# integrator and an explicit/implicit integrator comparison (§ integrator notes,
# items "compare integrators" + "conservative semi-implicit later extension").
#
# Target stiff mode: the parallel WHISTLER. Linearizing electron-MHD about a
# uniform guide field B0 x̂ with density n0, the transverse perturbation
# b = b_y + i b_z (parallel propagation, ∂x only) obeys
#   ∂t b = +i c ∂xx b,     c = B0 / n0   (Ω_ci-normalized; d_i = 1),
# i.e. per Fourier mode  ∂t b_k = −i ω_k b_k  with the whistler ω_k = c k²
# (ω ∝ k² — the dispersive branch that forces tiny explicit time steps).
#
# A Crank–Nicolson (trapezoidal, θ=½) update of ∂t b = −iω b is the per-mode
# multiplier
#   b^{n+1}_k = M_k b^n_k,   M_k = (1 − iω_k dt/2)/(1 + iω_k dt/2),   |M_k| ≡ 1.
# Because |M_k| = 1 for EVERY ω_k and dt, CN conserves the field energy
# Σ|b_k|² to round-off and is UNCONDITIONALLY stable (no whistler CFL). Forward
# Euler gives b^{n+1}_k = (1 − iω_k dt) b^n_k with |·| = √(1+ω²dt²) > 1 — it
# always grows. This is the conservative semi-implicit method; it is verified
# against the explicit scheme in test_semiimplicit.jl.

"""
    cn_multiplier(ω, dt) -> Complex

Crank–Nicolson (trapezoidal) per-step amplification factor for `∂t b = −iω b`:
`(1 − iω dt/2)/(1 + iω dt/2)`. Its modulus is exactly 1, so the scheme conserves
quadratic invariants and is unconditionally stable; the phase is the rational
(2,2)-Padé approximant of `e^{−iω dt}` (error `O((ω dt)³)` per step).
"""
@inline cn_multiplier(ω::Real, dt::Real) = (1 - im * ω * dt / 2) / (1 + im * ω * dt / 2)

"Forward-Euler per-step factor for `∂t b = −iω b` (|·| = √(1+ω²dt²) > 1)."
@inline euler_multiplier(ω::Real, dt::Real) = 1 - im * ω * dt

# whistler angular frequencies ω_k = c k² on an n-point periodic grid of length L.
# The evolved field b = b_y + i b_z is COMPLEX (no conjugate-symmetry constraint)
# and ω = c k² is even in k, so the Nyquist bin m = n/2 (even n) gets its full
# ω = c·(π n/L)² — zeroing it is the odd-derivative (ik) convention, which does
# not apply to an even operator (cf. the laplacian! Nyquist handling, §3).
function _whistler_omega(n::Integer, L::Real, c::Real)
    ω = Vector{Float64}(undef, n)
    @inbounds for m = 0:n-1
        mp = m <= n ÷ 2 ? m : m - n
        k = 2π * mp / L
        ω[m+1] = c * k^2
    end
    return ω
end

"""
    run_whistler(; method=:cn, n=128, L=2π, B0=1.0, n0=1.0, dt=0.05, nsteps=400,
                   seed=1) -> (; energy0, energyN, energy_ratio, bhat, ω)

Evolve a random transverse whistler packet on a periodic grid with the
`:cn` (Crank–Nicolson, conservative semi-implicit) or `:euler` (forward Euler,
explicit) integrator and report the field energy `Σ|b_k|²` before/after.
"""
function run_whistler(;
    method::Symbol = :cn,
    n::Integer = 128,
    L::Real = 2π,
    B0::Real = 1.0,
    n0::Real = 1.0,
    dt::Real = 0.05,
    nsteps::Integer = 400,
    seed::Integer = 1,
    band::Union{Nothing,Integer} = nothing,
)
    n >= 1 || throw(ArgumentError("n must be positive"))
    isfinite(L) && L > 0 || throw(ArgumentError("L must be finite and positive"))
    isfinite(B0) || throw(ArgumentError("B0 must be finite"))
    isfinite(n0) && n0 > 0 || throw(ArgumentError("n0 must be finite and positive"))
    isfinite(dt) || throw(ArgumentError("dt must be finite"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    band === nothing || band >= 0 || throw(ArgumentError("band must be non-negative"))
    ω = _whistler_omega(n, L, B0 / n0)
    rng = MersenneTwister(Int(seed))
    bhat = [complex(randn(rng), randn(rng)) for _ = 1:n]
    # optional band-limit: keep only Fourier modes with |mode index| ≤ band
    # (the low-k, well-resolved part of the spectrum) so an explicit integrator
    # is itself stable and a CN-vs-explicit agreement comparison is meaningful.
    if band !== nothing
        @inbounds for m = 1:n
            mp = m - 1 <= n ÷ 2 ? m - 1 : m - 1 - n
            abs(mp) > band && (bhat[m] = 0)
        end
    end
    energy0 = sum(abs2, bhat)
    step! =
        method === :cn ? (m -> cn_multiplier(ω[m], dt)) :
        method === :euler ? (m -> euler_multiplier(ω[m], dt)) :
        throw(ArgumentError("method must be :cn or :euler"))
    @inbounds for _ = 1:nsteps, m = 1:n
        bhat[m] *= step!(m)
    end
    energyN = sum(abs2, bhat)
    return (; energy0, energyN, energy_ratio = energyN / energy0, bhat, ω)
end

"""
    compare_integrators_whistler(; dt_resolved=0.02, dt_stiff=0.5, kwargs...)
        -> (; agree_resolved, cn_ratio_stiff, euler_ratio_stiff)

Compare the conservative semi-implicit (CN) and explicit (Euler) integrators on
an identical whistler packet:

  • `agree_resolved` — max field difference between CN and Euler at a small,
    well-resolved `dt_resolved` (the two integrators agree where both are
    accurate);
  • `cn_ratio_stiff` — CN energy ratio at a large `dt_stiff` (≈1: stable +
    conservative);
  • `euler_ratio_stiff` — Euler energy ratio at the same `dt_stiff` (≫1:
    the explicit scheme is unstable on the stiff whistler).
"""
function compare_integrators_whistler(;
    dt_resolved::Real = 0.01,
    dt_stiff::Real = 0.5,
    nsteps::Integer = 200,
    band_resolved::Integer = 1,
    kwargs...,
)
    isfinite(dt_resolved) || throw(ArgumentError("dt_resolved must be finite"))
    isfinite(dt_stiff) || throw(ArgumentError("dt_stiff must be finite"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    band_resolved >= 0 || throw(ArgumentError("band_resolved must be non-negative"))
    # Agreement: on a band-limited (low-k, well-resolved) mode, CN and the
    # explicit integrator are both accurate and agree.
    rc = run_whistler(; method = :cn, dt = dt_resolved, nsteps, band = band_resolved, kwargs...)
    re = run_whistler(; method = :euler, dt = dt_resolved, nsteps, band = band_resolved, kwargs...)
    agree_resolved = maximum(abs, rc.bhat .- re.bhat) / sqrt(rc.energy0)
    # Stiffness: on the FULL whistler spectrum at a large dt, CN conserves energy
    # (ratio ≈ 1, unconditionally stable) while explicit Euler blows up.
    sc = run_whistler(; method = :cn, dt = dt_stiff, nsteps, kwargs...)
    se = run_whistler(; method = :euler, dt = dt_stiff, nsteps, kwargs...)
    return (; agree_resolved, cn_ratio_stiff = sc.energy_ratio, euler_ratio_stiff = se.energy_ratio)
end
