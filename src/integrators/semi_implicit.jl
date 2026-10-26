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
quadratic invariants to roundoff and is unconditionally stable. This is the
[1/1] Padé approximant of `exp(-z)`, with `z = im*ω*dt`, agreeing through `z²`
and having local error `O(z³)`.

Finite operands are combined without intermediate range loss. Floating results
follow Julia's promotion; integer-only inputs return `ComplexF64` (or
`Complex{BigFloat}` for `BigInt`). Rational-only results remain exact, widening
to `Rational{BigInt}` when the promoted integer type cannot hold the result.
Direct nonfinite operands retain Julia's original complex-expression behavior.
"""
function cn_multiplier(ω::Real, dt::Real)
    (!isfinite(ω) || !isfinite(dt)) && return (1 - im * ω * dt / 2) / (1 + im * ω * dt / 2)
    if (ω isa Rational || dt isa Rational) &&
       ω isa Union{Integer,Rational} &&
       dt isa Union{Integer,Rational}
        h = _exact_rational(ω) * _exact_rational(dt) / 2
        den = 1 + h^2
        re, imaginary = (1 - h^2) / den, -2h / den
        I = promote_type(Int, typeof(numerator(ω)), typeof(numerator(dt)))
        if I === BigInt || all(
            x -> typemin(I) <= x <= typemax(I),
            (numerator(re), denominator(re), numerator(imaginary), denominator(imaginary)),
        )
            return complex(Rational{I}(re), Rational{I}(imaginary))
        end
        return complex(re, imaginary)
    end
    T = typeof(float(zero(promote_type(typeof(ω), typeof(dt)))))
    if T <: Union{Float16,Float32,Float64}
        a, b = T(ω), T(dt)
        isfinite(a) && isfinite(b) && a == ω && b == dt && return cn_multiplier(a, b)
    end
    return _cn_complete(T, ω, dt)
end

function cn_multiplier(ω::BigFloat, dt::BigFloat)
    (!isfinite(ω) || !isfinite(dt)) && return (1 - im * ω * dt / 2) / (1 + im * ω * dt / 2)
    (iszero(ω) || iszero(dt)) && return _cn_complete(BigFloat, ω, dt)
    bits = Int128(precision(ω)) + precision(dt)
    power = Int128(exponent(ω)) + exponent(dt)
    low = min(Int128(2), 2power - 2bits + 4)
    high = max(Int128(2), 2power + 3) + 1
    minimum = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
    maximum = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
    # Exact intermediate widths depend on stored precision, never on a remote
    # absolute exponent. Compact ranges also need room for every intermediate.
    if abs(power) > bits || low < minimum - 1 || high + 1 > maximum
        return _cn_complete(BigFloat, ω, dt)
    end
    # x has at most `bits` bits and x² at most twice that. Aligning its lowest
    # possible bit with 4 gives exact numerators/denominator before division.
    x = BigFloat(; precision = Int(bits))
    squared = BigFloat(; precision = Int(2bits))
    numerator = BigFloat(; precision = Int(high - low + 1))
    denominator = BigFloat(; precision = Int(high - low + 1))
    re, imaginary = BigFloat(), BigFloat()
    mode = Base.MPFR.rounding_raw(BigFloat)
    ccall(
        (:mpfr_mul, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        x,
        ω,
        dt,
        mode,
    )
    ccall(
        (:mpfr_sqr, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        squared,
        x,
        mode,
    )
    ccall(
        (:mpfr_ui_sub, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Culong, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        numerator,
        4,
        squared,
        mode,
    )
    ccall(
        (:mpfr_add_ui, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Culong, Base.MPFR.MPFRRoundingMode),
        denominator,
        squared,
        4,
        mode,
    )
    ccall(
        (:mpfr_div, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        re,
        numerator,
        denominator,
        mode,
    )
    iszero(numerator) &&
        ccall((:mpfr_set_zero, Base.MPFR.libmpfr), Cvoid, (Ref{BigFloat}, Cint), re, 1)
    ccall(
        (:mpfr_mul_2si, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Clong, Base.MPFR.MPFRRoundingMode),
        x,
        x,
        2,
        mode,
    )
    ccall(
        (:mpfr_neg, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        x,
        x,
        mode,
    )
    ccall(
        (:mpfr_div, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        imaginary,
        x,
        denominator,
        mode,
    )
    return complex(re, imaginary)
end

@inline function cn_multiplier(ω::T, dt::T) where {T<:Union{Float16,Float32,Float64}}
    (!isfinite(ω) || !isfinite(dt)) && return (1 - im * ω * dt / 2) / (1 + im * ω * dt / 2)
    (iszero(ω) || iszero(dt)) && return complex(one(T), -copysign(zero(T), ω * dt))
    product = ω * dt
    magnitude = abs(product)
    if floatmin(T) * 2 <= magnitude <= sqrt(floatmax(T)) / 2
        if one(T) <= magnitude <= T(4)
            # The error-free product residual retains the real component near
            # quadrature, where rounding ω*dt to ±2 would erase it.
            W = T === Float16 ? Float32 : T
            a, b = W(ω), W(dt)
            p = a * b
            residual = fma(a, b, -p)
            denominator = 4 + p * p
            re = ((2 - p) - residual) * ((2 + p) + residual) / denominator
            imaginary = -4p / denominator
            return complex(T(re), T(imaginary))
        end
        return (1 - im * ω * dt / 2) / (1 + im * ω * dt / 2)
    end
    return _cn_complete(T, ω, dt)
end

# Keep exact significands separate from their absolute binary exponents.
# Base.decompose also covers integers and rationals; irrational constants have
# the same floating interpretation as the original arithmetic expression.
_cn_decompose(value::Real) = Base.decompose(value)
_cn_decompose(value::AbstractIrrational) = Base.decompose(float(value))
_cn_decompose(value::Real, ::Type{T}) where {T} = _cn_decompose(value)
_cn_decompose(value::AbstractIrrational, ::Type{T}) where {T} = Base.decompose(T(value))

function _cn_complete(::Type{T}, ω, dt) where {T<:AbstractFloat}
    negative = xor(signbit(ω), signbit(dt))
    (iszero(ω) || iszero(dt)) && return complex(one(T), negative ? zero(T) : -zero(T))
    an, ae, ad = _cn_decompose(ω, T)
    bn, be, bd = _cn_decompose(dt, T)
    n, d = abs(BigInt(an) * BigInt(bn)), abs(BigInt(ad) * BigInt(bd))
    power = Int128(ae) + Int128(be) - 1
    inverted = _compare_scaled(n, power, d, Int128(0)) > 0
    if inverted
        n, d, power = d, n, -power
    end
    # For 1/2 <= |t| <= 1 the exponent alignment is bounded by the stored
    # significand widths. This includes every cancellation near quadrature.
    if _compare_scaled(n, power, d, Int128(-1)) >= 0
        t = power >= 0 ? (n << Int(power)) // d : n // (d << Int(-power))
        denominator = 1 + t^2
        re = (1 - t^2) / denominator
        imaginary = 2t / denominator
        return complex(
            _from_exact(T, inverted ? -re : re),
            _from_exact(T, negative ? imaginary : -imaginary),
        )
    end
    # Here 0 < t < 1/2: the real magnitude lies in (3/5,1), so its
    # comparison scale stays near unity even when t has an enormous exponent.
    nn, dd = n^2, d^2
    real_compare = let nn = nn, dd = dd, power = power
        function (q, e)
            # Every queried threshold is in [0,1]. At one the exact result is
            # strictly smaller; this also avoids a negative numerator coefficient.
            _compare_scaled(q, e, BigInt(1), Int128(0)) >= 0 && return -1
            scale = BigInt(1) << Int(-e)
            return _compare_scaled(dd * (scale - q), Int128(0), nn * (scale + q), 2power)
        end
    end
    re = _round_exact_comparison(T, Int128(-1), real_compare, inverted)
    imaginary_compare = let nd = n * d, nn = nn, dd = dd, power = power
        (q, e) -> _compare_scaled_sum(nd, power + 1, dd * q, e, nn * q, 2power + e)
    end
    binade = Int128(ndigits(n; base = 2) - ndigits(d; base = 2)) + power
    _compare_scaled(n, power, d, binade) < 0 && (binade -= 1)
    imaginary_compare(BigInt(1), binade + 1) >= 0 && (binade += 1)
    imaginary = _round_exact_comparison(T, binade, imaginary_compare, !negative)
    return complex(re, imaginary)
end

"Forward-Euler per-step factor for `∂t b = −iω b` (|·| = √(1+ω²dt²) > 1)."
@inline euler_multiplier(ω::Real, dt::Real) = 1 - im * ω * dt

# whistler angular frequencies ω_k = c k² on an n-point periodic grid of length L.
# The evolved field b = b_y + i b_z is COMPLEX (no conjugate-symmetry constraint)
# and ω = c k² is even in k, so the Nyquist bin m = n/2 (even n) gets its full
# ω = c·(π n/L)² — zeroing it is the odd-derivative (ik) convention, which does
# not apply to an even operator (cf. the laplacian! Nyquist handling, §3).
_whistler_omega(n::Integer, L::Real, c::Real) = _whistler_omega(n, L, c, 1)

function _whistler_omega(n::Integer, L::Real, B0::Real, n0::Real)
    n = Int(n)
    ω = Vector{Float64}(undef, n)
    ordinary =
        B0 isa Union{Float16,Float32,Float64,Base.BitInteger} &&
        n0 isa Union{Float16,Float32,Float64,Base.BitInteger} &&
        L isa Union{Float16,Float32,Float64,Base.BitInteger,AbstractIrrational} &&
        Float64(B0) == B0 &&
        Float64(n0) == n0 &&
        (L isa AbstractIrrational || Float64(L) == L)
    # Frequencies are Float64 even when stored physical parameters are narrower.
    # Convert their exact hardware values before the coefficient division.
    c = ordinary ? Float64(B0) / Float64(n0) : 0.0
    @inbounds for m = 0:n-1
        mp = m <= n ÷ 2 ? m : m - n
        if iszero(mp) || iszero(B0)
            ω[m+1] = copysign(0.0, B0)
            continue
        end
        k = ordinary ? 2π * mp / L : 0.0
        squared = k^2
        value = c * squared
        if ordinary && all(x -> isfinite(x) && abs(x) >= floatmin(Float64), (c, k, squared, value))
            ω[m+1] = value
        else
            ratio, power = BigInt(1) // BigInt(1), Int128(0)
            for (factor, count) in ((B0, 1), (n0, -1), (Float64(2π), 2), (mp, 2), (L, -2))
                numerator, exponent, denominator = _cn_decompose(factor)
                ratio *= (BigInt(numerator) // BigInt(denominator))^count
                power += Int128(exponent) * count
            end
            ω[m+1] = _from_exact(Float64, ratio, power)
        end
        isfinite(ω[m+1]) ||
            throw(ArgumentError("whistler frequency is not representable as finite Float64"))
    end
    return ω
end

"""
    run_whistler(; method=:cn, n=128, L=2π, B0=1.0, n0=1.0, dt=0.05, nsteps=400,
                   seed=1) -> (; energy0, energyN, energy_ratio, bhat, ω)

Evolve a random transverse whistler packet on a periodic grid with the
`:cn` (Crank–Nicolson, conservative semi-implicit) or `:euler` (forward Euler,
explicit) integrator and report the field energy `Σ|b_k|²` before/after.
Frequencies combine `B0/n0 * (2π*m/L)^2` before exceptional rounding; a complete
frequency outside finite Float64 range raises `ArgumentError`.
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
    # Mode indices include negative Fourier modes; keep their arithmetic signed
    # when the caller supplies an unsigned grid count.
    n = Int(n)
    ω = _whistler_omega(n, L, B0, n0)
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
