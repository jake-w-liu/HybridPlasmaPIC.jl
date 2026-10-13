# espic.jl — electrostatic particle-in-cell (full kinetic electrons + immobile
# ion background), the verifiable first target of the Phase-12 full-PIC extension.
# Kept SEPARATE from the hybrid loop. Electrons q=−1, m=1; ε0=1 ⇒ ω_pe = 1 at
# n0 = 1. Spectral (FFT) Poisson; reuses ParticleSet, the Boris kick
# (B=0 ⇒ pure E acceleration), and CIC deposit/gather.
#
# Verified oracles: cold Langmuir oscillation ω = ω_pe, and the cold two-stream
# instability growth rate γ = ω_pe/(2√2) at k v0 = √(3/8) ω_pe.
# ponytail: electrostatic is the cleanest exact oracle; the electromagnetic
# extension (Yee/PSTD fields + charge-conserving Esirkepov current, oracle
# ω²=ω_pe²+c²k²) is a further phase — the existing current! is NOT charge-
# conserving, so EM needs new deposition, not a reuse.

"""
    Electrostatic1D(g, Nparticles; n0=1.0)

State for a 1D electrostatic PIC on periodic grid `g`: electric field, electron
number density, a per-particle gather buffer, and the immobile-ion background
density `n0` (charge `ρ = n0 − n_e`).
"""
struct Electrostatic1D{T,G}
    g::G
    n0::T
    E::Vector{T}
    ne::Vector{T}
    Ep::Vector{T}
    primed::Base.RefValue{Bool}        # false until the one-time leapfrog velocity priming
    last_dt::Base.RefValue{T}          # preceding drift interval sets the velocity half level
end

function Electrostatic1D(g::FourierGrid{1,T}, Nparticles::Integer; n0 = 1.0) where {T}
    n = g.n[1]
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    Electrostatic1D{T,typeof(g)}(
        g,
        n0T,
        zeros(T, n),
        zeros(T, n),
        zeros(T, Np),
        Ref(false),
        Ref(zero(T)),
    )
end

"""
    ElectrostaticPIC(g, Nparticles; n0=1.0)

Dimension-parametric electrostatic PIC on a periodic `FourierGrid{D}`. The state
stores a 3-component electric field tuple over the `D`-dimensional mesh, electron
number density, per-particle gather buffers, and an immobile ion background
density `n0` (charge `ρ = n0 - n_e`). Components `D+1:3` of `E` are zero.
"""
struct ElectrostaticPIC{D,T,G}
    g::G
    n0::T
    E::NTuple{3,Array{T,D}}
    ne::Array{T,D}
    Ep::NTuple{3,Vector{T}}
    primed::Base.RefValue{Bool}        # false until the one-time leapfrog velocity priming
    last_dt::Base.RefValue{T}          # preceding drift interval sets the velocity half level
end

function ElectrostaticPIC(g::FourierGrid{D,T}, Nparticles::Integer; n0 = 1.0) where {D,T}
    1 <= D <= 3 || throw(ArgumentError("ElectrostaticPIC supports D = 1, 2, or 3"))
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    return ElectrostaticPIC{D,T,typeof(g)}(
        g,
        n0T,
        ntuple(_ -> zeros(T, g.n), 3),
        zeros(T, g.n),
        ntuple(_ -> zeros(T, Np), 3),
        Ref(false),
        Ref(zero(T)),
    )
end

function _resize_espic_particle_workspaces!(es::Electrostatic1D, N::Integer)
    length(es.Ep) == N || resize!(es.Ep, N)
    return es
end

function _resize_espic_particle_workspaces!(es::ElectrostaticPIC, N::Integer)
    for c = 1:3
        length(es.Ep[c]) == N || resize!(es.Ep[c], N)
    end
    return es
end

"Solve −∂²φ = ρ/ε0, E = −∂φ, spectrally: Ê_m = −i ρ̂_m/(ε0 k_m), Ê_0 = 0."
function poisson_E!(es::Electrostatic1D{T}) where {T}
    g = es.g
    n = g.n[1]
    cb = g.cbuf
    @inbounds for i = 1:n
        cb[i] = Complex{T}(es.n0 - es.ne[i])
    end
    g.plan * cb
    k = g.kvec[1]
    @inbounds for m = 1:n
        km = k[m]
        cb[m] = km == 0 ? zero(Complex{T}) : (-im * cb[m]) / km
    end
    g.iplan * cb
    @inbounds for i = 1:n
        es.E[i] = real(cb[i])
    end
    return es.E
end

# Longitudinal Fourier modes share exact-zero classification and one cold
# complete-expression projection. FFT input/output range remains owned by the
# caller: this path starts with the represented finite Fourier coefficients.
@inline function _spectral_wavevector(g::FourierGrid{D,T}, I::CartesianIndex{D}) where {D,T}
    return ntuple(d -> @inbounds(g.kvec[d][I[d]]), Val(D))
end

# A subnormal intermediate can lose bits that a later large factor restores.
# Keep the original operation order when every intermediate retains its range.
@inline _spectral_product_lost_range(a::T, b::T, p::T) where {T} =
    !isfinite(p) || (!iszero(a) && !iszero(b) && abs(p) < floatmin(T))
@inline _spectral_quotient_lost_range(a::T, b::T, p::T) where {T} =
    !isfinite(p) || (!iszero(a) && abs(p) < floatmin(T))
@inline _spectral_mode_is_null(k) = all(iszero, k)
function _spectral_norm2(k::NTuple{D,T}) where {D,T}
    total = zero(T)
    lost = false
    for kd in k
        square = kd^2
        lost |= _spectral_product_lost_range(kd, kd, square)
        total += square
    end
    return total, lost || !isfinite(total)
end
_spectral_exact_complex(z::Complex) = complex(_exact_rational(real(z)), _exact_rational(imag(z)))
function _spectral_exact_projection(
    a::NTuple{D,Complex{T}},
    k::NTuple{D,T},
    source::Complex{T},
    divisor::T,
) where {D,T}
    # The homogeneous zero mode needs no exact allocations, even at extreme k.
    iszero(source) && all(iszero, a) && return a
    kr = map(_exact_rational, k)
    ar = map(_spectral_exact_complex, a)
    k2 = sum(x -> x * x, kr)
    # Source/divisor and the original components remain unrounded until the
    # complete longitudinal correction is formed. In particular, dt is not
    # consumed by an overflowing/underflowing target before reaching here.
    defect =
        _spectral_exact_complex(source) / _exact_rational(divisor) -
        im * sum(kr[d] * ar[d] for d = 1:D)
    return map(ar, kr) do ad, kd
        result = ad - im * kd * defect / k2
        complex(_from_exact(T, real(result)), _from_exact(T, imag(result)))
    end
end
# Preserve callers' existing nonfinite-input behavior; the exact path accepts
# only finite represented operands and a nonzero divisor.
@inline _spectral_projection_eligible(a, k, source, divisor) =
    all(isfinite, a) &&
    all(isfinite, k) &&
    isfinite(source) &&
    isfinite(divisor) &&
    !iszero(divisor)
function _spectral_project_mode(
    a::NTuple{D,Complex{T}},
    k::NTuple{D,T},
    source::Complex{T},
    divisor::T,
) where {D,T}
    _spectral_mode_is_null(k) && return a
    k2, lost = _spectral_norm2(k)
    divergence = zero(Complex{T})
    for d = 1:D
        term = Complex{T}(0, k[d]) * a[d]
        lost |=
            _spectral_product_lost_range(k[d], imag(a[d]), real(term)) ||
            _spectral_product_lost_range(k[d], real(a[d]), imag(term))
        divergence += term
        lost |= !isfinite(divergence)
    end
    target = source / divisor
    lost |=
        _spectral_quotient_lost_range(real(source), divisor, real(target)) ||
        _spectral_quotient_lost_range(imag(source), divisor, imag(target))
    defect = target - divergence
    lost |= !isfinite(defect)
    coefficients = map(kd -> (-im * kd / k2), k)
    for d = 1:D
        lost |= _spectral_quotient_lost_range(k[d], k2, imag(coefficients[d]))
    end
    corrections = map(coefficient -> coefficient * defect, coefficients)
    for d = 1:D
        lost |=
            _spectral_product_lost_range(
                imag(coefficients[d]),
                imag(defect),
                real(corrections[d]),
            ) ||
            _spectral_product_lost_range(imag(coefficients[d]), real(defect), imag(corrections[d]))
    end
    result = map(+, a, corrections)
    lost |= !all(isfinite, result)
    cold = lost && _spectral_projection_eligible(a, k, source, divisor)
    return cold ? _spectral_exact_projection(a, k, source, divisor) : result
end
function _spectral_poisson_mode(k::NTuple{D,T}, source::Complex{T}) where {D,T}
    zeros_mode = map(_ -> zero(Complex{T}), k)
    _spectral_mode_is_null(k) && return zeros_mode
    k2, lost = _spectral_norm2(k)
    numerator = map(kd -> (-im * kd * source), k)
    result = map(value -> value / k2, numerator)
    for d = 1:D
        lost |=
            _spectral_product_lost_range(k[d], imag(source), real(numerator[d])) ||
            _spectral_product_lost_range(k[d], real(source), imag(numerator[d]))
        lost |=
            _spectral_quotient_lost_range(real(numerator[d]), k2, real(result[d])) ||
            _spectral_quotient_lost_range(imag(numerator[d]), k2, imag(result[d]))
    end
    cold = lost && _spectral_projection_eligible(zeros_mode, k, source, one(T))
    return cold ? _spectral_exact_projection(zeros_mode, k, source, one(T)) : result
end

"""
    poisson_E!(es::ElectrostaticPIC)

Solve `-∇²φ = ρ`, `E = -∇φ` spectrally in 1D, 2D, or 3D.
"""
function poisson_E!(es::ElectrostaticPIC{D,T}) where {D,T}
    g = es.g
    cb = g.cbuf
    @inbounds for I in CartesianIndices(cb)
        cb[I] = Complex{T}(es.n0 - es.ne[I])
    end
    g.plan * cb
    for c = 1:3
        fill!(es.E[c], zero(T))
    end
    buffers = (g.tbuf, g.abuf, g.sbuf)
    @inbounds for I in CartesianIndices(cb)
        values = _spectral_poisson_mode(_spectral_wavevector(g, I), cb[I])
        for d = 1:D
            buffers[d][I] = values[d]
        end
    end
    for d = 1:D
        g.iplan * buffers[d]
        for I in CartesianIndices(cb)
            es.E[d][I] = real(buffers[d][I])
        end
    end
    return es.E
end

@inline function _require_espic_electrons(e::ParticleSet{D,T}) where {D,T}
    q = _require_finite_real("electron charge q", e.q, T)
    m = _require_finite_real("electron mass m", e.m, T)
    q == -one(T) ||
        throw(ArgumentError("electrostatic PIC requires electron ParticleSet with q = -1"))
    m == one(T) ||
        throw(ArgumentError("electrostatic PIC requires electron ParticleSet with m = 1"))
    return nothing
end

"Deposit electron density and solve the field (call once after loading)."
function init_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}) where {T}
    _require_espic_electrons(e)
    _resize_espic_particle_workspaces!(es, nparticles(e))
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    es.primed[] = false                 # re-arm the one-time priming for the next step
    es.last_dt[] = zero(T)
    return es
end

function init_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}) where {D,T}
    _require_espic_electrons(e)
    _resize_espic_particle_workspaces!(es, nparticles(e))
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    es.primed[] = false                 # re-arm the one-time priming for the next step
    es.last_dt[] = zero(T)
    return es
end

"""
    step_espic!(es, electrons, dt)

Advance one leapfrog step: gather E → Boris kick (B=0) → drift → periodic wrap →
deposit density → solve Poisson. When `dt` changes, the velocity kick spans
`(last_dt + dt)/2` between the preceding and next half levels; positions advance
by `dt`. A zero timestep leaves both the state and its time staggering unchanged.
"""
function step_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}, dt::Real) where {T}
    _require_espic_electrons(e)
    dtT = _validated_nonnegative_dt(T, dt; name = "step_espic!")
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _resize_espic_particle_workspaces!(es, nparticles(e))
    g = es.g
    L = g.L[1]
    qm = -one(T)
    gather_scalar!(es.Ep, es.E, e, g, CIC())
    vx = e.v[1]
    xx = e.x[1]
    kick_dt = es.primed[] ? es.last_dt[] / 2 + dtT / 2 : dtT
    if !es.primed[]                     # prime once: loaded v is physical v^0 → v^{-1/2} (2nd order)
        @inbounds for p in eachindex(e.weight)
            vx[p] -= (dtT / 2) * qm * es.Ep[p]      # a^0 = qm E (B=0)
        end
        es.primed[] = true
    end
    @inbounds for p in eachindex(e.weight)
        nx, _, _ = boris_kick(
            vx[p],
            zero(T),
            zero(T),
            es.Ep[p],
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            qm,
            kick_dt,
        )
        vx[p] = nx
        xx[p] += dtT * nx
    end
    apply_periodic!(e, (zero(T),), (L,))
    density!(es.ne, e, g, CIC())
    poisson_E!(es)
    es.last_dt[] = dtT
    return es
end

function step_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}, dt::Real) where {D,T}
    _require_espic_electrons(e)
    dtT = _validated_nonnegative_dt(T, dt; name = "step_espic!")
    iszero(dtT) && return es            # dt=0 no-op: do not consume the one-time priming
    _resize_espic_particle_workspaces!(es, nparticles(e))
    g = es.g
    gather_vector!(es.Ep, es.E, e, g, CIC())
    qm = e.q / e.m
    vx, vy, vz = e.v
    kick_dt = es.primed[] ? es.last_dt[] / 2 + dtT / 2 : dtT
    if !es.primed[]                     # prime once: loaded v is physical v^0 → v^{-1/2} (2nd order)
        h = dtT / 2
        @inbounds for p in eachindex(e.weight)
            vx[p] -= h * qm * es.Ep[1][p]           # a^0 = qm E (B=0)
            vy[p] -= h * qm * es.Ep[2][p]
            vz[p] -= h * qm * es.Ep[3][p]
        end
        es.primed[] = true
    end
    @inbounds for p in eachindex(e.weight)
        nx, ny, nz = boris_kick(
            vx[p],
            vy[p],
            vz[p],
            es.Ep[1][p],
            es.Ep[2][p],
            es.Ep[3][p],
            zero(T),
            zero(T),
            zero(T),
            qm,
            kick_dt,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            e.x[d][p] += dtT * e.v[d][p]
        end
    end
    apply_periodic!(e, ntuple(_ -> zero(T), D), g.L)
    density!(es.ne, e, g, CIC())
    poisson_E!(es)
    es.last_dt[] = dtT
    return es
end

"Electrostatic field energy ∫ ½ E² dx."
field_energy(es::Electrostatic1D) = oftype(es.n0, 0.5) * sum(abs2, es.E) * prod(es.g.dx)

function field_energy(es::ElectrostaticPIC)
    s = sum(sum(abs2, es.E[c]) for c = 1:3)
    return oftype(es.n0, 0.5) * s * prod(es.g.dx)
end
