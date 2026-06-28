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
end

function Electrostatic1D(g::FourierGrid{1,T}, Nparticles::Integer; n0 = 1.0) where {T}
    n = g.n[1]
    Np = _particle_length(Nparticles)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    Electrostatic1D{T,typeof(g)}(g, n0T, zeros(T, n), zeros(T, n), zeros(T, Np))
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
    )
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

function _spectral_k2(g::FourierGrid{D,T}, I::CartesianIndex{D}) where {D,T}
    s = zero(T)
    idx = Tuple(I)
    @inbounds for d = 1:D
        s += g.kvec[d][idx[d]]^2
    end
    return s
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
    @inbounds for d = 1:D
        tb = g.tbuf
        for I in CartesianIndices(cb)
            idx = Tuple(I)
            k2 = _spectral_k2(g, I)
            kd = g.kvec[d][idx[d]]
            tb[I] = k2 == 0 ? zero(Complex{T}) : (-im * kd * cb[I]) / k2
        end
        g.iplan * tb
        for I in CartesianIndices(cb)
            es.E[d][I] = real(tb[I])
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
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    return es
end

function init_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}) where {D,T}
    _require_espic_electrons(e)
    density!(es.ne, e, es.g, CIC())
    poisson_E!(es)
    return es
end

"""
    step_espic!(es, electrons, dt)

Advance one leapfrog step: gather E → Boris kick (B=0) → drift → periodic wrap →
deposit density → solve Poisson.
"""
function step_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}, dt::Real) where {T}
    _require_espic_electrons(e)
    dtT = _validated_nonnegative_dt(T, dt; name = "step_espic!")
    g = es.g
    L = g.L[1]
    qm = -one(T)
    gather_scalar!(es.Ep, es.E, e, g, CIC())
    vx = e.v[1]
    xx = e.x[1]
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
            dtT,
        )
        vx[p] = nx
        xx[p] += dtT * nx
    end
    apply_periodic!(e, (zero(T),), (L,))
    density!(es.ne, e, g, CIC())
    poisson_E!(es)
    return es
end

function step_espic!(es::ElectrostaticPIC{D,T}, e::ParticleSet{D,T}, dt::Real) where {D,T}
    _require_espic_electrons(e)
    dtT = _validated_nonnegative_dt(T, dt; name = "step_espic!")
    g = es.g
    gather_vector!(es.Ep, es.E, e, g, CIC())
    qm = e.q / e.m
    vx, vy, vz = e.v
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
            dtT,
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
    return es
end

"Electrostatic field energy ∫ ½ E² dx."
field_energy(es::Electrostatic1D) = oftype(es.n0, 0.5) * sum(abs2, es.E) * prod(es.g.dx)

function field_energy(es::ElectrostaticPIC)
    s = sum(sum(abs2, es.E[c]) for c = 1:3)
    return oftype(es.n0, 0.5) * s * prod(es.g.dx)
end
