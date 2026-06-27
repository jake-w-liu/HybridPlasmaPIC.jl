# espic.jl — minimal 1D electrostatic particle-in-cell (full kinetic electrons +
# immobile ion background), the verifiable first target of the Phase-12 full-PIC
# extension. Kept SEPARATE from the hybrid loop. Electrons q=−1, m=1; ε0=1 ⇒
# ω_pe = 1 at n0 = 1. Spectral (FFT) Poisson; reuses ParticleSet, the Boris kick
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
    Electrostatic1D{T,typeof(g)}(g, T(n0), zeros(T, n), zeros(T, n), zeros(T, Nparticles))
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

"Deposit electron density and solve the field (call once after loading)."
function init_espic!(es::Electrostatic1D{T}, e::ParticleSet{1,T}) where {T}
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
    g = es.g
    L = g.L[1]
    qm = -one(T)
    dtT = T(dt)
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

"Electrostatic field energy ∫ ½ E² dx."
field_energy(es::Electrostatic1D) = oftype(es.n0, 0.5) * sum(abs2, es.E) * prod(es.g.dx)
