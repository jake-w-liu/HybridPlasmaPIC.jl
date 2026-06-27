# Loading.jl — particle loaders (from particles.jl)

"""
    load_uniform!(ps, rng, lo::NTuple{D}, hi::NTuple{D})

Random uniform particle positions in the box `[lo, hi)`.
"""
function load_uniform!(ps::ParticleSet{D,T}, rng, lo::NTuple{D}, hi::NTuple{D}) where {D,T}
    @inbounds for d = 1:D
        L = T(hi[d] - lo[d])
        l = T(lo[d])
        xd = ps.x[d]
        for p in eachindex(xd)
            xd[p] = l + L * rand(rng, T)
        end
    end
    return ps
end

"""
    load_lattice!(ps, lo::NTuple{D}, hi::NTuple{D}, counts::NTuple{D,Int})

Place `prod(counts)` particles on a regular cell-centered D-dimensional lattice
filling `[lo, hi)` — the quiet/low-noise spatial load. Requires
`prod(counts) == N`.
"""
function load_lattice!(
    ps::ParticleSet{D,T},
    lo::NTuple{D},
    hi::NTuple{D},
    counts::NTuple{D,Int},
) where {D,T}
    prod(counts) == nparticles(ps) || throw(ArgumentError("N must equal prod(counts)"))
    @inbounds for (p, I) in enumerate(CartesianIndices(counts))
        t = Tuple(I)
        for d = 1:D
            ps.x[d][p] = T(lo[d]) + T(hi[d] - lo[d]) * (t[d] - T(0.5)) / counts[d]
        end
    end
    return ps
end

"""
    load_lattice_1d!(ps, lo, hi)

Evenly spaced 1-D positions (cell-centered) — the quiet/low-noise spatial load.
"""
load_lattice_1d!(ps::ParticleSet{1,T}, lo, hi) where {T} =
    load_lattice!(ps, (T(lo),), (T(hi),), (nparticles(ps),))

"""
    set_density_weight!(ps, n0, g)

Set every particle weight so a uniform load deposits number density `n0`
(w = n0·∏L / N). Without this the deposited density is `nppc/ΔV`, which silently
rescales every 1/n term in Ohm's law.
"""
function set_density_weight!(ps::ParticleSet{D,T}, n0, g::FourierGrid{D,T}) where {D,T}
    N = nparticles(ps)
    N > 0 || throw(ArgumentError("no particles"))
    fill!(ps.weight, T(n0) * prod(g.L) / N)
    return ps
end

"""
    load_maxwellian!(ps, rng, u0::NTuple{3}, vth::NTuple{3})

Drifting (bi-)Maxwellian velocities: vᶜ = u0ᶜ + vthᶜ·𝒩(0,1). Anisotropic `vth`
gives a bi-Maxwellian; equal `vth` an isotropic Maxwellian; `vth=0` a cold beam.
`vthᶜ = √(k_B Tᶜ/m)`.
"""
function load_maxwellian!(ps::ParticleSet{D,T}, rng, u0::NTuple{3}, vth::NTuple{3}) where {D,T}
    @inbounds for c = 1:3
        vc = ps.v[c]
        u = T(u0[c])
        s = T(vth[c])
        for p in eachindex(vc)
            vc[p] = u + s * randn(rng, T)
        end
    end
    return ps
end

"""
    load_quiet_velocities!(ps, rng, u0::NTuple{3}, vth::NTuple{3})

Quiet-start velocities using mirrored pairs: particle `i` and `i+N/2` get
`u0 ± g`, so the thermal momentum cancels exactly (Σ v = N·u0 to roundoff).
Requires an even particle count.
"""
function load_quiet_velocities!(
    ps::ParticleSet{D,T},
    rng,
    u0::NTuple{3},
    vth::NTuple{3},
) where {D,T}
    N = nparticles(ps)
    iseven(N) || throw(ArgumentError("quiet start needs an even particle count"))
    half = N ÷ 2
    @inbounds for c = 1:3
        vc = ps.v[c]
        u = T(u0[c])
        s = T(vth[c])
        for i = 1:half
            g = s * randn(rng, T)
            vc[i] = u + g
            vc[i+half] = u - g
        end
    end
    return ps
end

# ---------------------------------------------------------------- Boris mover

# One Boris velocity rotation given local E, B (Birdsall & Langdon). Pure
# rotation preserves |v| exactly; reduces to v += qm·E·dt over a full step when B=0.
