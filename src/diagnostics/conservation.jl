# budgets.jl — energy & momentum budget diagnostics (§6/§29).
#
# Thin, composable wrappers over the existing conserved-total diagnostics
# (kinetic_energy, magnetic_energy, electron_internal_energy, total_momentum,
# electric_work) plus the local electromagnetic-work and Ohmic-heating rates.
# Everything is in Ω_ci-normalized units; volume integrals use prod(g.dx).

# ---------------------------------------------------------------- energy budget

"""
    energy_budget(ps, B, n, closure, g)
        -> (kinetic, magnetic, electron_internal, total)

Total conserved energies of a hybrid state: ion kinetic `Σ_p ½ m w |v|²`
([`kinetic_energy`](@ref)), magnetic `∫ ½|B|² dV` ([`magnetic_energy`](@ref)),
and electron internal `∫ p_e/(γ−1) dV` ([`electron_internal_energy`](@ref)).
`total` is their sum.

For an isothermal closure (γ=1) the electron internal energy has no closed
invariant, so `electron_internal` (and hence `total`) is `NaN` — use a
polytropic closure (γ≠1) for an energy invariant.

`B` is the 3-tuple of magnetic-field arrays (e.g. `stepper.fields.B`), `n` the
electron/ion density array, `closure` the `ElectronClosure`, and `g` the
`FourierGrid`.
"""
function energy_budget(
    ps::ParticleSet{D,T},
    B::NTuple{3,<:Array{T,D}},
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_tuple(:B, B, g)
    _require_grid_array(:n, n, g)
    ek = kinetic_energy(ps)
    em = magnetic_energy(B, g)
    ei = electron_internal_energy(n, closure, g)
    return (kinetic = ek, magnetic = em, electron_internal = ei, total = ek + em + ei)
end

# ---------------------------------------------------------------- momentum budget

"""
    momentum_budget(ps, B, g) -> (particle::NTuple{3}, total::NTuple{3})

Momentum budget of a hybrid state. `particle` is the total ion momentum
`Σ_p m w v` ([`total_momentum`](@ref)). `total` adds the electromagnetic field
momentum `∫ (E×B) dV`.

In the hybrid (Darwin / quasi-neutral) approximation there is no displacement
current and the electron mass is neglected, so the electromagnetic field
momentum is not a dynamical reservoir; it is taken to be zero here and
`total ≡ particle`. The `B` and `g` arguments are accepted for interface
symmetry and future extension (e.g. an explicit-E momentum reservoir).
"""
function momentum_budget(
    ps::ParticleSet{D,T},
    B::NTuple{3,<:Array{T,D}},
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_tuple(:B, B, g)
    p = total_momentum(ps)                 # (px, py, pz)
    # Field momentum omitted (no displacement current in the hybrid model).
    field = (zero(T), zero(T), zero(T))
    total = (p[1] + field[1], p[2] + field[2], p[3] + field[3])
    return (particle = p, total = total)
end

# ---------------------------------------------------------------- J·E work rate

"""
    jdotE_density(J, E) -> Array

Per-cell electromagnetic work rate `J·E = Σ_c J_c E_c` (the local rate of
energy transfer from the field to the particles). `J` and `E` are 3-tuples of
grid arrays of identical shape. Integrating the result against `prod(g.dx)`
reproduces [`electric_work`](@ref)`(J, E, g)`.
"""
function jdotE_density(J::NTuple{3,<:AbstractArray{T}}, E::NTuple{3,<:AbstractArray{T}}) where {T}
    _require_same_axes(:JdotE_input, (J..., E...))
    out = similar(J[1])
    @inbounds for I in eachindex(out, J[1], J[2], J[3], E[1], E[2], E[3])
        out[I] = J[1][I] * E[1][I] + J[2][I] * E[2][I] + J[3][I] * E[3][I]
    end
    return out
end

# ---------------------------------------------------------------- resistive heating

"""
    resistive_dissipation(J, η, g) -> Real

Volume-integrated Ohmic heating rate `∫ η |J|² dV = η Σ_g |J_g|² · prod(g.dx)`,
the rate at which resistivity converts field energy into electron internal
energy. `J` is the 3-tuple of current-density arrays, `η` the (scalar)
resistivity, `g` the `FourierGrid`.
"""
function resistive_dissipation(
    J::NTuple{3,<:AbstractArray{T}},
    η::Real,
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_tuple(:J, J, g)
    ηT = _require_finite_nonnegative_real("η", η, T)
    s = zero(T)
    @inbounds for I in eachindex(J[1], J[2], J[3])
        s += J[1][I] * J[1][I] + J[2][I] * J[2][I] + J[3][I] * J[3][I]
    end
    return ηT * s * prod(g.dx)
end

# --- particle momentum / electric work (from diagnostics.jl) ---
function total_momentum(ps::ParticleSet{D,T}) where {D,T}
    px = zero(T)
    py = zero(T)
    pz = zero(T)
    vx, vy, vz = ps.v
    w = ps.weight
    @inbounds for p in eachindex(w)
        px += w[p] * vx[p]
        py += w[p] * vy[p]
        pz += w[p] * vz[p]
    end
    m = ps.m
    return (m * px, m * py, m * pz)
end

"Total momentum summed over several species."
function total_momentum(species::AbstractVector{<:ParticleSet})
    p = total_momentum(first(species))
    for i = 2:length(species)
        q = total_momentum(species[i])
        p = (p[1] + q[1], p[2] + q[2], p[3] + q[3])
    end
    return p
end

# ---------------------------------------------------------------- electromagnetic work

"Volume-integrated electric work ∫ J·E dV (J, E are 3-tuples of grid arrays)."
function electric_work(
    J::NTuple{3,<:AbstractArray{T}},
    E::NTuple{3,<:AbstractArray{T}},
    g::FourierGrid,
) where {T}
    _require_grid_tuple(:J, J, g)
    _require_grid_tuple(:E, E, g)
    s = zero(T)
    @inbounds for I in eachindex(J[1], J[2], J[3], E[1], E[2], E[3])
        s += J[1][I] * E[1][I] + J[2][I] * E[2][I] + J[3][I] * E[3][I]
    end
    return s * prod(g.dx)
end

# ---------------------------------------------------------------- temperatures
