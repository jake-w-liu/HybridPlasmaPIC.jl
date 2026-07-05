# budgets.jl — energy & momentum budget diagnostics (§6/§29).
#
# Thin, composable wrappers over the existing conserved-total diagnostics
# (kinetic_energy, magnetic_energy, electron_internal_energy, total_momentum,
# electric_work) plus the local electromagnetic-work and Ohmic-heating rates.
# Everything is in Ω_ci-normalized units; volume integrals use prod(g.dx).

# ---------------------------------------------------------------- energy budget

"""
    energy_budget(ps, B, n, closure, g; de2 = 0.0)
        -> (kinetic, magnetic, electron_internal, electron_inertia, total)

Total conserved energies of a hybrid state: ion kinetic `Σ_p ½ m w |v|²`
([`kinetic_energy`](@ref)), magnetic `∫ ½|B|² dV` ([`magnetic_energy`](@ref)),
and the electron closure energy ([`electron_internal_energy`](@ref)) — the
internal energy `∫ p_e/(γ−1) dV` for a polytropic closure (γ≠1), the electron
free energy `T_e ∫ n ln n dV` (its exact γ→1 limit) for an isothermal closure,
and the gyrotropic internal energy `∫ (p_⊥ + p_∥/2) dV` for a CGL closure.
`total` is their sum, and at `η = 0` it is conserved by the continuum model —
exactly for the scalar closures, and modulo the anisotropic battery term
`∝ (p_⊥ − p_∥)` for CGL.

With finite electron mass (`HybridModel(...; de2 > 0)`), pass the same `de2`
here: the electron-inertia reservoir `∫ d_e² |J|²/2 dV` (`J = ∇×B` via the
spectral curl; physically the electron flow kinetic energy `½ n m_e |u_e|²` on
the `n ≈ 1` background the inertia filter assumes) is added as
`electron_inertia`. Without it, `total` in a multi-D `de2 > 0` run appears
non-conserved by up to the fluctuation magnetic energy as energy exchanges
between `B` and the electron flow. The default `de2 = 0.0` contributes
exactly `0` (no behavior change).

With resistivity `η > 0`, `total` decays at the Ohmic rate
([`resistive_dissipation`](@ref)); that loss is not deposited in any tracked
reservoir — see the [`resistive_dissipation`](@ref) docstring.

`B` is the 3-tuple of magnetic-field arrays (e.g. `stepper.fields.B`), `n` the
electron/ion density array, `closure` the `ElectronClosure`, and `g` the
`FourierGrid`.
"""
function energy_budget(
    ps::ParticleSet{D,T},
    B::NTuple{3,<:Array{T,D}},
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T};
    de2::Real = 0.0,
) where {D,T}
    _require_grid_tuple(:B, B, g)
    _require_grid_array(:n, n, g)
    de2T = _require_finite_nonnegative_real("de2", de2, T)
    ek = kinetic_energy(ps)
    em = magnetic_energy(B, g)
    ei = electron_internal_energy(n, B, closure, g)
    ej = zero(T)
    if de2T > zero(T)
        J = ntuple(_ -> similar(B[1]), 3)
        curl!(J, B, g)
        s = zero(T)
        @inbounds for I in eachindex(J[1], J[2], J[3])
            s += J[1][I]^2 + J[2][I]^2 + J[3][I]^2
        end
        ej = T(0.5) * de2T * s * prod(g.dx)
    end
    return (
        kinetic = ek,
        magnetic = em,
        electron_internal = ei,
        electron_inertia = ej,
        total = ek + em + ei + ej,
    )
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

With resistivity `η > 0` the total ion momentum is **not** conserved:
`dP/dt = η ∫ n J dV` (the ions are pushed with the full Ohm `E` including
`ηJ`, and the massless electron fluid cannot hold the recoil), which is
nonzero — second order but secular — whenever density and current
fluctuations correlate. A steady drift of `particle` in a resistive run is
un-budgeted model physics, not a pusher bug; only `η = 0` runs conserve
`particle` to particle noise.
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

Per-cell electromagnetic work rate `J·E = Σ_c J_c E_c` — the local rate of
energy transfer from the field to the **plasma** (ions *plus* electron fluid),
since `J = ∇×B` is the total current. The ion (particle) share is
`J_i·E = n u_i·E`; the remainder `J_e·E = u_e·∇p_e − η n u_e·J` is electron
pressure work plus Ohmic exchange, so `∫ J·E dV` is generally **not** the ion
heating rate (the two differ by exactly `∫ u_e·∇p_e dV` at `η = 0`). `J` and
`E` are 3-tuples of grid arrays of identical shape. Integrating the result
against `prod(g.dx)` reproduces [`electric_work`](@ref)`(J, E, g)`.
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

Volume-integrated Ohmic dissipation rate `∫ η |J|² dV = η Σ_g |J_g|² · prod(g.dx)`,
the rate at which resistivity **removes** energy from the tracked budget
([`energy_budget`](@ref)`.total` decays at approximately this rate). Physically
this is electron heating, but the algebraic closures (isothermal / polytropic /
CGL) slave `p_e` to `(n, B)` and cannot receive it, so the energy leaves the
model. Do **not** also add it to the electron-internal ledger — that would
double-book the loss. `J` is the 3-tuple of current-density arrays, `η` the
(scalar) resistivity, `g` the `FourierGrid`.
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
