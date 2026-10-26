# budgets.jl — energy & momentum budget diagnostics (§6/§29).
#
# Thin, composable wrappers over the existing conserved-total diagnostics
# (kinetic_energy, magnetic_energy, electron_internal_energy, total_momentum,
# electric_work) plus the local electromagnetic-work and Ohmic-heating rates.
# Everything is in Ω_ci-normalized units; volume integrals use prod(g.dx).

# ---------------------------------------------------------------- energy budget

function _electron_inertia_energy(J, de2::T, dx) where {T}
    if T === Float32 && all(a -> eltype(a) === Float32, J)
        value = _field_energy_float32((), J, Float64(de2), dx)
        value === nothing || return value
    end
    s = zero(T)
    minimum_j2 = T(Inf)
    @inbounds for I in eachindex(J[1], J[2], J[3])
        x, y, z = J[1][I], J[2][I], J[3][I]
        x2, y2, z2 = x^2, y^2, z^2
        s += x2 + y2 + z2
        minimum_j2 = min(
            minimum_j2,
            ifelse(iszero(x), T(Inf), x2),
            ifelse(iszero(y), T(Inf), y2),
            ifelse(iszero(z), T(Inf), z2),
        )
    end
    half_coefficient = T(0.5) * de2
    unscaled = half_coefficient * s
    volume = prod(dx)
    ordinary = unscaled * volume
    (T === Float32 || T === Float64) || return ordinary
    # A zero reduction plus untouched nonzero-square sentinel proves every
    # current component is zero; an underflowed nonzero square clears it.
    if iszero(s) && isinf(minimum_j2) && isfinite(de2) && all(isfinite, dx)
        return iszero(ordinary) ? ordinary : zero(T)
    end
    lost =
        minimum_j2 < floatmin(T) ||
        _spectral_product_lost_range(T(0.5), de2, half_coefficient) ||
        _spectral_product_lost_range(half_coefficient, s, unscaled) ||
        _spectral_product_lost_range(unscaled, volume, ordinary)
    for first = 1:length(dx)
        partial = one(T)
        for last = first:length(dx)
            product = partial * dx[last]
            lost |= _spectral_product_lost_range(partial, dx[last], product)
            partial = product
        end
    end
    if !lost && !iszero(s) && !iszero(de2) && !iszero(volume)
        # The shared positive-square bound already encloses the reduction and
        # volume; the coefficient contributes its outward power-of-two bound.
        upper = _field_energy_upper_exponent(J, (), s, dx) + exponent(abs(de2)) + 1
        lost = upper > exponent(floatmax(T))
    end
    if lost && isfinite(de2) && all(isfinite, dx) && all(a -> all(isfinite, a), J)
        value = _field_square_sum_exact(J) * _exact_rational(de2) * prod(_exact_rational, dx) / 2
        return _from_exact(T, value)
    end
    return ordinary
end

"""
    energy_budget(ps, B, n, closure, g; de2 = 0.0, pe = nothing)
    energy_budget(species::AbstractVector, B, n, closure, g; de2 = 0.0, pe = nothing)
        -> (kinetic, magnetic, electron_internal, electron_inertia, total)

Energy reservoirs of a hybrid state: ion kinetic `Σ_p ½ m w |v|²`
([`kinetic_energy`](@ref)), magnetic `∫ ½|B|² dV` ([`magnetic_energy`](@ref)),
and the electron closure energy ([`electron_internal_energy`](@ref)) — the
internal energy `∫ p_e/(γ−1) dV` for a polytropic closure (γ≠1), the electron
free energy `T_e ∫ n ln n dV` for an isothermal closure,
and the gyrotropic internal energy `∫ (p_⊥ + p_∥/2) dV` for a CGL closure.
For [`EvolvingElectrons`](@ref), pass the stored state as `pe=stepper.fields.pe`;
its internal energy is `∫ pe/(γ−1) dV`. Missing stored pressure is an error.
Algebraic closures reject an alternate `pe` argument.
`total` is their sum. For periodic ideal massless scalar-closure evolution
(`η = ηH = de2 = 0`) with inactive density floors, the continuum model conserves
this total. This is not an exact discrete-conservation guarantee for the explicit
particle/field integrator. CGL additionally has the anisotropic battery term
`∝ (p_⊥ − p_∥)`.

With finite electron mass (`HybridModel(...; de2 > 0)`), pass the same `de2`
here: the quadratic reservoir `∫ d_e² |J|²/2 dV` (`J = ∇×B` via the spectral
curl) is added as `electron_inertia`. It closes the ideal field-only balance
with fixed ions and uniform density. The implemented reduced field filter is
not the full electron momentum equation: for moving ions its balance retains
`d_e² ∫ J_i·∂tJ dV`, where `J_i = n u_i`. Adding this reservoir therefore does
not establish general total-energy conservation. The default `de2 = 0.0`
contributes exactly `0`.

For finite Float32 currents, inertia integration retains the computed spectral
current, validated `de2` and individual cell spacings through final rounding.
Float64 retains its ordinary summation order and complete-input range repair.
This does not change the spectral differentiation. With `de2 > 0`, concurrent budget calls
require separate grids because the curl uses each grid's spectral workspace.

For scalar closures with `de2 = 0` and inactive floors, dissipative Ohm terms
`R = ηJ − ηH∇²J` contribute `∫ (J_i − J)·R dV` to the continuum total rate.
The ion pusher receives the full electric field, including R, while the
algebraic electron pressure has no independently evolving heat reservoir.
Consequently [`resistive_dissipation`](@ref) alone is not the net tracked-energy
loss. See `docs/HYBRID_ENERGY_MODEL_PLAN.md` for the balance and witnesses.
[`TransportElectrons`](@ref) removes `ηJ` from the ion force and heats the
stored pressure by `η|J|²`, so its periodic continuum total is conserved; pass
`pe = stepper.fields.pe`. The explicit integrator closes it only to its order.

`B` is the 3-tuple of magnetic-field arrays (e.g. `stepper.fields.B`), `n` the
electron/ion density array, `closure` the `ElectronClosure`, and `g` the
`FourierGrid`. For a collection, particle energies are summed using each
population's mass and weights; shared magnetic, electron and inertia reservoirs
are counted once. A typed empty collection contributes zero particle energy.
The diagnostic uses supplied particle velocities without modifying or
time-centering them. Hybrid velocities are staggered relative to B and pe;
compare reservoirs at a common time when assessing energy convergence.
"""
function energy_budget(
    ps::_HybridPopulations{D,T},
    B::NTuple{3,<:Array{T,D}},
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T};
    de2::Real = 0.0,
    pe = nothing,
) where {D,T}
    _require_grid_tuple(:B, B, g)
    _require_grid_array(:n, n, g)
    de2T = _require_finite_nonnegative_real("de2", de2, T)
    _evolves_pressure(closure) &&
        !iszero(de2T) &&
        throw(ArgumentError("$(nameof(typeof(closure))) energy requires de2 = 0"))
    ek = kinetic_energy(ps)
    em = magnetic_energy(B, g)
    ei = electron_internal_energy(n, B, closure, g; pe)
    ej = zero(T)
    if de2T > zero(T)
        J = ntuple(_ -> similar(B[1]), 3)
        curl!(J, B, g)
        ej = _electron_inertia_energy(J, de2T, g.dx)
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
    momentum_budget(species::AbstractVector, B, g)

Momentum budget of a hybrid state. `particle` is the total ion momentum
`Σ_p m w v` ([`total_momentum`](@ref)), summed across populations for a collection.
A typed empty collection contributes zero momentum. `total` adds the electromagnetic field
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
`particle` to particle noise. [`TransportElectrons`](@ref) pushes ions with
the force field without `ηJ`, so this resistive drift term is absent there.
"""
function momentum_budget(
    ps::_HybridPopulations{D,T},
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

Volume-integrated Ohmic quadratic form
`∫ η |J|² dV = η Σ_g |J_g|² · prod(g.dx)`.

This is not generally the net loss of [`energy_budget`](@ref)`.total`. With
scalar algebraic pressure, `ηH = de2 = 0` and inactive density floors, the
continuum Ohmic contribution to that total is
`η ∫ J_i·J dV − ∫ η|J|² dV`, because the ion pusher also receives `ηJ`.
The algebraic closures have no separate heat reservoir. Adding this quadratic
form as electron heating alone would leave the ion-work term; a conservative
heating extension also needs a consistent ion recoil model.
[`TransportElectrons`](@ref) provides that pairing: its ions omit `ηJ` and its
stored pressure gains this heat, so the quadratic form is its Ohmic exchange.

`J` is the 3-tuple of current-density arrays, `η` the scalar resistivity, and
`g` the `FourierGrid`. The returned numerical quantity is unchanged by these
energy-accounting qualifications.
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

_empty_particle_momentum(::Type{<:ParticleSet{D,T}}) where {D,T} = (zero(T), zero(T), zero(T))
_empty_particle_momentum(::Type{<:ParticleSet}) = throw(
    ArgumentError("empty momentum collection must specify its particle dimension and precision"),
)

"Total momentum over populations; typed empty collections return a zero tuple."
function total_momentum(species::AbstractVector{P}) where {P<:ParticleSet}
    isempty(species) && return _empty_particle_momentum(P)
    p = total_momentum(first(species))
    for ps in Iterators.drop(species, 1)
        q = total_momentum(ps)
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
