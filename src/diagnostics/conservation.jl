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

function _budget_field_state(value::T, fields, dx, coefficient_operations = 0) where {T}
    if iszero(value) && all(a -> all(iszero, a), fields)
        return _energy_zero(value, Val(true))
    end
    T === Float32 && return _energy_rounded(value, Val(true))
    count = Int128(length(first(fields))) * length(fields) + length(dx) + 4 + coefficient_operations
    lower, upper = _energy_positive_interval(value, count, precision(T))
    return _energy_result(value, lower, upper, Val(true))
end

function _budget_electron_state(n, B, closure, g, pe)
    return _energy_unknown(electron_internal_energy(n, B, closure, g; pe), Val(true))
end
function _budget_electron_state(
    n::Array{T,D},
    B,
    closure::Union{IsothermalElectrons{C},PolytropicElectrons{C}},
    g,
    pe,
) where {T<:Union{Float32,Float64},D,C<:Union{Float32,Float64}}
    return _electron_internal_energy_checked(n, closure, g, pe, Val(true))
end
function _budget_electron_state(
    n::Array{T,D},
    B,
    closure::Union{EvolvingElectrons,_ScalarTransportClosure},
    g,
    pe,
) where {T<:Union{Float32,Float64},D}
    return _electron_internal_energy_checked(n, closure, g, pe, Val(true))
end
function _budget_electron_state(
    n::Array{T,D},
    B,
    closure::CGLElectrons{C},
    g,
    pe,
) where {T<:Union{Float32,Float64},D,C<:Union{Float32,Float64}}
    pe === nothing ||
        throw(ArgumentError("algebraic electron energy does not accept an alternate pe"))
    return _cgl_energy_checked(n, B, closure, g.dx, Val(true))
end

@inline function _budget_total_certified(total::T, states) where {T<:Union{Float32,Float64}}
    isfinite(total) || return false
    lower = upper = 0.0
    active = false
    for state in states
        state === nothing && continue
        iszero(state.lower) && iszero(state.upper) && continue
        isfinite(state.lower) && isfinite(state.upper) || return false
        # Exact zero intervals need no addition allowance; the first active
        # interval is already outward and can initialize the complete sum.
        lower = active ? prevfloat(lower + state.lower) : state.lower
        upper = active ? nextfloat(upper + state.upper) : state.upper
        active = true
    end
    active || return iszero(total)
    minimum, maximum = Float64(floatmin(T)), Float64(floatmax(T))
    return total > 0 ? minimum <= lower <= upper <= maximum :
           total < 0 && -maximum <= lower <= upper <= -minimum
end

_budget_closure_eligible(::ElectronClosure) = false
_budget_closure_eligible(c::IsothermalElectrons{<:Union{Float32,Float64}}) = isfinite(c.Te)
_budget_closure_eligible(c::PolytropicElectrons{<:Union{Float32,Float64}}) =
    all(isfinite, (c.pe0, c.n0, c.γ)) && c.pe0 >= 0 && c.n0 > 0 && c.γ > 0
_budget_closure_eligible(c::CGLElectrons{<:Union{Float32,Float64}}) =
    isfinite(c.Cperp) && isfinite(c.Cpar)
_budget_closure_eligible(::Union{EvolvingElectrons,_ScalarTransportClosure}) = true

# Only an inconclusive coarse whole-budget certificate requests this work.
_budget_refine_electron_state(state, closure, dx) = state
function _budget_refine_electron_state(state, closure::PolytropicElectrons, dx)
    hasproperty(state, :minimum_density) || return state
    lower = _polytropic_energy_lower_bound(state.minimum_density, state.active_count, closure, dx)
    return signbit(state.value) ? merge(state, (; upper = min(state.upper, -lower))) :
           merge(state, (; lower = max(state.lower, lower)))
end

function _budget_complete_total(
    ps,
    B,
    n::Array{T},
    closure,
    g,
    de2,
    pe,
    J,
    ordinary,
    electron,
) where {T<:Union{Float32,Float64}}
    _budget_closure_eligible(closure) || return ordinary
    species = ps isa ParticleSet ? (ps,) : ps
    finite_particles = all(species) do population
        isfinite(population.m) &&
            all(isfinite, population.weight) &&
            all(a -> all(isfinite, a), population.v)
    end
    finite_particles === true &&
        all(isfinite, n) &&
        all(isfinite, g.dx) &&
        all(a -> all(isfinite, a), B) &&
        (J === nothing || all(a -> all(isfinite, a), J)) || return ordinary
    volume = prod(_exact_rational, g.dx)
    quadratic = _kinetic_energy_rational(species) + _field_square_sum_exact(B) * volume / 2
    J === nothing || (quadratic += _exact_rational(de2) * _field_square_sum_exact(J) * volume / 2)
    if closure isa Union{EvolvingElectrons,_ScalarTransportClosure}
        return _from_exact(
            T,
            quadratic + _stored_pressure_energy_rational(pe, electron.gamma, g.dx),
        )
    elseif closure isa IsothermalElectrons
        return _isothermal_exact_energy(n, closure.Te, g.dx, quadratic)
    elseif closure isa PolytropicElectrons
        if closure.γ == 1
            coefficient = _exact_rational(closure.pe0) / _exact_rational(closure.n0)
            return _isothermal_exact_energy(n, coefficient, g.dx, quadratic)
        end
        return _polytropic_exact_energy(n, closure, g.dx, quadratic)
    elseif closure isa CGLElectrons
        return _cgl_energy_exact(n, B, closure, g.dx, ordinary, quadratic)
    end
    return ordinary
end

@inline function _energy_budget_finish(
    ps,
    B,
    n::Array{T},
    closure,
    g,
    de2,
    pe,
    J,
    components,
    states,
) where {T}
    ek, em, ei = components
    ej = J === nothing ? zero(T) : _electron_inertia_energy(J, de2, g.dx)
    total = ek + em + ei + ej
    if states !== nothing && total isa T && ei isa T
        inertia = J === nothing ? _energy_zero(ej, Val(true)) : _budget_field_state(ej, J, g.dx, 3)
        if !_budget_total_certified(total, (states..., inertia))
            electron = _budget_refine_electron_state(states[3], closure, g.dx)
            if electron === states[3] ||
               !_budget_total_certified(total, (states[1], states[2], electron, inertia))
                total = _budget_complete_total(ps, B, n, closure, g, de2, pe, J, total, electron)
            end
        end
    end
    return (
        kinetic = ek,
        magnetic = em,
        electron_internal = ei,
        electron_inertia = ej,
        total = total,
    )
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
`total` retains the ordinary associated sum when complete reservoir bounds certify
normal range and sign. Range-sensitive finite Float32/64 built-in states combine
original reservoir expressions before final rounding; `total` can therefore be
finite while reported components are infinite, or differ from their rounded sum.
True complete overflow, nonfinite/custom arithmetic and the CGL singularity remain
visible. This changes no model law or exact-conservation guarantee.
For periodic ideal massless scalar-closure evolution
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
[`TransportElectrons`](@ref) removes R from the ion force and adds positive
Ohmic/hyperresistive heat, pairing the periodic magnetic loss; pass stored
`pe = stepper.fields.pe`. The explicit integrator closes the total to its order.

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
    supported = T === Float32 || T === Float64
    kinetic = magnetic = electron = nothing
    if supported
        kinetic =
            ps isa ParticleSet ? _kinetic_energy_checked((ps,), T, Val(true), Val(true)) :
            _kinetic_energy_checked(ps, T, Val(false), Val(true))
        ek = kinetic.value
    else
        ek = kinetic_energy(ps)
    end
    em = magnetic_energy(B, g)
    if supported
        magnetic = _budget_field_state(em, B, g.dx)
        electron = _budget_electron_state(n, B, closure, g, pe)
        ei = electron.value
    else
        ei = electron_internal_energy(n, B, closure, g; pe)
    end
    components = (ek, em, ei)
    states = supported ? (kinetic, magnetic, electron) : nothing
    if de2T > zero(T)
        J = ntuple(_ -> similar(B[1]), 3)
        curl!(J, B, g)
        # Keep concrete current and no-current completion call sites separate;
        # both consume the same single computed spectral-current candidate.
        return _energy_budget_finish(ps, B, n, closure, g, de2T, pe, J, components, states)
    end
    return _energy_budget_finish(ps, B, n, closure, g, de2T, pe, nothing, components, states)
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
the force field without `ηJ − ηH∇²J`, so dissipative drift terms are absent there.
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

"""
    electron_conduction_diagnostics(stepper, species)

Instantaneous field-aligned conduction budget at the current integer-time
pressure and B, using a private fresh deposition of the supplied ions. Returns
exact rational integrals of the computed, unrounded conduction source, prepared
channel dissipation and measured adjoint defect. This is neither the last RK
stage's history nor the rounded complete pressure source. FFTW and the two
preparations can produce a nonzero defect. Public state and RNG are unchanged;
the grid's usual nonconcurrent FFT scratch contract still applies.
"""
function electron_conduction_diagnostics(
    st::HybridStepper{D,T},
    species::_HybridPopulations{D,T},
) where {D,T}
    closure, g = st.model.closure, st.g
    closure isa FieldAlignedConductingElectrons ||
        throw(ArgumentError("conduction diagnostics require FieldAlignedConductingElectrons"))
    _validate_pressure_step_state(st, species)
    density = zeros(T, g.n)
    if species isa ParticleSet
        density!(density, species, g, st.shape)
    else
        _validate_hybrid_species(species)
        number_density = similar(density)
        for population in species
            density!(number_density, population, g, st.shape)
            # Match the density row of _accumulate_charge_moments!, including
            # its stored per-population accumulation order, without computing
            # unused momentum or allocating a full public-field workspace.
            @. density =
                getproperty(_charge_moment_parts(density, population.q, number_density), :value)
        end
    end
    _validate_transport_density(density, T(st.model.nfloor), "conduction diagnostic")
    Q = Rational{BigInt}
    source = zeros(Q, g.n)
    quotient = zeros(Q, g.n)
    coefficient = _exact_rational(_pressure_heat_gamma(T(closure.γ)))
    dissipation = zero(Q)
    if !iszero(_electron_conductivity(closure, T))
        _conduction_quotient!(quotient, st.fields.pe, density, g)
        dissipation = _field_aligned_conduction!(
            source,
            quotient,
            coefficient,
            closure,
            st.fields.B,
            _field_aligned_storage(g, closure),
            g;
            budget = true,
        )
    end
    volume = prod(_exact_rational, g.dx)
    variance_rate = volume * sum(I -> quotient[I] * source[I], eachindex(source))
    channel_dissipation = volume * dissipation
    return (
        scope = :instantaneous,
        time = st.time[],
        step = st.step[],
        pressure_time_offset = zero(T),
        magnetic_time_offset = zero(T),
        internal_energy_rate = volume * sum(source) / (_exact_rational(T(closure.γ)) - 1),
        variance_rate = variance_rate,
        channel_dissipation = channel_dissipation,
        adjoint_defect = variance_rate + coefficient * channel_dissipation,
    )
end

"""
    electron_ion_exchange_diagnostics(stepper::CollisionalHybridStepper, species)

Return an independently owned numeric receipt for the last completed reduced
electron-ion thermal source, or `nothing` immediately after initialization.
Recheck its physical coefficients, actual stored endpoint reservoirs, rounding
certificates and moment-entropy bound. This is source-split history, not an
instantaneous collision rate or a conservation claim for the field integrator.
The diagnostic performs private deposition/bins without drawing RNG values.
NGP returns its policy-1 per-node receipt; CIC/TSC return policy-2 ordered
overlap edge history, mesh/source defects and telescoping native budgets.
"""
function electron_ion_exchange_diagnostics(st::CollisionalHybridStepper, species::AbstractVector)
    _acquire_collisional!(st)
    try
        _validate_collisional_binding(st, species)
        _validate_collisional_boundary(st)
        st.electron_ion === nothing && throw(
            ArgumentError(
                "electron-ion exchange diagnostics require an enabled thermal configuration",
            ),
        )
        _validate_electron_ion_report(
            st.thermal_report,
            st.electron_ion,
            st.inner,
            st.species,
            st.inner.fields,
            st.source_time,
            st.source_count,
            st.dt,
            st.epsilon0,
            st.weight_to_number,
        )
        return deepcopy(st.thermal_report)
    finally
        st.busy[] = false
    end
end

function electron_conduction_diagnostics(st::CollisionalHybridStepper, species::AbstractVector)
    _acquire_collisional!(st)
    try
        _validate_collisional_binding(st, species)
        _validate_collisional_boundary(st)
        return electron_conduction_diagnostics(st.inner, st.species)
    finally
        st.busy[] = false
    end
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
energy-accounting qualifications. Finite Float32/Float64 inputs retain the
complete integral through final rounding when an intermediate loses range;
ordinary inputs retain the existing summation order.
"""
function resistive_dissipation(
    J::NTuple{3,<:AbstractArray{T}},
    η::Real,
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_tuple(:J, J, g)
    ηT = _require_finite_nonnegative_real("η", η, T)
    s = zero(T)
    minimum_j2 = T(Inf)
    @inbounds for I in eachindex(J[1], J[2], J[3])
        x, y, z = J[1][I], J[2][I], J[3][I]
        x2, y2, z2 = x * x, y * y, z * z
        s += x2 + y2 + z2
        minimum_j2 = min(
            minimum_j2,
            ifelse(iszero(x), T(Inf), x2),
            ifelse(iszero(y), T(Inf), y2),
            ifelse(iszero(z), T(Inf), z2),
        )
    end
    weighted = ηT * s
    volume = prod(g.dx)
    result = weighted * volume
    (T === Float32 || T === Float64) || return result
    lost =
        minimum_j2 < floatmin(T) ||
        _spectral_product_lost_range(ηT, s, weighted) ||
        _spectral_product_lost_range(weighted, volume, result)
    # A spacing product may lose range before a later factor restores it.
    # Contiguous subproducts cover either product association in D <= 3.
    for first = 1:length(g.dx)
        partial = one(T)
        for last = first:length(g.dx)
            product = partial * g.dx[last]
            lost |= _spectral_product_lost_range(partial, g.dx[last], product)
            partial = product
        end
    end
    if !lost && !iszero(s) && !iszero(ηT) && !iszero(volume)
        upper = _field_energy_upper_exponent(J, (), s, g.dx) + exponent(ηT) + 1
        lost = upper > exponent(floatmax(T))
    end
    if lost && all(isfinite, g.dx) && all(a -> all(isfinite, a), J)
        integral = _field_square_sum_exact(J) * _exact_rational(ηT) * prod(_exact_rational, g.dx)
        return _from_exact(T, integral)
    end
    return result
end

# --- particle momentum / electric work (from diagnostics.jl) ---
# Exact original-factor totals for collective diagnostics. The ordinary serial
# methods below retain their native summation and rounding order.
function _momentum_rational(species)
    px = py = pz = zero(Rational{BigInt})
    for ps in species
        mass = _exact_rational(ps.m)
        for p in eachindex(ps.weight)
            factor = mass * _exact_rational(ps.weight[p])
            iszero(factor) && continue
            px += factor * _exact_rational(ps.v[1][p])
            py += factor * _exact_rational(ps.v[2][p])
            pz += factor * _exact_rational(ps.v[3][p])
        end
    end
    return px, py, pz
end

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
