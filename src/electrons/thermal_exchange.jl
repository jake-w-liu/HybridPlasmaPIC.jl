# Reduced scalar Maxwellian electron/ion thermal moment exchange. The wrapper
# owns clocks, RNG and publication; this owner prepares one complete source.
const _ElectronIonQ = Rational{BigInt}

"""
    ElectronIonThermalExchange(; electron_mass, coulomb_log, log_bits=nothing)

Opt in to two-way stored-electron-pressure/ion-thermal-energy exchange through
`CollisionalHybridStepper(...; electron_ion=request)`. `electron_mass` is in
the same normalized mass units as the ions; `coulomb_log` is a nonnegative
scalar or one value per ion species. Grid binding owns Float32/64 coefficients.
This Maxwellian moment approximation omits drift friction and electron inertia.
See `docs/ELECTRON_ION_THERMAL_EXCHANGE.md` for geometry, split and rejection rules.
"""
struct ElectronIonThermalExchange{M,L}
    electron_mass::M
    coulomb_log::L
    log_bits::Union{Nothing,Int}
end

function ElectronIonThermalExchange(; electron_mass::Real, coulomb_log, log_bits = nothing)
    isfinite(electron_mass) && electron_mass > 0 ||
        throw(ArgumentError("electron_mass must be finite and positive"))
    if coulomb_log isa Real
        isfinite(coulomb_log) && coulomb_log >= 0 ||
            throw(ArgumentError("electron-ion Coulomb logarithms must be finite and nonnegative"))
    elseif coulomb_log isa AbstractVector{<:Real}
        axes(coulomb_log, 1) == Base.OneTo(length(coulomb_log)) ||
            throw(ArgumentError("electron-ion Coulomb logarithms require one-based axes"))
        all(x -> isfinite(x) && x >= 0, coulomb_log) ||
            throw(ArgumentError("electron-ion Coulomb logarithms must be finite and nonnegative"))
    else
        throw(ArgumentError("electron-ion Coulomb logarithms must be a real scalar or vector"))
    end
    if log_bits !== nothing
        log_bits isa Integer && !(log_bits isa Bool) && 0 < log_bits <= typemax(Int) ||
            throw(ArgumentError("log_bits must be a positive non-Bool integer fitting Int"))
        log_bits = Int(log_bits)
    end
    return ElectronIonThermalExchange(deepcopy(electron_mass), deepcopy(coulomb_log), log_bits)
end

struct _ElectronIonThermalConfig{T,N}
    electron_mass::T
    logarithm::Vector{T}
    binding::NTuple{N,T}
    log_bits::Int
    policy::Int
end

mutable struct _ElectronIonThermalWork{D,T}
    bins::CoulombCollisionWorkspace{T}
    density::Array{T,D}
    number_density::Array{T,D}
    pressure::Array{T,D}
    pending::Union{Nothing,NamedTuple}
end

_bind_electron_ion(::Nothing, inner, species) = nothing
function _bind_electron_ion(request::ElectronIonThermalExchange, inner::HybridStepper{D,T}, species) where {D,T}
    _evolves_pressure(inner.model.closure) && inner.shape isa NGP && iszero(inner.model.de2) ||
        throw(ArgumentError("electron-ion thermal exchange requires stored scalar pressure, NGP and de2 = 0"))
    T <: Union{Float32,Float64} ||
        throw(ArgumentError("electron-ion thermal exchange requires Float32 or Float64"))
    T(inner.model.closure.γ) == T(5 / 3) ||
        throw(ArgumentError("electron-ion thermal exchange requires grid gamma = T(5/3)"))
    mass = _require_finite_positive_real("electron_mass", request.electron_mass, T)
    input = request.coulomb_log
    input isa Real || length(input) == length(species) ||
        throw(DimensionMismatch("electron-ion Coulomb logarithms require one value per ion species"))
    logs = T[_electron_ion_log(input isa Real ? input : input[a], T) for a in eachindex(species)]
    bits = request.log_bits === nothing ? precision(T) + 16 : request.log_bits
    precision(T) + 16 <= bits <= 4096 ||
        throw(ArgumentError("log_bits must be between precision(T)+16 and 4096"))
    return _ElectronIonThermalConfig(mass, logs, Tuple(logs), bits, 1)
end
_bind_electron_ion(request, inner, species) =
    throw(ArgumentError("electron_ion must be nothing or ElectronIonThermalExchange"))

function _electron_ion_log(value, ::Type{T}) where {T}
    value isa Real && isfinite(value) && value >= 0 ||
        throw(ArgumentError("electron-ion Coulomb logarithms must be finite and nonnegative"))
    result = _require_finite_nonnegative_real("electron-ion Coulomb logarithm", value, T)
    !iszero(value) && iszero(result) &&
        throw(ArgumentError("nonzero electron-ion Coulomb logarithms must survive grid conversion"))
    return result
end

_electron_ion_work(::Nothing, inner, species) = nothing
function _electron_ion_work(::_ElectronIonThermalConfig, inner::HybridStepper{D,T}, species) where {D,T}
    return _ElectronIonThermalWork(
        CoulombCollisionWorkspace(species, inner.g), zeros(T, inner.g.n),
        zeros(T, inner.g.n), zeros(T, inner.g.n), nothing,
    )
end

_electron_ion_active(::Nothing) = false
_electron_ion_active(config::_ElectronIonThermalConfig) = any(!iszero, config.logarithm)
_electron_ion_configuration(::Nothing) = nothing
_electron_ion_configuration(config::_ElectronIonThermalConfig) = (
    policy = config.policy, electron_mass = config.electron_mass,
    logarithm = config.binding, log_bits = config.log_bits,
)

function _validate_electron_ion_config(config::_ElectronIonThermalConfig{T}, inner, species) where {T}
    config.policy == 1 && isfinite(config.electron_mass) && config.electron_mass > 0 ||
        throw(ArgumentError("invalid electron-ion thermal policy or electron mass"))
    length(config.logarithm) == length(species) == length(config.binding) &&
        all(a -> isequal(config.logarithm[a], config.binding[a]) && isfinite(config.logarithm[a]) && config.logarithm[a] >= 0, eachindex(config.logarithm)) ||
        throw(ArgumentError("electron-ion Coulomb logarithms are fixed"))
    precision(T) + 16 <= config.log_bits <= 4096 ||
        throw(ArgumentError("invalid electron-ion entropy precision"))
    _evolves_pressure(inner.model.closure) && inner.shape isa NGP &&
        T(inner.model.closure.γ) == T(5 / 3) && iszero(inner.model.de2) ||
        throw(ArgumentError("electron-ion thermal model and shape are fixed"))
    return nothing
end

# Match the canonical nearest-node stencil, including its ties-up and periodic
# policies. The existing ion-ion floor-cell bins deliberately stay unchanged.
function _electron_ion_node(ps::ParticleSet{D}, g, p) where {D}
    index = ntuple(d -> mod(_stencil1d(NGP(), _particle_cell_position(ps, g, d, p))[1], g.n[d]) + 1, D)
    return LinearIndices(g.n)[index...]
end

function _electron_ion_density!(work, species, inner)
    fill!(work.density, zero(eltype(work.density)))
    for ps in species
        density!(work.number_density, ps, inner.g, inner.shape)
        @. work.density = getproperty(_charge_moment_parts(work.density, ps.q, work.number_density), :value)
    end
    _validate_transport_density(work.density, eltype(work.density)(inner.model.nfloor), "electron-ion exchange")
    return nothing
end

function _electron_ion_stats(ps, bins, cell, velocity = ps.v)
    Q = _ElectronIonQ
    weight = kinetic = zero(Q)
    momentum = (zero(Q), zero(Q), zero(Q))
    count = bins.offset[cell + 1] - bins.offset[cell]
    mass = _exact_rational(ps.m)
    for k in bins.offset[cell]:(bins.offset[cell + 1] - 1)
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        weight += w
        v = ntuple(c -> _exact_rational(velocity[c][p]), 3)
        momentum = ntuple(c -> momentum[c] + mass * w * v[c], 3)
        kinetic += mass * w * sum(x -> x * x, v) / 2
    end
    mean = iszero(weight) ? (zero(Q), zero(Q), zero(Q)) : ntuple(c -> momentum[c] / (mass * weight), 3)
    thermal = kinetic - mass * weight * sum(x -> x * x, mean) / 2
    thermal >= 0 || throw(ArgumentError("electron-ion thermal energy must be nonnegative"))
    return (; weight, mean, kinetic, momentum, thermal, capacity = 3weight / 2, count)
end

# Fixed work in the new coefficient owner; the old adaptive pi regime keeps its
# exact same loop and comparisons through the shared literal refinement row.
function _electron_ion_pi_bounds(::Type{T}) where {T}
    a = b = zero(_ElectronIonQ)
    ta, tb = BigInt(1) // BigInt(5), BigInt(1) // BigInt(239)
    low = high = zero(_ElectronIonQ)
    for n = 0:(cld(precision(T) + 18, 4) - 1)
        a, b, ta, tb, low, high = _coulomb_pi_refine(a, b, ta, tb, n)
    end
    return low, high
end

function _electron_ion_root_scale(::Type{T}, value::_ElectronIonQ, mode) where {T}
    value > 0 || throw(ArgumentError("electron-ion scaled root requires positive energy"))
    n, d = numerator(value), denominator(value)
    power = ndigits(n; base = 2) - ndigits(d; base = 2)
    if power >= 0 ? n < d << power : n << -power < d
        power -= 1
    end
    exponent = fld(power, 2)
    reduced = _exact_binary_scale(value, -2exponent)
    mantissa = _sqrt_from_exact(T, reduced, Int128(0), mode, false)
    return _exact_binary_scale(mantissa, exponent)
end

function _electron_ion_beta(::Type{T}, K::_ElectronIonQ, S::_ElectronIonQ, pi_bounds) where {T}
    iszero(K) && return (zero(K), zero(K), zero(K))
    K > 0 && S > 0 || throw(ArgumentError("singular mixed cold electron-ion rate"))
    pi_low, pi_high = pi_bounds
    low = _electron_ion_root_scale(T, K^2 / (18pi_high^3 * S^3), RoundDown)
    high = _electron_ion_root_scale(T, K^2 / (18pi_low^3 * S^3), RoundUp)
    b_low, b_high = low / (1 + low), high / (1 + high)
    u = BigInt(1) // (BigInt(1) << precision(T))
    0 < b_low <= b_high < 1 && b_high - b_low <= 8u * b_low ||
        throw(ArgumentError("electron-ion coefficient accuracy is unresolved at grid precision"))
    return b_low, (b_low + b_high) / 2, b_high
end

function _electron_ion_targets(Ce, Te, statistics, beta)
    numerator = Ce * Te
    denominator = Ce
    for a in eachindex(statistics)
        s = statistics[a]
        iszero(s.weight) && continue
        alpha = s.capacity * beta[a][2]
        numerator += alpha * s.thermal / s.capacity
        denominator += alpha
    end
    electron_temperature = numerator / denominator
    targets = _ElectronIonQ[
        iszero(s.weight) ? zero(_ElectronIonQ) :
        (1 - beta[a][2]) * s.thermal + beta[a][2] * s.capacity * electron_temperature
        for (a, s) in enumerate(statistics)
    ]
    return electron_temperature, targets
end

# A conservative half-spacing at each finite endpoint; no tolerance fitted to
# measured physical work. The zero quantum also covers signed-zero publication.
function _electron_ion_radius(value::T) where {T}
    isfinite(value) || throw(ArgumentError("electron-ion output must be finite"))
    iszero(value) && return _exact_rational(nextfloat(zero(T)))
    q = _exact_rational(value)
    before, after = prevfloat(value), nextfloat(value)
    left = isfinite(before) ? q - _exact_rational(before) : _exact_rational(after) - q
    right = isfinite(after) ? _exact_rational(after) - q : left
    return max(left, right) / 2
end

function _electron_ion_publish!(bins, ps, cell, old, target, rng)
    T = eltype(ps.weight)
    Q = _ElectronIonQ
    if target == old.thermal
        return (scale = one(Q), seed_energy = old.thermal, frame_energy = old.kinetic,
                momentum_bound = (zero(Q), zero(Q), zero(Q)), kinetic_bound = zero(Q))
    end
    seed_mean, seed_energy = old.mean, old.thermal
    if iszero(seed_energy)
        old.count >= 2 || throw(ArgumentError("a single positive marker cannot receive momentum-preserving heat"))
        for k in bins.offset[cell]:(bins.offset[cell + 1] - 1), c = 1:3
            p = bins.index[k]
            draw = randn(rng, T)
            isfinite(draw) || throw(ArgumentError("electron-ion cold seed must be finite"))
            bins.velocity[c][p] = draw
        end
        seed = _electron_ion_stats(ps, bins, cell, bins.velocity)
        seed_mean, seed_energy = seed.mean, seed.thermal
        seed_energy > 0 || throw(ArgumentError("electron-ion cold seed has zero weighted thermal energy"))
    end
    scale = _electron_ion_root_scale(T, target / seed_energy, RoundNearest)
    u = BigInt(1) // (BigInt(1) << precision(T))
    abs(scale^2 * seed_energy - target) <= (2u + u^2) * target ||
        throw(ArgumentError("electron-ion scale accuracy is unresolved"))
    momentum_bound = (zero(Q), zero(Q), zero(Q))
    kinetic_bound = zero(Q)
    mass = _exact_rational(ps.m)
    for k in bins.offset[cell]:(bins.offset[cell + 1] - 1)
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        for c = 1:3
            seed = iszero(old.thermal) ? bins.velocity[c][p] : ps.v[c][p]
            frame = old.mean[c] + scale * (_exact_rational(seed) - seed_mean[c])
            value = _from_exact(T, frame)
            radius = _electron_ion_radius(value)
            abs(_exact_rational(value) - frame) <= radius ||
                throw(ArgumentError("electron-ion component publication exceeded its rounding bound"))
            bins.velocity[c][p] = value
            momentum_bound = ntuple(d -> momentum_bound[d] + (c == d ? mass * w * radius : zero(Q)), 3)
            kinetic_bound += mass * w * (abs(frame) * radius + radius^2 / 2)
        end
    end
    frame_energy = mass * old.weight * sum(x -> x * x, old.mean) / 2 + scale^2 * seed_energy
    return (; scale, seed_energy, frame_energy, momentum_bound, kinetic_bound)
end

function _electron_ion_entropy(capacities, before, after, bits)
    Q = _ElectronIonQ
    all(>(zero(Q)), after) || throw(ArgumentError("electron-ion active reservoirs must publish positive energy"))
    any(iszero, before) && return :cold_positive, (zero(Q), zero(Q))
    terms = Dict{BigInt,Q}()
    function add!(base, coefficient)
        (base == 1 || iszero(coefficient)) && return nothing
        value = get(terms, base, zero(Q)) + coefficient
        iszero(value) ? delete!(terms, base) : (terms[base] = value)
        return nothing
    end
    for a in eachindex(capacities)
        ratio = after[a] / before[a]
        add!(numerator(ratio), capacities[a])
        add!(denominator(ratio), -capacities[a])
    end
    isempty(terms) && return :resolved, (zero(Q), zero(Q))
    logtwo = _isothermal_atanh_log_bounds(BigInt(1) // BigInt(3), bits)
    low = high = zero(Q)
    for (base, coefficient) in terms
        lo, hi = _isothermal_integer_log_bounds(base, bits, logtwo)
        low += coefficient * (coefficient > 0 ? lo : hi)
        high += coefficient * (coefficient > 0 ? hi : lo)
    end
    low >= 0 || throw(ArgumentError("electron-ion actual entropy decrease or unresolved entropy sign"))
    return :resolved, (low, high)
end

# Bounds used when checking an incoming receipt without re-running the source.
# |frame| <= |stored|+rho caps the saved, tighter frame-based energy bound.
function _electron_ion_endpoint_bounds(ps, bins, cell)
    Q = _ElectronIonQ
    momentum = (zero(Q), zero(Q), zero(Q))
    kinetic = zero(Q)
    mass = _exact_rational(ps.m)
    for k in bins.offset[cell]:(bins.offset[cell + 1] - 1), c = 1:3
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        radius = _electron_ion_radius(ps.v[c][p])
        momentum = ntuple(d -> momentum[d] + (c == d ? mass * w * radius : zero(Q)), 3)
        kinetic += mass * w * (abs(_exact_rational(ps.v[c][p])) * radius + 3radius^2 / 2)
    end
    return momentum, kinetic
end

function _electron_ion_coefficients(config, statistics, species, Ce, Te, volume, h, epsilon, number_scale, pi_bounds, ::Type{T}) where {T}
    Q = _ElectronIonQ
    active = [!iszero(config.logarithm[a]) && !iszero(statistics[a].weight) for a in eachindex(species)]
    all_cold = iszero(Te) && all(a -> !active[a] || iszero(statistics[a].thermal), eachindex(species))
    beta = NTuple{3,Q}[]
    me = _exact_rational(config.electron_mass)
    for a in eachindex(species)
        s, ps = statistics[a], species[a]
        if !active[a] || all_cold
            push!(beta, (zero(Q), zero(Q), zero(Q)))
            continue
        end
        mass, charge = _exact_rational(ps.m), _exact_rational(ps.q)
        S = Te / me + s.thermal / (s.capacity * mass)
        K = _exact_rational(h) * (Ce / s.capacity) * _exact_rational(number_scale) *
            (s.weight / volume) * charge^2 * _exact_rational(config.logarithm[a]) /
            (_exact_rational(epsilon)^2 * me * mass)
        push!(beta, _electron_ion_beta(T, K, S, pi_bounds))
    end
    return active, beta, all_cold
end

function _prepare_electron_ion!(config, work, inner, species, h, epsilon, number_scale, rng)
    T = eltype(work.pressure)
    Q = _ElectronIonQ
    _electron_ion_density!(work, species, inner)
    for a in eachindex(species)
        ps = species[a]
        _bin_coulomb_particles!(work.bins.populations[a], ps, p -> _electron_ion_node(ps, inner.g, p))
    end
    copyto!(work.pressure, inner.fields.pe)
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    pi_bounds = _electron_ion_pi_bounds(T)
    plans = NamedTuple[]
    for cell in eachindex(work.pressure)
        statistics = [_electron_ion_stats(species[a], work.bins.populations[a], cell) for a in eachindex(species)]
        Ne = _exact_rational(work.density[cell]) * volume
        Ce = Ne / gamma
        old_pressure = inner.fields.pe[cell]
        isfinite(old_pressure) && old_pressure >= 0 || throw(ArgumentError("electron-ion pressure must be finite and nonnegative"))
        Te = _exact_rational(old_pressure) / _exact_rational(work.density[cell])
        active, beta, all_cold = _electron_ion_coefficients(config, statistics, species, Ce, Te, volume, h, epsilon, number_scale, pi_bounds, T)
        electron_target, targets = _electron_ion_targets(Ce, Te, statistics, beta)
        # Every cell's rate and sparse checks precede every thermal draw.
        for a in eachindex(species)
            iszero(statistics[a].thermal) && targets[a] > 0 && statistics[a].count < 2 &&
                throw(ArgumentError("a single positive marker cannot receive momentum-preserving heat"))
        end
        push!(plans, (; cell, statistics, Ce, Te, old_pressure, active, beta,
            all_cold, electron_target, targets))
    end
    cells = NamedTuple[]
    for plan in plans
        (; cell, statistics, Ce, Te, old_pressure, active, beta,
            all_cold, electron_target, targets) = plan
        Ne = _exact_rational(work.density[cell]) * volume
        ions = NamedTuple[]
        change = zero(Q)
        before_temperature, after_temperature = Q[Te], Q[]
        capacities = Q[Ce]
        for a in eachindex(species)
            ps, bins, old = species[a], work.bins.populations[a], statistics[a]
            materialized = _electron_ion_publish!(bins, ps, cell, old, targets[a], rng)
            actual = _electron_ion_stats(ps, bins, cell, bins.velocity)
            all(c -> abs(actual.momentum[c] - old.momentum[c]) <= materialized.momentum_bound[c], 1:3) &&
                abs(actual.kinetic - materialized.frame_energy) <= materialized.kinetic_bound ||
                throw(ArgumentError("electron-ion actual momentum or kinetic-energy rounding certificate failed"))
            active[a] && !all_cold && actual.thermal <= 0 &&
                throw(ArgumentError("electron-ion target lost positive ion thermal energy"))
            change += actual.kinetic - old.kinetic
            if active[a]
                push!(capacities, old.capacity)
                push!(before_temperature, old.thermal / old.capacity)
            end
            push!(ions, (; before = old, after = actual, beta = beta[a], target = targets[a], materialized))
        end
        exact_pressure = _exact_rational(old_pressure) - gamma * change / volume
        exact_pressure >= 0 || throw(ArgumentError("electron-ion actual heat exceeds stored electron energy"))
        pressure = iszero(change) ? old_pressure : _from_exact(T, exact_pressure)
        any(active) && !all_cold && pressure <= 0 &&
            throw(ArgumentError("electron-ion target lost positive electron thermal energy"))
        radius = _electron_ion_radius(pressure)
        abs(_exact_rational(pressure) - exact_pressure) <= radius ||
            throw(ArgumentError("electron-ion pressure publication exceeded its rounding bound"))
        work.pressure[cell] = pressure
        energy_defect = change + volume * (_exact_rational(pressure) - _exact_rational(old_pressure)) / gamma
        energy_bound = volume * radius / gamma
        abs(energy_defect) <= energy_bound || throw(ArgumentError("electron-ion actual combined-energy certificate failed"))
        push!(after_temperature, _exact_rational(pressure) / _exact_rational(work.density[cell]))
        for a in eachindex(species)
            active[a] && push!(after_temperature, ions[a].after.thermal / statistics[a].capacity)
        end
        unchanged = all(a -> ions[a].before == ions[a].after, eachindex(ions)) && isequal(old_pressure, pressure)
        entropy, entropy_bounds = unchanged ? (:unchanged, (zero(Q), zero(Q))) :
            _electron_ion_entropy(capacities, before_temperature, after_temperature, config.log_bits)
        push!(cells, (; cell, density = work.density[cell], pressure_before = old_pressure,
            pressure_after = pressure, ions, electron_target, entropy, entropy_bounds,
            energy_defect, energy_bound,
            charge_defect = Ne - sum(a -> _exact_rational(species[a].q) * statistics[a].weight, eachindex(species))))
    end
    next_time, next_step = _validated_step_clock(inner.time[], inner.step[], h)
    work.pending = (policy = 1, time = next_time, step = next_step, dt = h, cells = cells)
    # Every target and certificate is prepared before publishing any velocity/pe.
    for a in eachindex(species), c = 1:3
        copyto!(species[a].v[c], work.bins.populations[a].velocity[c])
    end
    copyto!(inner.fields.pe, work.pressure)
    return nothing
end

const _ELECTRON_ION_REPORT_KEYS = (:policy, :time, :step, :dt, :cells)
const _ELECTRON_ION_CELL_KEYS = (
    :cell, :density, :pressure_before, :pressure_after, :ions, :electron_target,
    :entropy, :entropy_bounds, :energy_defect, :energy_bound, :charge_defect,
)
const _ELECTRON_ION_ION_KEYS = (:before, :after, :beta, :target, :materialized)
const _ELECTRON_ION_STATS_KEYS = (:weight, :mean, :kinetic, :momentum, :thermal, :capacity, :count)
const _ELECTRON_ION_MATERIALIZED_KEYS = (:scale, :seed_energy, :frame_energy, :momentum_bound, :kinetic_bound)

function _electron_ion_record(value, names)
    value isa NamedTuple && keys(value) == names ||
        throw(ArgumentError("electron-ion receipt has an invalid field inventory"))
    return nothing
end

_electron_ion_finite(value) = value isa _ElectronIonQ && denominator(value) > 0

function _electron_ion_validate_statistics(s, mass, current)
    Q = _ElectronIonQ
    _electron_ion_record(s, _ELECTRON_ION_STATS_KEYS)
    all(x -> _electron_ion_finite(x) && x >= 0, (s.weight, s.kinetic, s.thermal, s.capacity)) &&
        s.mean isa NTuple{3,Q} && all(_electron_ion_finite, s.mean) &&
        s.momentum isa NTuple{3,Q} && all(_electron_ion_finite, s.momentum) && s.count isa Int && s.count >= 0 ||
        throw(ArgumentError("electron-ion receipt statistics have invalid scalar types or signs"))
    s.weight == current.weight && s.count == current.count && s.capacity == 3s.weight / 2 &&
        all(c -> s.momentum[c] == mass * s.weight * s.mean[c], 1:3) &&
        s.kinetic == s.thermal + mass * s.weight * sum(x -> x * x, s.mean) / 2 ||
        throw(ArgumentError("electron-ion receipt statistics violate reservoir identities"))
    if iszero(s.weight)
        s == current || throw(ArgumentError("empty electron-ion receipt reservoir must be zero"))
    end
    return nothing
end

function _electron_ion_validate_materialized(row, ps, bins, cell, ::Type{T}) where {T}
    Q = _ElectronIonQ
    m, old, actual = row.materialized, row.before, row.after
    _electron_ion_record(m, _ELECTRON_ION_MATERIALIZED_KEYS)
    all(x -> _electron_ion_finite(x) && x >= 0, (m.scale, m.seed_energy, m.frame_energy, m.kinetic_bound)) &&
        m.momentum_bound isa NTuple{3,Q} && all(x -> _electron_ion_finite(x) && x >= 0, m.momentum_bound) &&
        _electron_ion_finite(row.target) && row.target >= 0 ||
        throw(ArgumentError("electron-ion receipt materialization has invalid scalar types or signs"))
    if row.target == old.thermal
        m == (scale = one(Q), seed_energy = old.thermal, frame_energy = old.kinetic,
            momentum_bound = (zero(Q), zero(Q), zero(Q)), kinetic_bound = zero(Q)) && actual == old ||
            throw(ArgumentError("zero-transfer electron-ion receipt must retain reservoir statistics"))
        return nothing
    end
    row.target > 0 && m.seed_energy > 0 && m.scale > 0 ||
        throw(ArgumentError("nontrivial electron-ion receipt requires positive target and scale"))
    if iszero(old.thermal)
        old.count >= 2 || throw(ArgumentError("cold electron-ion receipt needs two positive markers"))
    else
        m.seed_energy == old.thermal || throw(ArgumentError("warm electron-ion receipt must scale its original energy"))
    end
    m.scale == _electron_ion_root_scale(T, row.target / m.seed_energy, RoundNearest) ||
        throw(ArgumentError("electron-ion receipt scale differs from canonical publication"))
    u = BigInt(1) // (BigInt(1) << precision(T))
    abs(m.scale^2 * m.seed_energy - row.target) <= (2u + u^2) * row.target ||
        throw(ArgumentError("electron-ion receipt scale is outside its accuracy bound"))
    mass = _exact_rational(ps.m)
    m.frame_energy == mass * old.weight * sum(x -> x * x, old.mean) / 2 + m.scale^2 * m.seed_energy ||
        throw(ArgumentError("electron-ion receipt frame energy is inconsistent"))
    momentum_cap, kinetic_cap = _electron_ion_endpoint_bounds(ps, bins, cell)
    all(c -> m.momentum_bound[c] == momentum_cap[c] &&
        abs(actual.momentum[c] - old.momentum[c]) <= m.momentum_bound[c], 1:3) &&
        m.kinetic_bound <= kinetic_cap && abs(actual.kinetic - m.frame_energy) <= m.kinetic_bound ||
        throw(ArgumentError("electron-ion receipt exceeds actual endpoint rounding bounds"))
    return nothing
end

# Validation uses current endpoint particles and pressure, not a collision/RNG
# replay. Fresh derived bins/deposition are private and never published.
function _validate_electron_ion_report(report, config, inner, species, fields, time, step, dt, epsilon, number_scale)
    if step == 0
        report === nothing || throw(ArgumentError("unadvanced electron-ion state cannot carry a report"))
        return nothing
    end
    _electron_ion_record(report, _ELECTRON_ION_REPORT_KEYS)
    T = eltype(fields.pe)
    report.policy isa Int && report.policy == 1 && report.time isa T &&
        isequal(report.time, time) && report.step isa Int && report.step == step &&
        report.dt isa T && isequal(report.dt, dt) && report.cells isa Vector{NamedTuple} ||
        throw(ArgumentError("electron-ion receipt does not match its completed source clock"))
    if !_electron_ion_active(config)
        isempty(report.cells) || throw(ArgumentError("inactive electron-ion state must have an empty receipt"))
        return nothing
    end
    length(report.cells) == length(fields.pe) || throw(DimensionMismatch("electron-ion receipt cell count differs"))
    work = _electron_ion_work(config, inner, species)
    _electron_ion_density!(work, species, inner)
    for a in eachindex(species)
        ps = species[a]
        _bin_coulomb_particles!(work.bins.populations[a], ps, p -> _electron_ion_node(ps, inner.g, p))
    end
    Q = _ElectronIonQ
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    pi_bounds = _electron_ion_pi_bounds(T)
    for cell in eachindex(fields.pe)
        receipt = report.cells[cell]
        _electron_ion_record(receipt, _ELECTRON_ION_CELL_KEYS)
        receipt.cell isa Int && receipt.cell == cell && receipt.density isa T &&
            isequal(receipt.density, work.density[cell]) && receipt.pressure_before isa T &&
            isfinite(receipt.pressure_before) && receipt.pressure_before >= 0 &&
            receipt.pressure_after isa T && isequal(receipt.pressure_after, fields.pe[cell]) &&
            receipt.ions isa Vector{NamedTuple} && length(receipt.ions) == length(species) &&
            _electron_ion_finite(receipt.electron_target) && receipt.electron_target >= 0 &&
            all(_electron_ion_finite, (receipt.energy_defect, receipt.energy_bound, receipt.charge_defect)) &&
            receipt.entropy_bounds isa NTuple{2,Q} && all(_electron_ion_finite, receipt.entropy_bounds) ||
            throw(ArgumentError("electron-ion cell receipt has invalid endpoint values"))
        old = NamedTuple[]
        for a in eachindex(species)
            row, ps = receipt.ions[a], species[a]
            _electron_ion_record(row, _ELECTRON_ION_ION_KEYS)
            current = _electron_ion_stats(ps, work.bins.populations[a], cell)
            _electron_ion_validate_statistics(row.before, _exact_rational(ps.m), current)
            _electron_ion_validate_statistics(row.after, _exact_rational(ps.m), current)
            row.after == current && row.beta isa NTuple{3,Q} && all(_electron_ion_finite, row.beta) ||
                throw(ArgumentError("electron-ion receipt differs from stored endpoint particles"))
            _electron_ion_validate_materialized(row, ps, work.bins.populations[a], cell, T)
            push!(old, row.before)
        end
        Ne = _exact_rational(receipt.density) * volume
        Ce, Te = Ne / gamma, _exact_rational(receipt.pressure_before) / _exact_rational(receipt.density)
        active, beta, all_cold = _electron_ion_coefficients(config, old, species, Ce, Te, volume, dt, epsilon, number_scale, pi_bounds, T)
        any(active) && !all_cold && receipt.pressure_after <= 0 &&
            throw(ArgumentError("electron-ion receipt lost positive electron thermal energy"))
        electron_target, targets = _electron_ion_targets(Ce, Te, old, beta)
        receipt.electron_target == electron_target &&
            all(a -> receipt.ions[a].beta == beta[a] && receipt.ions[a].target == targets[a], eachindex(species)) ||
            throw(ArgumentError("electron-ion receipt differs from its physical coefficient/target map"))
        change = sum((row.after.kinetic - row.before.kinetic for row in receipt.ions); init = zero(Q))
        exact_pressure = _exact_rational(receipt.pressure_before) - gamma * change / volume
        expected_pressure = iszero(change) ? receipt.pressure_before : _from_exact(T, exact_pressure)
        exact_pressure >= 0 && isequal(expected_pressure, receipt.pressure_after) ||
            throw(ArgumentError("electron-ion receipt does not pay actual ion kinetic change"))
        defect = change + volume * (_exact_rational(receipt.pressure_after) - _exact_rational(receipt.pressure_before)) / gamma
        bound = volume * _electron_ion_radius(receipt.pressure_after) / gamma
        receipt.energy_defect == defect && receipt.energy_bound == bound && abs(defect) <= bound &&
            receipt.charge_defect == Ne - sum(a -> _exact_rational(species[a].q) * old[a].weight, eachindex(species)) ||
            throw(ArgumentError("electron-ion receipt energy or charge identity is inconsistent"))
        capacities, before, after = Q[Ce], Q[Te], Q[_exact_rational(receipt.pressure_after) / _exact_rational(receipt.density)]
        for a in eachindex(species)
            active[a] || continue
            !all_cold && receipt.ions[a].after.thermal <= 0 &&
                throw(ArgumentError("electron-ion receipt lost positive ion thermal energy"))
            push!(capacities, old[a].capacity)
            push!(before, old[a].thermal / old[a].capacity)
            push!(after, receipt.ions[a].after.thermal / old[a].capacity)
        end
        unchanged = all(row -> row.before == row.after, receipt.ions) && isequal(receipt.pressure_before, receipt.pressure_after)
        entropy, entropy_bounds = unchanged ? (:unchanged, (zero(Q), zero(Q))) :
            _electron_ion_entropy(capacities, before, after, config.log_bits)
        receipt.entropy === entropy && receipt.entropy_bounds == entropy_bounds ||
            throw(ArgumentError("electron-ion actual entropy receipt is inconsistent"))
    end
    return nothing
end
