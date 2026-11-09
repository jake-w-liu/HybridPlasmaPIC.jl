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

struct _ElectronIonThermalConfig{T,N,S}
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
function _bind_electron_ion(
    request::ElectronIonThermalExchange,
    inner::HybridStepper{D,T},
    species,
) where {D,T}
    _evolves_pressure(inner.model.closure) &&
        inner.shape isa Union{NGP,CIC,TSC} &&
        iszero(inner.model.de2) || throw(
        ArgumentError(
            "electron-ion thermal exchange requires stored scalar pressure, NGP/CIC/TSC and de2 = 0",
        ),
    )
    T <: Union{Float32,Float64} ||
        throw(ArgumentError("electron-ion thermal exchange requires Float32 or Float64"))
    T(inner.model.closure.γ) == T(5 / 3) ||
        throw(ArgumentError("electron-ion thermal exchange requires grid gamma = T(5/3)"))
    mass = _require_finite_positive_real("electron_mass", request.electron_mass, T)
    input = request.coulomb_log
    input isa Real ||
        length(input) == length(species) ||
        throw(
            DimensionMismatch("electron-ion Coulomb logarithms require one value per ion species"),
        )
    logs = T[_electron_ion_log(input isa Real ? input : input[a], T) for a in eachindex(species)]
    bits = request.log_bits === nothing ? precision(T) + 16 : request.log_bits
    precision(T) + 16 <= bits <= 4096 ||
        throw(ArgumentError("log_bits must be between precision(T)+16 and 4096"))
    return _ElectronIonThermalConfig{T,length(logs),typeof(inner.shape)}(
        mass,
        logs,
        Tuple(logs),
        bits,
        _electron_ion_policy(inner.shape),
    )
end
_bind_electron_ion(request, inner, species) =
    throw(ArgumentError("electron_ion must be nothing or ElectronIonThermalExchange"))

function _electron_ion_log(value, ::Type{T}) where {T}
    value isa Real && isfinite(value) && value >= 0 ||
        throw(ArgumentError("electron-ion Coulomb logarithms must be finite and nonnegative"))
    result = _require_finite_nonnegative_real("electron-ion Coulomb logarithm", value, T)
    !iszero(value) &&
        iszero(result) &&
        throw(ArgumentError("nonzero electron-ion Coulomb logarithms must survive grid conversion"))
    return result
end

_electron_ion_work(::Nothing, inner, species) = nothing
function _electron_ion_work(
    ::_ElectronIonThermalConfig,
    inner::HybridStepper{D,T},
    species,
) where {D,T}
    return _ElectronIonThermalWork(
        CoulombCollisionWorkspace(species, inner.g),
        zeros(T, inner.g.n),
        zeros(T, inner.g.n),
        zeros(T, inner.g.n),
        nothing,
    )
end

_electron_ion_active(::Nothing) = false
_electron_ion_active(config::_ElectronIonThermalConfig) = any(!iszero, config.logarithm)
_electron_ion_configuration(::Nothing) = nothing
_electron_ion_configuration(config::_ElectronIonThermalConfig{T,N,NGP}) where {T,N} = (
    policy = config.policy,
    electron_mass = config.electron_mass,
    logarithm = config.binding,
    log_bits = config.log_bits,
)

_electron_ion_policy(::NGP) = 1
_electron_ion_policy(::Union{CIC,TSC}) = 2
_electron_ion_shape(::Type{CIC}) = :cic
_electron_ion_shape(::Type{TSC}) = :tsc
function _electron_ion_configuration(config::_ElectronIonThermalConfig{T,N,S}) where {T,N,S}
    return (
        policy = config.policy,
        electron_mass = config.electron_mass,
        logarithm = config.binding,
        log_bits = config.log_bits,
        shape = _electron_ion_shape(S),
        frame = :periodic_floor_cell,
        incidence = :exact_stored_shape_products,
        order = :species_frame_node,
        publication = :current_frame_causal_node,
    )
end

function _validate_electron_ion_config(
    config::_ElectronIonThermalConfig{T,N,S},
    inner,
    species,
) where {T,N,S}
    inner.shape isa S &&
        config.policy == _electron_ion_policy(inner.shape) &&
        isfinite(config.electron_mass) &&
        config.electron_mass > 0 ||
        throw(ArgumentError("invalid electron-ion thermal policy or electron mass"))
    length(config.logarithm) == length(species) == length(config.binding) && all(
        a ->
            isequal(config.logarithm[a], config.binding[a]) &&
                isfinite(config.logarithm[a]) &&
                config.logarithm[a] >= 0,
        eachindex(config.logarithm),
    ) || throw(ArgumentError("electron-ion Coulomb logarithms are fixed"))
    precision(T) + 16 <= config.log_bits <= 4096 ||
        throw(ArgumentError("invalid electron-ion entropy precision"))
    _evolves_pressure(inner.model.closure) &&
        inner.shape isa Union{NGP,CIC,TSC} &&
        T(inner.model.closure.γ) == T(5 / 3) &&
        iszero(inner.model.de2) ||
        throw(ArgumentError("electron-ion thermal model and shape are fixed"))
    return nothing
end

# Match the canonical nearest-node stencil, including its ties-up and periodic
# policies. The existing ion-ion floor-cell bins deliberately stay unchanged.
function _electron_ion_node(ps::ParticleSet{D}, g, p) where {D}
    index =
        ntuple(d -> mod(_stencil1d(NGP(), _particle_cell_position(ps, g, d, p))[1], g.n[d]) + 1, D)
    return LinearIndices(g.n)[index...]
end

function _electron_ion_density!(work, species, inner)
    fill!(work.density, zero(eltype(work.density)))
    for ps in species
        density!(work.number_density, ps, inner.g, inner.shape)
        @. work.density =
            getproperty(_charge_moment_parts(work.density, ps.q, work.number_density), :value)
    end
    _validate_transport_density(
        work.density,
        eltype(work.density)(inner.model.nfloor),
        "electron-ion exchange",
    )
    return nothing
end

function _electron_ion_stats(ps, bins, cell, velocity = ps.v)
    Q = _ElectronIonQ
    weight = kinetic = zero(Q)
    momentum = (zero(Q), zero(Q), zero(Q))
    count = bins.offset[cell+1] - bins.offset[cell]
    mass = _exact_rational(ps.m)
    for k = bins.offset[cell]:(bins.offset[cell+1]-1)
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        weight += w
        v = ntuple(c -> _exact_rational(velocity[c][p]), 3)
        momentum = (
            momentum[1] + mass * w * v[1],
            momentum[2] + mass * w * v[2],
            momentum[3] + mass * w * v[3],
        )
        kinetic += mass * w * sum(x -> x * x, v) / 2
    end
    mean = iszero(weight) ? (zero(Q), zero(Q), zero(Q)) :
        (momentum[1] / (mass * weight), momentum[2] / (mass * weight), momentum[3] / (mass * weight))
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
    for n = 0:(cld(precision(T) + 18, 4)-1)
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
        (1 - beta[a][2]) * s.thermal + beta[a][2] * s.capacity * electron_temperature for
        (a, s) in enumerate(statistics)
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

function _electron_ion_publish!(bins, ps, cell, old, target, rng, source_velocity = ps.v)
    T = eltype(ps.weight)
    Q = _ElectronIonQ
    if target == old.thermal
        return (
            scale = one(Q),
            seed_energy = old.thermal,
            frame_energy = old.kinetic,
            momentum_bound = (zero(Q), zero(Q), zero(Q)),
            kinetic_bound = zero(Q),
        )
    end
    seed_mean, seed_energy = old.mean, old.thermal
    if iszero(seed_energy)
        old.count >= 2 ||
            throw(ArgumentError("a single positive marker cannot receive momentum-preserving heat"))
        for k = bins.offset[cell]:(bins.offset[cell+1]-1), c = 1:3
            p = bins.index[k]
            draw = randn(rng, T)
            isfinite(draw) || throw(ArgumentError("electron-ion cold seed must be finite"))
            bins.velocity[c][p] = draw
        end
        seed = _electron_ion_stats(ps, bins, cell, bins.velocity)
        seed_mean, seed_energy = seed.mean, seed.thermal
        seed_energy > 0 ||
            throw(ArgumentError("electron-ion cold seed has zero weighted thermal energy"))
    end
    scale = _electron_ion_root_scale(T, target / seed_energy, RoundNearest)
    u = BigInt(1) // (BigInt(1) << precision(T))
    abs(scale^2 * seed_energy - target) <= (2u + u^2) * target ||
        throw(ArgumentError("electron-ion scale accuracy is unresolved"))
    momentum_bound = (zero(Q), zero(Q), zero(Q))
    kinetic_bound = zero(Q)
    mass = _exact_rational(ps.m)
    for k = bins.offset[cell]:(bins.offset[cell+1]-1)
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        for c = 1:3
            seed = iszero(old.thermal) ? bins.velocity[c][p] : source_velocity[c][p]
            frame = old.mean[c] + scale * (_exact_rational(seed) - seed_mean[c])
            value = _from_exact(T, frame)
            radius = _electron_ion_radius(value)
            abs(_exact_rational(value) - frame) <= radius || throw(
                ArgumentError("electron-ion component publication exceeded its rounding bound"),
            )
            bins.velocity[c][p] = value
            momentum_bound = (
                momentum_bound[1] + (c == 1 ? mass * w * radius : zero(Q)),
                momentum_bound[2] + (c == 2 ? mass * w * radius : zero(Q)),
                momentum_bound[3] + (c == 3 ? mass * w * radius : zero(Q)),
            )
            kinetic_bound += mass * w * (abs(frame) * radius + radius^2 / 2)
        end
    end
    frame_energy = mass * old.weight * sum(x -> x * x, old.mean) / 2 + scale^2 * seed_energy
    return (; scale, seed_energy, frame_energy, momentum_bound, kinetic_bound)
end

function _electron_ion_entropy(capacities, before, after, bits)
    Q = _ElectronIonQ
    all(>(zero(Q)), after) ||
        throw(ArgumentError("electron-ion active reservoirs must publish positive energy"))
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
    low >= 0 ||
        throw(ArgumentError("electron-ion actual entropy decrease or unresolved entropy sign"))
    return :resolved, (low, high)
end

# Bounds used when checking an incoming receipt without re-running the source.
# |frame| <= |stored|+rho caps the saved, tighter frame-based energy bound.
function _electron_ion_endpoint_bounds(ps, bins, cell)
    Q = _ElectronIonQ
    momentum = (zero(Q), zero(Q), zero(Q))
    kinetic = zero(Q)
    mass = _exact_rational(ps.m)
    for k = bins.offset[cell]:(bins.offset[cell+1]-1), c = 1:3
        p = bins.index[k]
        w = _exact_rational(ps.weight[p])
        radius = _electron_ion_radius(ps.v[c][p])
        momentum = (
            momentum[1] + (c == 1 ? mass * w * radius : zero(Q)),
            momentum[2] + (c == 2 ? mass * w * radius : zero(Q)),
            momentum[3] + (c == 3 ? mass * w * radius : zero(Q)),
        )
        kinetic += mass * w * (abs(_exact_rational(ps.v[c][p])) * radius + 3radius^2 / 2)
    end
    return momentum, kinetic
end

function _electron_ion_coefficients(
    config,
    statistics,
    species,
    Ce,
    Te,
    volume,
    h,
    epsilon,
    number_scale,
    pi_bounds,
    ::Type{T},
    incidence_weights = nothing,
) where {T}
    Q = _ElectronIonQ
    active =
        [!iszero(config.logarithm[a]) && !iszero(statistics[a].weight) for a in eachindex(species)]
    all_cold =
        iszero(Te) && all(a -> !active[a] || iszero(statistics[a].thermal), eachindex(species))
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
        K =
            _exact_rational(h) *
            (Ce / s.capacity) *
            _exact_rational(number_scale) *
            ((incidence_weights === nothing ? s.weight : incidence_weights[a]) / volume) *
            charge^2 *
            _exact_rational(config.logarithm[a]) / (_exact_rational(epsilon)^2 * me * mass)
        push!(beta, _electron_ion_beta(T, K, S, pi_bounds))
    end
    return active, beta, all_cold
end

const _ELECTRON_ION_OVERLAP_KEYS = (
    :policy,
    :time,
    :step,
    :dt,
    :shape,
    :frame_policy,
    :edge_order,
    :density,
    :initial_pressure,
    :initial_frames,
    :incidence,
    :edges,
    :frame_partition,
    :charge_defect,
    :energy_defect,
    :energy_bound,
    :momentum_defect,
    :momentum_bound,
)
const _ELECTRON_ION_OVERLAP_FRAME_KEYS = (:species, :frame, :statistics)
const _ELECTRON_ION_OVERLAP_INCIDENCE_KEYS = (:species, :frame, :node, :weight)
const _ELECTRON_ION_OVERLAP_EDGE_KEYS = (
    :species,
    :frame,
    :node,
    :weight,
    :beta,
    :cold,
    :pressure_before,
    :pressure_after,
    :before,
    :after,
    :target,
    :electron_target,
    :materialized,
    :entropy,
    :entropy_bounds,
    :energy_defect,
    :energy_bound,
)

function _electron_ion_overlap_empty_report(
    config::_ElectronIonThermalConfig{T,N,S},
    time,
    step,
    dt,
) where {T,N,S}
    Q = _ElectronIonQ
    return (
        policy = 2,
        time = time,
        step = step,
        dt = dt,
        shape = _electron_ion_shape(S),
        frame_policy = :periodic_floor_cell,
        edge_order = :species_frame_node,
        density = nothing,
        initial_pressure = nothing,
        initial_frames = NamedTuple[],
        incidence = NamedTuple[],
        edges = NamedTuple[],
        frame_partition = Q[],
        charge_defect = Q[],
        energy_defect = zero(Q),
        energy_bound = zero(Q),
        momentum_defect = (zero(Q), zero(Q), zero(Q)),
        momentum_bound = (zero(Q), zero(Q), zero(Q)),
    )
end

function _electron_ion_overlap_layout!(work, species, inner)
    _electron_ion_density!(work, species, inner)
    strides = _particle_cell_strides(inner.g)
    frames = NamedTuple[]
    incidence = Dict{NTuple{3,Int},_ElectronIonQ}()
    ncells = work.bins.ncells
    for a in eachindex(species)
        ps, bins = species[a], work.bins.populations[a]
        _bin_coulomb_particles!(bins, ps, inner.g, strides)
        for frame = 1:ncells
            push!(
                frames,
                (
                    species = a,
                    frame = frame,
                    statistics = _electron_ion_stats(ps, bins, frame, bins.velocity),
                ),
            )
        end
        _foreach_exact_particle_contribution(
            ps,
            (ps.weight,),
            inner.g,
            inner.shape,
            Val(:periodic),
        ) do index, weight, p
            _electron_ion_finite(weight) && weight >= 0 || throw(
                ArgumentError("thermal overlap requires finite nonnegative stored shape factors"),
            )
            iszero(weight) && return nothing
            node, frame = LinearIndices(inner.g.n)[index...], bins.cell[p]
            key = (a, frame, node)
            incidence[key] = get(incidence, key, zero(_ElectronIonQ)) + weight
            return nothing
        end
    end
    entries = NamedTuple[]
    for (a, frame, node) in sort!(collect(keys(incidence)))
        push!(
            entries,
            (species = a, frame = frame, node = node, weight = incidence[(a, frame, node)]),
        )
    end
    partition = _ElectronIonQ[-row.statistics.weight for row in frames]
    volume = prod(_exact_rational, inner.g.dx)
    defect = map(x -> _exact_rational(x) * volume, work.density)
    for edge in entries
        f = (edge.species - 1) * ncells + edge.frame
        partition[f] += edge.weight
        defect[edge.node] -= _exact_rational(species[edge.species].q) * edge.weight
    end
    return frames, entries, partition, defect
end

# Components identify only the declared entirely-cold identity. They do not
# introduce a solve, field scheduler, cutoff or omission of a mixed cold edge.
function _electron_ion_overlap_plans(
    config,
    frames,
    incidence,
    density,
    pressure,
    species,
    inner,
    h,
    epsilon,
    number_scale,
)
    ncells = length(density)
    parent = collect(1:Base.checked_add(ncells, length(frames)))
    rank = zeros(Int, length(parent))
    function root(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        return vertex
    end
    for edge in incidence
        iszero(config.logarithm[edge.species]) && continue
        f = (edge.species - 1) * ncells + edge.frame
        left, right = root(ncells + f), root(edge.node)
        left == right && continue
        if rank[left] < rank[right]
            parent[left] = right
        else
            parent[right] = left
            rank[left] == rank[right] && (rank[left] += 1)
        end
    end
    hot = falses(length(parent))
    for j in eachindex(pressure)
        isfinite(pressure[j]) && pressure[j] >= 0 ||
            throw(ArgumentError("thermal overlap pressure must be finite and nonnegative"))
        pressure[j] > 0 && (hot[root(j)] = true)
    end
    for f in eachindex(frames)
        frames[f].statistics.thermal > 0 && (hot[root(ncells + f)] = true)
    end
    T = eltype(density)
    Q = _ElectronIonQ
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    pi_bounds = _electron_ion_pi_bounds(T)
    plans = NamedTuple[]
    for edge in incidence
        a, j = edge.species, edge.node
        iszero(config.logarithm[a]) && continue
        f = (a - 1) * ncells + edge.frame
        old = frames[f].statistics
        Ce = _exact_rational(density[j]) * volume / gamma
        Te = _exact_rational(pressure[j]) / _exact_rational(density[j])
        cold = !hot[root(j)]
        if cold
            beta = (zero(Q), zero(Q), zero(Q))
        else
            (iszero(Te) && iszero(old.thermal)) &&
                throw(ArgumentError("singular mixed cold thermal overlap edge"))
            iszero(old.thermal) &&
                old.count < 2 &&
                throw(
                    ArgumentError(
                        "a single positive marker cannot receive momentum-preserving heat",
                    ),
                )
            local_config =
                (electron_mass = config.electron_mass, logarithm = (config.logarithm[a],))
            _, coefficients, _ = _electron_ion_coefficients(
                local_config,
                (old,),
                (species[a],),
                Ce,
                Te,
                volume,
                h,
                epsilon,
                number_scale,
                pi_bounds,
                T,
                (edge.weight,),
            )
            beta = coefficients[1]
        end
        push!(plans, merge(edge, (; beta, cold, Ce)))
    end
    return plans
end

function _electron_ion_overlap_entropy(
    before,
    after,
    Ce,
    density,
    pressure_before,
    pressure_after,
    bits,
)
    Q = _ElectronIonQ
    unchanged = before == after && isequal(pressure_before, pressure_after)
    unchanged && return :unchanged, (zero(Q), zero(Q))
    return _electron_ion_entropy(
        (Ce, before.capacity),
        (
            _exact_rational(pressure_before) / _exact_rational(density),
            before.thermal / before.capacity,
        ),
        (
            _exact_rational(pressure_after) / _exact_rational(density),
            after.thermal / before.capacity,
        ),
        bits,
    )
end

function _prepare_electron_ion_overlap!(config, work, inner, species, h, epsilon, number_scale, rng)
    T = eltype(work.pressure)
    Q = _ElectronIonQ
    frames, incidence, partition, charge_defect =
        _electron_ion_overlap_layout!(work, species, inner)
    copyto!(work.pressure, inner.fields.pe)
    plans = _electron_ion_overlap_plans(
        config,
        frames,
        incidence,
        work.density,
        work.pressure,
        species,
        inner,
        h,
        epsilon,
        number_scale,
    )
    initial_pressure = copy(work.pressure)
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    edges = NamedTuple[]
    energy_defect = energy_bound = zero(Q)
    momentum_defect = momentum_bound = (zero(Q), zero(Q), zero(Q))
    for plan in plans
        a, f, j = plan.species, plan.frame, plan.node
        ps, bins = species[a], work.bins.populations[a]
        before = _electron_ion_stats(ps, bins, f, bins.velocity)
        pressure_before = work.pressure[j]
        Te = _exact_rational(pressure_before) / _exact_rational(work.density[j])
        electron_target, targets = _electron_ion_targets(plan.Ce, Te, (before,), (plan.beta,))
        target = targets[1]
        materialized = _electron_ion_publish!(bins, ps, f, before, target, rng, bins.velocity)
        after = _electron_ion_stats(ps, bins, f, bins.velocity)
        all(
            c -> abs(after.momentum[c] - before.momentum[c]) <= materialized.momentum_bound[c],
            1:3,
        ) && abs(after.kinetic - materialized.frame_energy) <= materialized.kinetic_bound ||
            throw(ArgumentError("thermal overlap native frame certificate failed"))
        change = after.kinetic - before.kinetic
        exact_pressure = _exact_rational(pressure_before) - gamma * change / volume
        exact_pressure >= 0 ||
            throw(ArgumentError("thermal overlap actual heat exceeds its causal electron energy"))
        pressure_after = iszero(change) ? pressure_before : _from_exact(T, exact_pressure)
        radius = _electron_ion_radius(pressure_after)
        abs(_exact_rational(pressure_after) - exact_pressure) <= radius ||
            throw(ArgumentError("thermal overlap pressure exceeded its publication bound"))
        !plan.cold &&
            (pressure_after <= 0 || after.thermal <= 0) &&
            throw(ArgumentError("thermal overlap target lost positive thermal energy"))
        defect =
            change +
            volume * (_exact_rational(pressure_after) - _exact_rational(pressure_before)) / gamma
        bound = volume * radius / gamma
        abs(defect) <= bound ||
            throw(ArgumentError("thermal overlap causal-node energy certificate failed"))
        entropy, entropy_bounds = _electron_ion_overlap_entropy(
            before,
            after,
            plan.Ce,
            work.density[j],
            pressure_before,
            pressure_after,
            config.log_bits,
        )
        work.pressure[j] = pressure_after
        energy_defect += defect
        energy_bound += bound
        momentum_defect = (
            momentum_defect[1] + after.momentum[1] - before.momentum[1],
            momentum_defect[2] + after.momentum[2] - before.momentum[2],
            momentum_defect[3] + after.momentum[3] - before.momentum[3],
        )
        momentum_bound = (
            momentum_bound[1] + materialized.momentum_bound[1],
            momentum_bound[2] + materialized.momentum_bound[2],
            momentum_bound[3] + materialized.momentum_bound[3],
        )
        push!(
            edges,
            (
                species = a,
                frame = f,
                node = j,
                weight = plan.weight,
                beta = plan.beta,
                cold = plan.cold,
                pressure_before = pressure_before,
                pressure_after = pressure_after,
                before = before,
                after = after,
                target = target,
                electron_target = electron_target,
                materialized = materialized,
                entropy = entropy,
                entropy_bounds = entropy_bounds,
                energy_defect = defect,
                energy_bound = bound,
            ),
        )
    end
    next_time, next_step = _validated_step_clock(inner.time[], inner.step[], h)
    report = _electron_ion_overlap_empty_report(config, next_time, next_step, h)
    work.pending = merge(
        report,
        (
            density = copy(work.density),
            initial_pressure = initial_pressure,
            initial_frames = frames,
            incidence = incidence,
            edges = edges,
            frame_partition = partition,
            charge_defect = charge_defect,
            energy_defect = energy_defect,
            energy_bound = energy_bound,
            momentum_defect = momentum_defect,
            momentum_bound = momentum_bound,
        ),
    )
    for a in eachindex(species), c = 1:3
        copyto!(species[a].v[c], work.bins.populations[a].velocity[c])
    end
    copyto!(inner.fields.pe, work.pressure)
    return nothing
end

function _prepare_electron_ion!(config, work, inner, species, h, epsilon, number_scale, rng)
    T = eltype(work.pressure)
    Q = _ElectronIonQ
    _electron_ion_density!(work, species, inner)
    for a in eachindex(species)
        ps = species[a]
        _bin_coulomb_particles!(
            work.bins.populations[a],
            ps,
            p -> _electron_ion_node(ps, inner.g, p),
        )
    end
    copyto!(work.pressure, inner.fields.pe)
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    pi_bounds = _electron_ion_pi_bounds(T)
    plans = NamedTuple[]
    for cell in eachindex(work.pressure)
        statistics = [
            _electron_ion_stats(species[a], work.bins.populations[a], cell) for
            a in eachindex(species)
        ]
        Ne = _exact_rational(work.density[cell]) * volume
        Ce = Ne / gamma
        old_pressure = inner.fields.pe[cell]
        isfinite(old_pressure) && old_pressure >= 0 ||
            throw(ArgumentError("electron-ion pressure must be finite and nonnegative"))
        Te = _exact_rational(old_pressure) / _exact_rational(work.density[cell])
        active, beta, all_cold = _electron_ion_coefficients(
            config,
            statistics,
            species,
            Ce,
            Te,
            volume,
            h,
            epsilon,
            number_scale,
            pi_bounds,
            T,
        )
        electron_target, targets = _electron_ion_targets(Ce, Te, statistics, beta)
        # Every cell's rate and sparse checks precede every thermal draw.
        for a in eachindex(species)
            iszero(statistics[a].thermal) &&
                targets[a] > 0 &&
                statistics[a].count < 2 &&
                throw(
                    ArgumentError(
                        "a single positive marker cannot receive momentum-preserving heat",
                    ),
                )
        end
        push!(
            plans,
            (;
                cell,
                statistics,
                Ce,
                Te,
                old_pressure,
                active,
                beta,
                all_cold,
                electron_target,
                targets,
            ),
        )
    end
    cells = NamedTuple[]
    for plan in plans
        (;
            cell,
            statistics,
            Ce,
            Te,
            old_pressure,
            active,
            beta,
            all_cold,
            electron_target,
            targets,
        ) = plan
        Ne = _exact_rational(work.density[cell]) * volume
        ions = NamedTuple[]
        change = zero(Q)
        before_temperature, after_temperature = Q[Te], Q[]
        capacities = Q[Ce]
        for a in eachindex(species)
            ps, bins, old = species[a], work.bins.populations[a], statistics[a]
            materialized = _electron_ion_publish!(bins, ps, cell, old, targets[a], rng)
            actual = _electron_ion_stats(ps, bins, cell, bins.velocity)
            all(
                c -> abs(actual.momentum[c] - old.momentum[c]) <= materialized.momentum_bound[c],
                1:3,
            ) && abs(actual.kinetic - materialized.frame_energy) <= materialized.kinetic_bound ||
                throw(
                    ArgumentError(
                        "electron-ion actual momentum or kinetic-energy rounding certificate failed",
                    ),
                )
            active[a] &&
                !all_cold &&
                actual.thermal <= 0 &&
                throw(ArgumentError("electron-ion target lost positive ion thermal energy"))
            change += actual.kinetic - old.kinetic
            if active[a]
                push!(capacities, old.capacity)
                push!(before_temperature, old.thermal / old.capacity)
            end
            push!(
                ions,
                (; before = old, after = actual, beta = beta[a], target = targets[a], materialized),
            )
        end
        exact_pressure = _exact_rational(old_pressure) - gamma * change / volume
        exact_pressure >= 0 ||
            throw(ArgumentError("electron-ion actual heat exceeds stored electron energy"))
        pressure = iszero(change) ? old_pressure : _from_exact(T, exact_pressure)
        any(active) &&
            !all_cold &&
            pressure <= 0 &&
            throw(ArgumentError("electron-ion target lost positive electron thermal energy"))
        radius = _electron_ion_radius(pressure)
        abs(_exact_rational(pressure) - exact_pressure) <= radius ||
            throw(ArgumentError("electron-ion pressure publication exceeded its rounding bound"))
        work.pressure[cell] = pressure
        energy_defect =
            change + volume * (_exact_rational(pressure) - _exact_rational(old_pressure)) / gamma
        energy_bound = volume * radius / gamma
        abs(energy_defect) <= energy_bound ||
            throw(ArgumentError("electron-ion actual combined-energy certificate failed"))
        push!(after_temperature, _exact_rational(pressure) / _exact_rational(work.density[cell]))
        for a in eachindex(species)
            active[a] && push!(after_temperature, ions[a].after.thermal / statistics[a].capacity)
        end
        unchanged =
            all(a -> ions[a].before == ions[a].after, eachindex(ions)) &&
            isequal(old_pressure, pressure)
        entropy, entropy_bounds =
            unchanged ? (:unchanged, (zero(Q), zero(Q))) :
            _electron_ion_entropy(
                capacities,
                before_temperature,
                after_temperature,
                config.log_bits,
            )
        push!(
            cells,
            (;
                cell,
                density = work.density[cell],
                pressure_before = old_pressure,
                pressure_after = pressure,
                ions,
                electron_target,
                entropy,
                entropy_bounds,
                energy_defect,
                energy_bound,
                charge_defect = Ne - sum(
                    a -> _exact_rational(species[a].q) * statistics[a].weight,
                    eachindex(species),
                ),
            ),
        )
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
    :cell,
    :density,
    :pressure_before,
    :pressure_after,
    :ions,
    :electron_target,
    :entropy,
    :entropy_bounds,
    :energy_defect,
    :energy_bound,
    :charge_defect,
)
const _ELECTRON_ION_ION_KEYS = (:before, :after, :beta, :target, :materialized)
const _ELECTRON_ION_STATS_KEYS = (:weight, :mean, :kinetic, :momentum, :thermal, :capacity, :count)
const _ELECTRON_ION_MATERIALIZED_KEYS =
    (:scale, :seed_energy, :frame_energy, :momentum_bound, :kinetic_bound)

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
        s.mean isa NTuple{3,Q} &&
        all(_electron_ion_finite, s.mean) &&
        s.momentum isa NTuple{3,Q} &&
        all(_electron_ion_finite, s.momentum) &&
        s.count isa Int &&
        s.count >= 0 ||
        throw(ArgumentError("electron-ion receipt statistics have invalid scalar types or signs"))
    s.weight == current.weight &&
        s.count == current.count &&
        s.capacity == 3s.weight / 2 &&
        all(c -> s.momentum[c] == mass * s.weight * s.mean[c], 1:3) &&
        s.kinetic == s.thermal + mass * s.weight * sum(x -> x * x, s.mean) / 2 ||
        throw(ArgumentError("electron-ion receipt statistics violate reservoir identities"))
    if iszero(s.weight)
        s == current || throw(ArgumentError("empty electron-ion receipt reservoir must be zero"))
    end
    return nothing
end

function _electron_ion_validate_materialized(
    row,
    ps,
    bins,
    cell,
    ::Type{T},
    caps = nothing,
) where {T}
    Q = _ElectronIonQ
    m, old, actual = row.materialized, row.before, row.after
    _electron_ion_record(m, _ELECTRON_ION_MATERIALIZED_KEYS)
    all(
            x -> _electron_ion_finite(x) && x >= 0,
            (m.scale, m.seed_energy, m.frame_energy, m.kinetic_bound),
        ) &&
        m.momentum_bound isa NTuple{3,Q} &&
        all(x -> _electron_ion_finite(x) && x >= 0, m.momentum_bound) &&
        _electron_ion_finite(row.target) &&
        row.target >= 0 || throw(
        ArgumentError("electron-ion receipt materialization has invalid scalar types or signs"),
    )
    if row.target == old.thermal
        m == (
            scale = one(Q),
            seed_energy = old.thermal,
            frame_energy = old.kinetic,
            momentum_bound = (zero(Q), zero(Q), zero(Q)),
            kinetic_bound = zero(Q),
        ) && actual == old || throw(
            ArgumentError("zero-transfer electron-ion receipt must retain reservoir statistics"),
        )
        return nothing
    end
    row.target > 0 && m.seed_energy > 0 && m.scale > 0 ||
        throw(ArgumentError("nontrivial electron-ion receipt requires positive target and scale"))
    if iszero(old.thermal)
        old.count >= 2 ||
            throw(ArgumentError("cold electron-ion receipt needs two positive markers"))
    else
        m.seed_energy == old.thermal ||
            throw(ArgumentError("warm electron-ion receipt must scale its original energy"))
    end
    m.scale == _electron_ion_root_scale(T, row.target / m.seed_energy, RoundNearest) ||
        throw(ArgumentError("electron-ion receipt scale differs from canonical publication"))
    u = BigInt(1) // (BigInt(1) << precision(T))
    abs(m.scale^2 * m.seed_energy - row.target) <= (2u + u^2) * row.target ||
        throw(ArgumentError("electron-ion receipt scale is outside its accuracy bound"))
    mass = _exact_rational(ps.m)
    m.frame_energy ==
    mass * old.weight * sum(x -> x * x, old.mean) / 2 + m.scale^2 * m.seed_energy ||
        throw(ArgumentError("electron-ion receipt frame energy is inconsistent"))
    momentum_cap, kinetic_cap =
        caps === nothing ? _electron_ion_endpoint_bounds(ps, bins, cell) : caps
    all(
            c ->
                (
                    caps === nothing ? m.momentum_bound[c] == momentum_cap[c] :
                    m.momentum_bound[c] <= momentum_cap[c]
                ) && abs(actual.momentum[c] - old.momentum[c]) <= m.momentum_bound[c],
            1:3,
        ) &&
        m.kinetic_bound <= kinetic_cap &&
        abs(actual.kinetic - m.frame_energy) <= m.kinetic_bound ||
        throw(ArgumentError("electron-ion receipt exceeds actual endpoint rounding bounds"))
    return nothing
end

# Historical intermediate frames have subsequently been rescaled. These caps
# follow rho(x)<=u*abs(x)+eta and weighted Cauchy, using that edge's actual K.
# Evolution still checks the tighter component/frame bounds before each edge.
function _electron_ion_overlap_caps(statistics, mass, ::Type{T}) where {T}
    Q = _ElectronIonQ
    W, K = statistics.weight, statistics.kinetic
    iszero(W) && return (zero(Q), zero(Q), zero(Q)), zero(Q)
    S = 2K / mass
    u = BigInt(1) // (BigInt(1) << precision(T))
    eta = _exact_rational(nextfloat(zero(T)))
    first = iszero(K) ? zero(Q) : _electron_ion_root_scale(T, W * S, RoundUp)
    second = iszero(K) ? zero(Q) : _electron_ion_root_scale(T, 3W * S, RoundUp)
    momentum = mass * (u * first + eta * W)
    kinetic = mass * ((u + 3u^2 / 2) * S + (1 + 3u) * eta * second + 9eta^2 * W / 2)
    return (momentum, momentum, momentum), kinetic
end

function _validate_electron_ion_overlap_roles(report)
    report === nothing && return nothing
    _electron_ion_record(report, _ELECTRON_ION_OVERLAP_KEYS)
    roles = AbstractArray[]
    for name in (
        :density,
        :initial_pressure,
        :initial_frames,
        :incidence,
        :edges,
        :frame_partition,
        :charge_defect,
    )
        array = getproperty(report, name)
        array === nothing && name in (:density, :initial_pressure) && continue
        array isa Array ||
            throw(ArgumentError("thermal overlap receipt roles require standard Array storage"))
        for other in roles
            _particle_arrays_mightalias(array, other) &&
                throw(ArgumentError("thermal overlap receipt roles must be independent"))
        end
        push!(roles, array)
    end
    return nothing
end

function _validate_electron_ion_overlap_report(
    report,
    config::_ElectronIonThermalConfig{T,N,S},
    inner,
    species,
    fields,
    time,
    step,
    dt,
    epsilon,
    number_scale,
) where {T,N,S}
    if step == 0
        report === nothing ||
            throw(ArgumentError("unadvanced thermal overlap cannot carry a report"))
        return nothing
    end
    Q = _ElectronIonQ
    _electron_ion_record(report, _ELECTRON_ION_OVERLAP_KEYS)
    report.policy isa Int &&
        report.policy == 2 &&
        report.time isa T &&
        isequal(report.time, time) &&
        report.step isa Int &&
        report.step == step &&
        report.dt isa T &&
        isequal(report.dt, dt) &&
        report.shape === _electron_ion_shape(S) &&
        report.frame_policy === :periodic_floor_cell &&
        report.edge_order === :species_frame_node ||
        throw(ArgumentError("thermal overlap receipt does not match policy, shape or source clock"))
    _validate_electron_ion_overlap_roles(report)
    if !_electron_ion_active(config)
        expected = _electron_ion_overlap_empty_report(config, time, step, dt)
        typeof(report) === typeof(expected) && isequal(report, expected) ||
            throw(ArgumentError("inactive thermal overlap requires its empty numeric receipt"))
        return nothing
    end
    report.density isa Array{T} &&
        size(report.density) == inner.g.n &&
        report.initial_pressure isa Array{T} &&
        size(report.initial_pressure) == inner.g.n &&
        all(x -> isfinite(x) && x >= 0, report.initial_pressure) &&
        report.initial_frames isa Vector{NamedTuple} &&
        report.incidence isa Vector{NamedTuple} &&
        report.edges isa Vector{NamedTuple} &&
        report.frame_partition isa Vector{Q} &&
        report.charge_defect isa Array{Q} &&
        size(report.charge_defect) == inner.g.n &&
        all(_electron_ion_finite, report.frame_partition) &&
        all(_electron_ion_finite, report.charge_defect) &&
        _electron_ion_finite(report.energy_defect) &&
        _electron_ion_finite(report.energy_bound) &&
        report.energy_bound >= 0 &&
        report.momentum_defect isa NTuple{3,Q} &&
        all(_electron_ion_finite, report.momentum_defect) &&
        report.momentum_bound isa NTuple{3,Q} &&
        all(x -> _electron_ion_finite(x) && x >= 0, report.momentum_bound) ||
        throw(ArgumentError("thermal overlap receipt has invalid array/scalar types or signs"))
    work = _electron_ion_work(config, inner, species)
    actual, incidence, partition, charge_defect =
        _electron_ion_overlap_layout!(work, species, inner)
    isequal(report.density, work.density) &&
        report.incidence == incidence &&
        report.frame_partition == partition &&
        report.charge_defect == charge_defect &&
        length(report.initial_frames) == length(actual) || throw(
        ArgumentError(
            "thermal overlap receipt differs from actual shape incidence or source inventory",
        ),
    )
    ncells = work.bins.ncells
    for f in eachindex(actual)
        initial = report.initial_frames[f]
        _electron_ion_record(initial, _ELECTRON_ION_OVERLAP_FRAME_KEYS)
        initial.species isa Int &&
            initial.species == actual[f].species &&
            initial.frame isa Int &&
            initial.frame == actual[f].frame ||
            throw(ArgumentError("thermal overlap initial frame order is invalid"))
        _electron_ion_validate_statistics(
            initial.statistics,
            _exact_rational(species[initial.species].m),
            actual[f].statistics,
        )
    end
    for edge in report.incidence
        _electron_ion_record(edge, _ELECTRON_ION_OVERLAP_INCIDENCE_KEYS)
        edge.species isa Int &&
            edge.frame isa Int &&
            edge.node isa Int &&
            _electron_ion_finite(edge.weight) &&
            edge.weight > 0 ||
            throw(ArgumentError("thermal overlap incidence has invalid scalar types or signs"))
    end
    plans = _electron_ion_overlap_plans(
        config,
        report.initial_frames,
        incidence,
        work.density,
        report.initial_pressure,
        species,
        inner,
        dt,
        epsilon,
        number_scale,
    )
    length(report.edges) == length(plans) ||
        throw(ArgumentError("thermal overlap receipt omits or adds active edges"))
    current = [row.statistics for row in report.initial_frames]
    pressure = copy(report.initial_pressure)
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    energy_defect = energy_bound = zero(Q)
    momentum_defect = momentum_bound = (zero(Q), zero(Q), zero(Q))
    for k in eachindex(plans)
        plan, row = plans[k], report.edges[k]
        _electron_ion_record(row, _ELECTRON_ION_OVERLAP_EDGE_KEYS)
        row.species isa Int &&
            row.species == plan.species &&
            row.frame isa Int &&
            row.frame == plan.frame &&
            row.node isa Int &&
            row.node == plan.node &&
            _electron_ion_finite(row.weight) &&
            row.weight == plan.weight &&
            row.beta isa NTuple{3,Q} &&
            all(_electron_ion_finite, row.beta) &&
            row.beta == plan.beta &&
            row.cold isa Bool &&
            row.cold === plan.cold &&
            row.pressure_before isa T &&
            isequal(row.pressure_before, pressure[plan.node]) &&
            row.pressure_after isa T &&
            isfinite(row.pressure_after) &&
            row.pressure_after >= 0 &&
            _electron_ion_finite(row.target) &&
            row.target >= 0 &&
            _electron_ion_finite(row.electron_target) &&
            row.electron_target >= 0 &&
            _electron_ion_finite(row.energy_defect) &&
            _electron_ion_finite(row.energy_bound) &&
            row.energy_bound >= 0 &&
            row.entropy_bounds isa NTuple{2,Q} &&
            all(_electron_ion_finite, row.entropy_bounds) || throw(
            ArgumentError("thermal overlap edge order, coefficients or numeric leaves are invalid"),
        )
        f = (plan.species - 1) * ncells + plan.frame
        ps, bins = species[plan.species], work.bins.populations[plan.species]
        mass = _exact_rational(ps.m)
        _electron_ion_validate_statistics(row.before, mass, actual[f].statistics)
        _electron_ion_validate_statistics(row.after, mass, actual[f].statistics)
        row.before == current[f] || throw(ArgumentError("thermal overlap frame chain is broken"))
        Te = _exact_rational(row.pressure_before) / _exact_rational(work.density[plan.node])
        electron_target, targets = _electron_ion_targets(plan.Ce, Te, (row.before,), (plan.beta,))
        row.target == targets[1] && row.electron_target == electron_target ||
            throw(ArgumentError("thermal overlap edge differs from its conservative pair target"))
        _electron_ion_validate_materialized(
            row,
            ps,
            bins,
            plan.frame,
            T,
            _electron_ion_overlap_caps(row.after, mass, T),
        )
        !plan.cold &&
            (row.pressure_after <= 0 || row.after.thermal <= 0) &&
            throw(ArgumentError("thermal overlap receipt lost positive thermal energy"))
        change = row.after.kinetic - row.before.kinetic
        exact_pressure = _exact_rational(row.pressure_before) - gamma * change / volume
        expected = iszero(change) ? row.pressure_before : _from_exact(T, exact_pressure)
        exact_pressure >= 0 && isequal(expected, row.pressure_after) || throw(
            ArgumentError("thermal overlap receipt does not pay actual heat at its causal node"),
        )
        defect =
            change +
            volume * (_exact_rational(row.pressure_after) - _exact_rational(row.pressure_before)) /
            gamma
        bound = volume * _electron_ion_radius(row.pressure_after) / gamma
        row.energy_defect == defect && row.energy_bound == bound && abs(defect) <= bound ||
            throw(ArgumentError("thermal overlap edge energy identity is inconsistent"))
        entropy, entropy_bounds = _electron_ion_overlap_entropy(
            row.before,
            row.after,
            plan.Ce,
            work.density[plan.node],
            row.pressure_before,
            row.pressure_after,
            config.log_bits,
        )
        row.entropy === entropy && row.entropy_bounds == entropy_bounds ||
            throw(ArgumentError("thermal overlap actual edge entropy is inconsistent"))
        pressure[plan.node] = row.pressure_after
        current[f] = row.after
        energy_defect += defect
        energy_bound += bound
        momentum_defect = (
            momentum_defect[1] + row.after.momentum[1] - row.before.momentum[1],
            momentum_defect[2] + row.after.momentum[2] - row.before.momentum[2],
            momentum_defect[3] + row.after.momentum[3] - row.before.momentum[3],
        )
        momentum_bound = (
            momentum_bound[1] + row.materialized.momentum_bound[1],
            momentum_bound[2] + row.materialized.momentum_bound[2],
            momentum_bound[3] + row.materialized.momentum_bound[3],
        )
    end
    all(f -> current[f] == actual[f].statistics, eachindex(actual)) &&
        isequal(pressure, fields.pe) ||
        throw(ArgumentError("thermal overlap trace does not reach the saved physical endpoint"))
    kinetic_change = sum(
        (
            actual[f].statistics.kinetic - report.initial_frames[f].statistics.kinetic for
            f in eachindex(actual)
        );
        init = zero(Q),
    )
    electron_change =
        volume * sum(
            (
                _exact_rational(fields.pe[j]) - _exact_rational(report.initial_pressure[j]) for
                j in eachindex(pressure)
            );
            init = zero(Q),
        ) / gamma
    energy_defect == kinetic_change + electron_change &&
        report.energy_defect == energy_defect &&
        report.energy_bound == energy_bound &&
        abs(energy_defect) <= energy_bound &&
        report.momentum_defect == momentum_defect &&
        report.momentum_bound == momentum_bound &&
        all(c -> abs(momentum_defect[c]) <= momentum_bound[c], 1:3) || throw(
        ArgumentError("thermal overlap global native energy/momentum budget does not telescope"),
    )
    return nothing
end

# Validation uses current endpoint particles and pressure, not a collision/RNG
# replay. Fresh derived bins/deposition are private and never published.
function _validate_electron_ion_report(
    report,
    config,
    inner,
    species,
    fields,
    time,
    step,
    dt,
    epsilon,
    number_scale,
)
    config.policy == 2 && return _validate_electron_ion_overlap_report(
        report,
        config,
        inner,
        species,
        fields,
        time,
        step,
        dt,
        epsilon,
        number_scale,
    )
    if step == 0
        report === nothing ||
            throw(ArgumentError("unadvanced electron-ion state cannot carry a report"))
        return nothing
    end
    _electron_ion_record(report, _ELECTRON_ION_REPORT_KEYS)
    T = eltype(fields.pe)
    report.policy isa Int &&
        report.policy == 1 &&
        report.time isa T &&
        isequal(report.time, time) &&
        report.step isa Int &&
        report.step == step &&
        report.dt isa T &&
        isequal(report.dt, dt) &&
        report.cells isa Vector{NamedTuple} ||
        throw(ArgumentError("electron-ion receipt does not match its completed source clock"))
    if !_electron_ion_active(config)
        isempty(report.cells) ||
            throw(ArgumentError("inactive electron-ion state must have an empty receipt"))
        return nothing
    end
    length(report.cells) == length(fields.pe) ||
        throw(DimensionMismatch("electron-ion receipt cell count differs"))
    work = _electron_ion_work(config, inner, species)
    _electron_ion_density!(work, species, inner)
    for a in eachindex(species)
        ps = species[a]
        _bin_coulomb_particles!(
            work.bins.populations[a],
            ps,
            p -> _electron_ion_node(ps, inner.g, p),
        )
    end
    Q = _ElectronIonQ
    volume = prod(_exact_rational, inner.g.dx)
    gamma = _exact_rational(T(inner.model.closure.γ)) - 1
    pi_bounds = _electron_ion_pi_bounds(T)
    for cell in eachindex(fields.pe)
        receipt = report.cells[cell]
        _electron_ion_record(receipt, _ELECTRON_ION_CELL_KEYS)
        receipt.cell isa Int &&
            receipt.cell == cell &&
            receipt.density isa T &&
            isequal(receipt.density, work.density[cell]) &&
            receipt.pressure_before isa T &&
            isfinite(receipt.pressure_before) &&
            receipt.pressure_before >= 0 &&
            receipt.pressure_after isa T &&
            isequal(receipt.pressure_after, fields.pe[cell]) &&
            receipt.ions isa Vector{NamedTuple} &&
            length(receipt.ions) == length(species) &&
            _electron_ion_finite(receipt.electron_target) &&
            receipt.electron_target >= 0 &&
            all(
                _electron_ion_finite,
                (receipt.energy_defect, receipt.energy_bound, receipt.charge_defect),
            ) &&
            receipt.entropy_bounds isa NTuple{2,Q} &&
            all(_electron_ion_finite, receipt.entropy_bounds) ||
            throw(ArgumentError("electron-ion cell receipt has invalid endpoint values"))
        old = NamedTuple[]
        for a in eachindex(species)
            row, ps = receipt.ions[a], species[a]
            _electron_ion_record(row, _ELECTRON_ION_ION_KEYS)
            current = _electron_ion_stats(ps, work.bins.populations[a], cell)
            _electron_ion_validate_statistics(row.before, _exact_rational(ps.m), current)
            _electron_ion_validate_statistics(row.after, _exact_rational(ps.m), current)
            row.after == current &&
                row.beta isa NTuple{3,Q} &&
                all(_electron_ion_finite, row.beta) ||
                throw(ArgumentError("electron-ion receipt differs from stored endpoint particles"))
            _electron_ion_validate_materialized(row, ps, work.bins.populations[a], cell, T)
            push!(old, row.before)
        end
        Ne = _exact_rational(receipt.density) * volume
        Ce, Te =
            Ne / gamma, _exact_rational(receipt.pressure_before) / _exact_rational(receipt.density)
        active, beta, all_cold = _electron_ion_coefficients(
            config,
            old,
            species,
            Ce,
            Te,
            volume,
            dt,
            epsilon,
            number_scale,
            pi_bounds,
            T,
        )
        any(active) &&
            !all_cold &&
            receipt.pressure_after <= 0 &&
            throw(ArgumentError("electron-ion receipt lost positive electron thermal energy"))
        electron_target, targets = _electron_ion_targets(Ce, Te, old, beta)
        receipt.electron_target == electron_target && all(
            a -> receipt.ions[a].beta == beta[a] && receipt.ions[a].target == targets[a],
            eachindex(species),
        ) || throw(
            ArgumentError("electron-ion receipt differs from its physical coefficient/target map"),
        )
        change =
            sum((row.after.kinetic - row.before.kinetic for row in receipt.ions); init = zero(Q))
        exact_pressure = _exact_rational(receipt.pressure_before) - gamma * change / volume
        expected_pressure =
            iszero(change) ? receipt.pressure_before : _from_exact(T, exact_pressure)
        exact_pressure >= 0 && isequal(expected_pressure, receipt.pressure_after) ||
            throw(ArgumentError("electron-ion receipt does not pay actual ion kinetic change"))
        defect =
            change +
            volume *
            (_exact_rational(receipt.pressure_after) - _exact_rational(receipt.pressure_before)) /
            gamma
        bound = volume * _electron_ion_radius(receipt.pressure_after) / gamma
        receipt.energy_defect == defect &&
            receipt.energy_bound == bound &&
            abs(defect) <= bound &&
            receipt.charge_defect ==
            Ne - sum(a -> _exact_rational(species[a].q) * old[a].weight, eachindex(species)) ||
            throw(ArgumentError("electron-ion receipt energy or charge identity is inconsistent"))
        capacities, before, after = Q[Ce],
        Q[Te],
        Q[_exact_rational(receipt.pressure_after)/_exact_rational(receipt.density)]
        for a in eachindex(species)
            active[a] || continue
            !all_cold &&
                receipt.ions[a].after.thermal <= 0 &&
                throw(ArgumentError("electron-ion receipt lost positive ion thermal energy"))
            push!(capacities, old[a].capacity)
            push!(before, old[a].thermal / old[a].capacity)
            push!(after, receipt.ions[a].after.thermal / old[a].capacity)
        end
        unchanged =
            all(row -> row.before == row.after, receipt.ions) &&
            isequal(receipt.pressure_before, receipt.pressure_after)
        entropy, entropy_bounds =
            unchanged ? (:unchanged, (zero(Q), zero(Q))) :
            _electron_ion_entropy(capacities, before, after, config.log_bits)
        receipt.entropy === entropy && receipt.entropy_bounds == entropy_bounds ||
            throw(ArgumentError("electron-ion actual entropy receipt is inconsistent"))
    end
    return nothing
end
