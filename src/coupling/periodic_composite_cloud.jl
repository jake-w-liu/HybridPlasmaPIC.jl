# Periodic Q1 moments integrate the same original free-space cloud. Image bands
# are exact integer polynomials; no count-proportional image loop or fit exists.
struct _PeriodicCompositeSource
    total::Vector{_CompositeQ}
    particle::Vector{_CompositeQ}
    background::Vector{_CompositeQ}
    background_density::_CompositeQ
    declared_background::_CompositeQ
    neutralizing_charge::_CompositeQ
end
struct _PeriodicCompositeTransport{D,T<:Union{Float32,Float64}}
    old_positions::Vector{NTuple{D,Vector{T}}}
    half_velocity::Vector{NTuple{D,Vector{T}}}
    winding::Vector{NTuple{D,Vector{BigInt}}}
    defect::Vector{NTuple{D,Vector{_CompositeQ}}}
    tiles_before::Vector{NTuple{D,BigInt}}
    tiles_after::Vector{NTuple{D,BigInt}}
end
function _pacc_choose(N::BigInt, k::Int)
    result = one(BigInt)
    for j = 1:k
        result = div(result * (N - j + 1), j)
    end
    return result
end
function _pacc_axis_moments(x, lo, hi, a, L)
    first = floor(BigInt, (lo - x - a) / L) + 1
    last = ceil(BigInt, (hi - x + a) / L) - 1
    result = fill(zero(_CompositeQ), 3)
    if first <= last
        cuts = BigInt[first, last+1]
        for endpoint in (lo, hi), knot in (-a, zero(a), a)
            cut = floor(BigInt, (knot - x + endpoint) / L) + 1
            first < cut <= last && push!(cuts, cut)
        end
        sort!(unique!(cuts))
        for band = 1:length(cuts)-1
            origin, N = cuts[band], cuts[band+1] - cuts[band]
            degree = Int(min(BigInt(3), N - 1))
            values = NTuple{3,_CompositeQ}[]
            for j = 0:degree
                moment = _acc_moments(x + (origin + j) * L, lo, hi, a, L)
                push!(values, (moment.interior..., moment.mass))
            end
            for order = 0:degree
                coefficient = _pacc_choose(N, order + 1)
                for c = 1:3
                    result[c] += coefficient * values[1][c]
                end
                if order < degree
                    values =
                        [ntuple(c -> values[j+1][c] - values[j][c], 3) for j = 1:length(values)-1]
                end
            end
        end
    end
    all(>=(zero(_CompositeQ)), result) && result[1] + result[2] == result[3] && result[3] <= 1 ||
        throw(ArgumentError("periodic composite axis moments lost positive mass"))
    pair = (result[1], result[2])
    return (interior = pair, full = pair, mass = result[3])
end
function _pacc_span(mesh::CartesianAMRHierarchy{D}, x, widths) where {D}
    intervals = ntuple(D) do d
        L, a = _exact_rational(mesh.L[d]), widths[d]
        if 2a >= L
            [(zero(_CompositeQ), L)]
        elseif x[d] - a < 0
            [(zero(_CompositeQ), x[d] + a), (x[d] - a + L, L)]
        elseif x[d] + a > L
            [(x[d] - a, L), (zero(_CompositeQ), x[d] + a - L)]
        else
            [(x[d] - a, x[d] + a)]
        end
    end
    leaves = Int[]
    for box in Iterators.product(intervals...)
        append!(leaves, _camr_query(mesh, ntuple(d -> box[d][1], D), ntuple(d -> box[d][2], D)))
    end
    return sort!(unique!(leaves))
end
function _pacc_source(mesh, species, particle, background, neutrality)
    neutrality in (:require, :neutralize) ||
        throw(ArgumentError("periodic neutrality must be :require or :neutralize"))
    charges = _pic_collection_charges(mesh, background, species)
    neutrality === :require &&
        !iszero(charges.net) &&
        throw(ArgumentError("periodic composite Poisson requires exact physical neutrality"))
    neutralizing = neutrality === :neutralize ? -charges.net : zero(_CompositeQ)
    density = _exact_rational(background) + neutralizing / prod(_exact_rational, mesh.L)
    bkg = _acc_background(mesh, density)
    total = particle + bkg
    iszero(sum(total; init = zero(_CompositeQ))) ||
        throw(ArgumentError("periodic composite source has a nonzero original null mode"))
    return _PeriodicCompositeSource(total, particle, bkg, density, charges.background, neutralizing)
end
function _pacc_transport(species, ::Val{D}, ::Type{T}) where {D,T}
    old = NTuple{D,Vector{T}}[map(copy, ps.x) for ps in species]
    half = NTuple{D,Vector{T}}[ntuple(d -> copy(ps.v[d]), D) for ps in species]
    winding = NTuple{D,Vector{BigInt}}[
        ntuple(_ -> [zero(BigInt) for _ in eachindex(ps.weight)], D) for ps in species
    ]
    defect = NTuple{D,Vector{_CompositeQ}}[
        ntuple(_ -> fill(zero(_CompositeQ), length(ps.weight)), D) for ps in species
    ]
    return _PeriodicCompositeTransport{D,T}(
        old,
        half,
        winding,
        defect,
        NTuple{D,BigInt}[],
        NTuple{D,BigInt}[],
    )
end
function _pacc_transport_arrays(transport::_PeriodicCompositeTransport)
    arrays = Any[
        transport.old_positions,
        transport.half_velocity,
        transport.winding,
        transport.defect,
        transport.tiles_before,
        transport.tiles_after,
    ]
    for roles in
        (transport.old_positions, transport.half_velocity, transport.winding, transport.defect),
        row in roles

        append!(arrays, row)
    end
    return arrays
end
function _pacc_crossings(mesh::CartesianAMRHierarchy{D}, old_positions, new, winding) where {D}
    by_leaf = [Int[] for _ in mesh.addresses]
    for (i, face) in enumerate(mesh.faces)
        face.left > 0 && face.right > 0 || continue
        push!(by_leaf[face.left], i)
        push!(by_leaf[face.right], i)
    end
    events = _AdaptiveCompositeCrossing{D}[]
    before_tiles, after_tiles = NTuple{D,BigInt}[], NTuple{D,BigInt}[]
    periods = map(_exact_rational, mesh.L)
    for (s, ps) in enumerate(new), p in eachindex(ps.weight)
        x = ntuple(d -> _exact_rational(old_positions[s][d][p]), D)
        y = ntuple(d -> _exact_rational(ps.x[d][p]) + winding[s][d][p] * periods[d], D)
        x == y && continue
        hits = Dict{_CompositeQ,NTuple{D,Bool}}()
        seams = _CompositeQ[0, 1]
        ranges = ntuple(D) do d
            if iszero(y[d] - x[d])
                (one(BigInt), zero(BigInt))
            else
                (
                    ceil(BigInt, min(x[d], y[d]) / periods[d]),
                    floor(BigInt, max(x[d], y[d]) / periods[d]),
                )
            end
        end
        plane_count = sum(range -> max(zero(BigInt), range[2] - range[1] + 1), ranges)
        plane_count + 2 <= typemax(Int) ||
            throw(ArgumentError("periodic transport seam count does not fit Int"))
        for d = 1:D
            delta = y[d] - x[d]
            iszero(delta) && continue
            for plane = ranges[d][1]:ranges[d][2]
                parameter = (plane * periods[d] - x[d]) / delta
                push!(seams, parameter)
                _pacc_hit!(hits, parameter, ntuple(j -> j == d, D))
            end
        end
        sort!(unique!(seams))
        for i = 1:length(seams)-1
            a, b = seams[i], seams[i+1]
            middle = (a + b) / 2
            tile = ntuple(d -> floor(BigInt, (x[d] + middle * (y[d] - x[d])) / periods[d]), D)
            left = ntuple(d -> x[d] + a * (y[d] - x[d]) - tile[d] * periods[d], D)
            right = ntuple(d -> x[d] + b * (y[d] - x[d]) - tile[d] * periods[d], D)
            for (local_parameter, axes) in _acc_segment_hits(mesh, left, right, by_leaf)
                parameter = a + (b - a) * local_parameter
                _pacc_hit!(hits, parameter, axes)
            end
        end
        parameters = sort!(unique!(vcat(_CompositeQ[0, 1], collect(keys(hits)))))
        owner, owner_tile = _camr_locate(mesh, x), ntuple(_ -> zero(BigInt), D)
        for i in eachindex(parameters)
            parameter = parameters[i]
            if i == length(parameters)
                tile = ntuple(d -> winding[s][d][p], D)
                point = ntuple(d -> _exact_rational(ps.x[d][p]), D)
            else
                middle = (parameter + parameters[i+1]) / 2
                lifted = ntuple(d -> x[d] + middle * (y[d] - x[d]), D)
                tile = ntuple(d -> floor(BigInt, lifted[d] / periods[d]), D)
                point = ntuple(d -> lifted[d] - tile[d] * periods[d], D)
            end
            destination = _camr_locate(mesh, point)
            if destination != owner || tile != owner_tile
                axes = get(hits, parameter, ntuple(_ -> false, D))
                any(axes) && destination != owner || throw(
                    ArgumentError("periodic transport changed tile/owner without a physical cut"),
                )
                Base.checked_add(length(events), 1)
                push!(events, _AdaptiveCompositeCrossing(s, p, parameter, owner, destination, axes))
                push!(before_tiles, owner_tile)
                push!(after_tiles, tile)
            end
            owner, owner_tile = destination, tile
        end
    end
    return events, before_tiles, after_tiles
end
function _pacc_transport_check(
    mesh::CartesianAMRHierarchy{D,T},
    species,
    transport,
    report,
) where {D,T}
    transport isa _PeriodicCompositeTransport{D,T} ||
        throw(ArgumentError("periodic transport receipt type differs"))
    all(
        rows -> length(rows) == length(species),
        (transport.old_positions, transport.half_velocity, transport.winding, transport.defect),
    ) || throw(ArgumentError("periodic transport population layout differs"))
    _ace_independent(_pacc_transport_arrays(transport))
    for (s, ps) in enumerate(species)
        n = length(ps.weight)
        for rows in
            (transport.old_positions, transport.half_velocity, transport.winding, transport.defect),
            array in rows[s]

            length(array) == n && axes(array, 1) == Base.OneTo(n) ||
                throw(DimensionMismatch("periodic transport component counts differ"))
        end
        for p in eachindex(ps.weight), d = 1:D
            old, half = transport.old_positions[s][d][p], transport.half_velocity[s][d][p]
            isfinite(old) &&
                isfinite(half) &&
                0 <= _exact_rational(old) < _exact_rational(mesh.L[d]) &&
                _sce_q_finite(transport.defect[s][d][p]) ||
                throw(ArgumentError("periodic transport operands/defect are invalid"))
            ideal = _particle_coordinate(
                _exact_rational(old),
                _exact_rational(half),
                _exact_rational(report.dt),
            )
            represented = _periodic_represented_lift(ideal, mesh.L[d])
            isequal(represented.wrapped, ps.x[d][p]) &&
                represented.winding == transport.winding[s][d][p] &&
                represented.defect == transport.defect[s][d][p] || throw(
                ArgumentError(
                    "periodic transport lift/rounding receipt differs from actual endpoint",
                ),
            )
        end
    end
    events, before, after =
        _pacc_crossings(mesh, transport.old_positions, species, transport.winding)
    event_values(xs) = [(e.population, e.particle, e.parameter, e.from, e.to, e.axes) for e in xs]
    isequal(event_values(events), event_values(report.crossings)) &&
        isequal(before, transport.tiles_before) &&
        isequal(after, transport.tiles_after) ||
        throw(ArgumentError("periodic tile/leaf crossing receipt differs from represented path"))
    return nothing
end

function _pacc_hit!(hits::Dict{_CompositeQ,NTuple{D,Bool}}, parameter, axes) where {D}
    # Reserve both interval endpoints in the eventual parameter vector as well
    # as a newly distinct hit before mutating private history storage.
    haskey(hits, parameter) || Base.checked_add(length(hits), 3)
    previous = get(hits, parameter, ntuple(_ -> false, D))
    hits[parameter] = ntuple(d -> previous[d] || axes[d], D)
    return nothing
end
