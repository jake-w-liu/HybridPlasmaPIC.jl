struct _AdaptiveCompositeSource
    total::Vector{_CompositeQ}
    inside::Vector{_CompositeQ}
    tail::Vector{_CompositeQ}
    background::Vector{_CompositeQ}
    strata::Vector{_CompositeQ}
end
struct _AdaptiveCompositeCrossing{D}
    population::Int
    particle::Int
    parameter::_CompositeQ
    from::Int
    to::Int
    axes::NTuple{D,Bool}
end
_acc_strata(::Val{D}) where {D} = [ntuple(d -> (k ÷ (3^(d-1))) % 3 - 1, D) for k = 0:(3^D-1)]
function _acc_moments(x, lo, hi, a, L)
    left = _yee_cloud_primitives(x - lo, a)[2]
    right = _yee_cloud_primitives(x - hi, a)[2]
    A = _yee_cloud_edge(x, lo, hi - lo, a)[1]
    interior = (left - A, A - right)
    full = (interior[1] + (iszero(lo) ? 1 - left : zero(_CompositeQ)), interior[2] + (hi == L ? right : zero(_CompositeQ)))
    all(>=(zero(_CompositeQ)), (interior..., full...)) || throw(ArgumentError("composite cloud moment lost positive charge partition"))
    return (interior = interior, full = full, mass = left - right)
end
function _acc_span(mesh::CartesianAMRHierarchy{D}, x, widths) where {D}
    mesh.periodic && return _pacc_span(mesh,x,widths)
    lo = ntuple(d -> max(zero(_CompositeQ), x[d] - widths[d]), D)
    hi = ntuple(d -> min(_exact_rational(mesh.L[d]), x[d] + widths[d]), D)
    return _camr_query(mesh, lo, hi)
end
function _acc_particle_vectors(mesh::CartesianAMRHierarchy{D}, species, clouds) where {D}
    count = length(mesh.masters)
    particle = fill(zero(_CompositeQ), count)
    inside = mesh.periodic ? nothing : fill(zero(_CompositeQ), count)
    strata = mesh.periodic ? nothing : fill(zero(_CompositeQ), 3^D)
    bits, strata_axes = _camr_bits(Val(D)), _acc_strata(Val(D))
    widths = map(c -> _exact_rational(c.half_width), clouds)
    for ps in species, p in eachindex(ps.weight)
        charge = _yee_macrocharge(ps, p, Val(true))
        if !mesh.periodic
            iszero(charge) && continue
        end
        x = ntuple(d -> _exact_rational(ps.x[d][p]), D)
        factors = if mesh.periodic
            nothing
        else
            ntuple(D) do d
                a, L = widths[d], _exact_rational(mesh.L[d])
                left, right = _yee_cloud_primitives(x[d], a)[2], _yee_cloud_primitives(x[d] - L, a)[2]
                (1 - left, left - right, right)
            end
        end
        (mesh.periodic || all(f -> all(>=(zero(_CompositeQ)), f) && sum(f) == 1, factors)) || throw(ArgumentError("composite cloud axis strata lost unit mass"))
        full_mass, inside_mass = zero(_CompositeQ), zero(_CompositeQ)
        if !mesh.periodic
            for (i, sigma) in enumerate(strata_axes)
                strata[i] += charge * prod(d -> factors[d][sigma[d]+2], 1:D)
            end
        end
        for k in _acc_span(mesh, x, widths)
            moments = ntuple(d -> _acc_axis_moments(mesh, x[d], mesh.lower[k][d], mesh.upper[k][d], widths[d], _exact_rational(mesh.L[d])), D)
            for (a, raw) in enumerate(mesh.corners[k])
                unit_full = prod(d -> moments[d].full[bits[a][d]+1], 1:D)
                unit_inside = mesh.periodic ? zero(_CompositeQ) : prod(d -> moments[d].interior[bits[a][d]+1], 1:D)
                full_mass += unit_full
                inside_mass += unit_inside
                full, volume = charge * unit_full, charge * unit_inside
                for (master, coefficient) in mesh.rows[raw]
                    particle[master] += coefficient * full
                    mesh.periodic || (inside[master] += coefficient * volume)
                end
            end
        end
        full_mass == 1 && (mesh.periodic || inside_mass == prod(d -> factors[d][2], 1:D)) || throw(ArgumentError("composite individual cloud lost physical mass"))
    end
    total_charge = sum((_yee_macrocharge(ps, p, Val(true)) for ps in species for p in eachindex(ps.weight)); init = zero(_CompositeQ))
    sum(particle; init = zero(_CompositeQ)) == total_charge || throw(ArgumentError("composite physical cloud lost total charge"))
    if !mesh.periodic
        sum(strata; init = zero(_CompositeQ)) == total_charge || throw(ArgumentError("composite physical strata lost total charge"))
        mid = findfirst(sigma -> all(iszero, sigma), strata_axes)
        sum(inside; init = zero(_CompositeQ)) == strata[mid] || throw(ArgumentError("composite inside cloud charge differs from physical integral"))
    end
    return (particle = particle,inside = inside,strata = strata)
end
_acc_axis_moments(mesh,x,lo,hi,a,L) = mesh.periodic ? _pacc_axis_moments(x,lo,hi,a,L) : _acc_moments(x,lo,hi,a,L)
function _acc_background(mesh,background)
    bkg = fill(zero(_CompositeQ),length(mesh.masters))
    rho = _exact_rational(background)
    D = length(mesh.L)
    for k in eachindex(mesh.addresses), raw in mesh.corners[k]
        value = rho * mesh.volumes[k] / (1 << D)
        for (master, coefficient) in mesh.rows[raw]
            bkg[master] += value * coefficient
        end
    end
    sum(bkg; init = zero(_CompositeQ)) == rho * prod(_exact_rational, mesh.L) || throw(ArgumentError("composite background volume normalization differs"))
    return bkg
end
function _acc_source(mesh::CartesianAMRHierarchy{D},species,clouds,background;periodic_neutrality = :require) where {D}
    charge = _acc_particle_vectors(mesh,species,clouds)
    mesh.periodic && return _pacc_source(mesh,species,charge.particle,background,periodic_neutrality)
    bkg = _acc_background(mesh,background)
    tail = charge.particle-charge.inside
    all(i -> iszero(tail[i]),mesh.free) || throw(ArgumentError("composite wall tail reached an interior master"))
    return _AdaptiveCompositeSource(charge.particle+bkg,charge.inside,tail,bkg,charge.strata)
end
function _acc_gather(mesh::CartesianAMRHierarchy{D}, raw, x, clouds) where {D}
    result = fill(zero(_CompositeQ), 3)
    full_mass = zero(_CompositeQ)
    bits = _camr_bits(Val(D))
    widths = map(c -> _exact_rational(c.half_width), clouds)
    for k in _acc_span(mesh, x, widths)
        moments = ntuple(d -> _acc_axis_moments(mesh, x[d], mesh.lower[k][d], mesh.upper[k][d], widths[d], _exact_rational(mesh.L[d])), D)
        full_mass += sum(bit -> prod(d -> moments[d].full[bit[d]+1], 1:D), bits; init = zero(_CompositeQ))
        for (a, vertex) in enumerate(mesh.corners[k]), d = 1:D
            derivative = (iszero(bits[a][d]) ? -one(_CompositeQ) : one(_CompositeQ)) * moments[d].mass / (mesh.upper[k][d] - mesh.lower[k][d])
            for j = 1:D
                j == d || (derivative *= moments[j].full[bits[a][j]+1])
            end
            result[d] -= raw[vertex] * derivative
        end
    end
    full_mass == 1 || throw(ArgumentError("composite force cloud support lost unit mass"))
    return Tuple(result)
end
function _acc_bins(species, mesh)
    count = length(mesh.addresses)
    counts = zeros(Int, count)
    locations = Vector{Vector{Int}}()
    total = 0
    for ps in species
        total = Base.checked_add(total, length(ps.weight))
        row = Int[]
        for p in eachindex(ps.weight)
            k = _camr_locate(mesh, map(d -> _exact_rational(ps.x[d][p]), eachindex(ps.x)))
            push!(row, k)
            counts[k] = Base.checked_add(counts[k], 1)
        end
        push!(locations, row)
    end
    offsets = Vector{Int}(undef, Base.checked_add(count, 1))
    offsets[1] = 1
    for k = 1:count
        offsets[k+1] = Base.checked_add(offsets[k], counts[k])
    end
    slots = copy(offsets[1:count])
    populations, particles = Vector{Int}(undef, total), Vector{Int}(undef, total)
    for (s, row) in enumerate(locations), (p, k) in enumerate(row)
        populations[slots[k]], particles[slots[k]] = s, p
        slots[k] += 1
    end
    return _CompositeParticleBins(offsets, populations, particles)
end
function _acc_validate_bins(bins, species, mesh)
    expected = _acc_bins(species, mesh)
    all(name -> isequal(getfield(bins, name), getfield(expected, name)), fieldnames(_CompositeParticleBins)) || throw(ArgumentError("composite particle ownership bins changed"))
    return nothing
end
function _acc_crossings(old, new, mesh::CartesianAMRHierarchy{D}) where {D}
    by_leaf = [Int[] for _ in mesh.addresses]
    for (i, face) in enumerate(mesh.faces)
        iszero(face.left) || push!(by_leaf[face.left], i)
        iszero(face.right) || push!(by_leaf[face.right], i)
    end
    events = _AdaptiveCompositeCrossing{D}[]
    for (s, (before, after)) in enumerate(zip(old, new)), p in eachindex(before.weight)
        x, y = ntuple(d -> _exact_rational(before.x[d][p]), D), ntuple(d -> _exact_rational(after.x[d][p]), D)
        x == y && continue
        hits = _acc_segment_hits(mesh,x,y,by_leaf)
        parameters = sort!(unique!(vcat(_CompositeQ[0,1], collect(keys(hits)))))
        owner = _camr_locate(mesh, x)
        for i in eachindex(parameters)
            parameter = parameters[i]
            destination = if i == length(parameters)
                _camr_locate(mesh, y)
            else
                middle = (parameter + parameters[i+1]) / 2
                _camr_locate(mesh, ntuple(d -> x[d] + middle * (y[d] - x[d]), D))
            end
            if destination != owner
                axes = get(hits, parameter, ntuple(_ -> false, D))
                any(axes) || throw(ArgumentError("composite crossing changed owner without an interface"))
                push!(events, _AdaptiveCompositeCrossing(s, p, parameter, owner, destination, axes))
            end
            owner = destination
        end
    end
    return events
end

function _acc_segment_hits(mesh::CartesianAMRHierarchy{D},x,y,by_leaf) where {D}
    faces = Set{Int}()
    for k in _camr_segment(mesh, x, y)
        union!(faces, by_leaf[k])
    end
    hits = Dict{_CompositeQ,NTuple{D,Bool}}()
    for i in sort!(collect(faces))
        face = mesh.faces[i]
        face.left > 0 && face.right > 0 || continue
        d, delta = face.axis, y[face.axis] - x[face.axis]
        iszero(delta) && continue
        parameter = (face.lower[d] - x[d]) / delta
        0 <= parameter <= 1 || continue
        point = ntuple(j -> x[j] + parameter * (y[j] - x[j]), D)
        all(j -> face.lower[j] <= point[j] <= face.upper[j], 1:D) || continue
        previous = get(hits, parameter, ntuple(_ -> false, D))
        hits[parameter] = ntuple(j -> previous[j] || j == d, D)
    end
    return hits
end
