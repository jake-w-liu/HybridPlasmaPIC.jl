_ace_manual_selection() = (
    kind = :manual,
    refine_threshold = nothing,
    coarsen_threshold = nothing,
    max_level = 0,
    buffer = 0,
    refine_tags = 0,
    coarsened = 0,
    balance_refined = 0,
)
function _ace_selection_check(selection, ::Type{T}) where {T}
    names = (
        :kind,
        :refine_threshold,
        :coarsen_threshold,
        :max_level,
        :buffer,
        :refine_tags,
        :coarsened,
        :balance_refined,
    )
    selection isa NamedTuple && keys(selection) == names ||
        throw(ArgumentError("adaptive composite selection schema differs"))
    all(
        i -> i isa Int && i >= 0,
        (
            selection.max_level,
            selection.buffer,
            selection.refine_tags,
            selection.coarsened,
            selection.balance_refined,
        ),
    ) || throw(ArgumentError("adaptive composite selection counts are invalid"))
    if selection.kind === :manual
        isequal(selection, _ace_manual_selection()) ||
            throw(ArgumentError("manual adaptive composite selection has tag metadata"))
    elseif selection.kind === :multilevel_curvature
        hi, lo = selection.refine_threshold, selection.coarsen_threshold
        hi isa T &&
            lo isa T &&
            isfinite(hi) &&
            isfinite(lo) &&
            0 <= lo <= hi &&
            (iszero(hi) || lo < hi) ||
            throw(ArgumentError("adaptive composite hysteresis thresholds are invalid"))
    else
        throw(ArgumentError("unsupported adaptive composite selection kind"))
    end
    return nothing
end
function _aca_denominators(nbase, maximum_level)
    denominators = [nbase]
    for level = 1:maximum_level
        push!(denominators, map(i -> Base.checked_mul(2, i), denominators[end]))
    end
    return denominators
end
function _aca_leaves(n::NTuple{D,Int}, refined) where {D}
    leaves = Tuple{Int,NTuple{D,Int}}[]
    function visit(a)
        if a in refined
            foreach(visit, _camr_children(a))
        else
            push!(leaves, a)
        end
        return nothing
    end
    for I in _camr_indices(n)
        visit((0, I))
    end
    return sort!(leaves)
end
function _aca_geometry(mesh, refined, denominators)
    addresses = _aca_leaves(mesh.nbase, refined)
    bounds = [_camr_bounds(a, mesh.L, denominators) for a in addresses]
    return (
        nbase = mesh.nbase,
        L = mesh.L,
        periodic = mesh.periodic,
        refined = Dict(a => i for (i, a) in enumerate(sort!(collect(refined)))),
        leaf_index = Dict(a => i for (i, a) in enumerate(addresses)),
        addresses = addresses,
        lower = first.(bounds),
        upper = last.(bounds),
    )
end
function _aca_ensure!(refined, address)
    ancestors = typeof(address)[]
    a = address
    while a[1] > 0
        a = _camr_parent(a)
        push!(ancestors, a)
    end
    for parent in reverse(ancestors)
        push!(refined, parent)
    end
    push!(refined, address)
    return nothing
end
function _aca_proposal(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    refine_threshold::T,
    coarsen_threshold::T,
    max_level::Integer,
    buffer::Integer,
) where {D,T}
    selection = (
        kind = :multilevel_curvature,
        refine_threshold = refine_threshold,
        coarsen_threshold = coarsen_threshold,
        max_level = 0,
        buffer = 0,
        refine_tags = 0,
        coarsened = 0,
        balance_refined = 0,
    )
    _ace_selection_check(selection, T)
    0 <= max_level <= typemax(Int) && 0 <= buffer <= typemax(Int) || throw(
        ArgumentError("adaptive composite max_level/buffer must be nonnegative Int-sized integers"),
    )
    maximum_level, margin = Int(max_level), Int(buffer)
    maximum_level >= maximum(first, st.mesh.addresses) ||
        throw(ArgumentError("adaptive composite max_level is below the existing deepest leaf"))
    denominators = _aca_denominators(st.mesh.nbase, maximum_level)
    oldmesh = st.mesh
    raw = _acq_raw_potential(oldmesh, st._state.phi)
    means = Dict{Tuple{Int,NTuple{D,Int}},NTuple{3,T}}()
    function mean(a)
        return get!(means, a) do
            l, u = _camr_bounds(a, oldmesh.L, denominators)
            value = map(v -> _from_exact(T, v), _acq_box_field(oldmesh, raw, l, u))
            all(isfinite, value) ||
                throw(ArgumentError("adaptive composite virtual-cell field average is nonfinite"))
            value
        end
    end
    function high(a, threshold)
        centre = mean(a)
        sizes = denominators[a[1]+1]
        for d = 1:D
            oldmesh.periodic || 1 < a[2][d] < sizes[d] || continue
            left = (
                a[1],
                ntuple(
                    j ->
                        j == d ? (oldmesh.periodic && a[2][j] == 1 ? sizes[j] : a[2][j] - 1) :
                        a[2][j],
                    D,
                ),
            )
            right = (
                a[1],
                ntuple(
                    j ->
                        j == d ? (oldmesh.periodic && a[2][j] == sizes[j] ? 1 : a[2][j] + 1) :
                        a[2][j],
                    D,
                ),
            )
            lm, rm = mean(left), mean(right)
            for c = 1:D
                refine_flags(T[lm[c], centre[c], rm[c]], threshold)[2] && return true
            end
        end
        return false
    end
    refined = Set{Tuple{Int,NTuple{D,Int}}}(st.mesh.refinements)
    coarsened = 0
    changed = true
    while changed
        changed = false
        for a in sort!(collect(refined); rev = true)
            children = _camr_children(a)
            any(child -> child in refined, children) && continue
            high(a, coarsen_threshold) && continue
            any(child -> high(child, coarsen_threshold), children) && continue
            delete!(refined, a)
            coarsened = Base.checked_add(coarsened, 1)
            changed = true
        end
    end
    tags = [
        a for a in _aca_leaves(oldmesh.nbase, refined) if
        a[1] < maximum_level && high(a, refine_threshold)
    ]
    for a in tags
        sizes = denominators[a[1]+1]
        ranges = ntuple(
            d ->
                oldmesh.periodic ? _aca_periodic_buffer(a[2][d], sizes[d], margin) :
                ((a[2][d]-min(margin, a[2][d] - 1)):(a[2][d]+min(margin, sizes[d] - a[2][d]))),
            D,
        )
        total = 1
        for range in ranges
            total = Base.checked_mul(total, length(range))
        end
        for I in CartesianIndices(ntuple(d -> length(ranges[d]), D))
            _aca_ensure!(refined, (a[1], ntuple(d -> ranges[d][I[d]], D)))
        end
    end
    balanced = 0
    while true
        geometry = _aca_geometry(oldmesh, refined, denominators)
        extra = Set{Tuple{Int,NTuple{D,Int}}}()
        for k in eachindex(geometry.addresses),
            j in _camr_neighbors(geometry, geometry.lower[k], geometry.upper[k])

            if geometry.addresses[k][1] + 1 < geometry.addresses[j][1]
                push!(extra, geometry.addresses[k])
            end
        end
        isempty(extra) && break
        for a in sort!(collect(extra))
            a[1] < maximum_level || throw(
                ArgumentError("adaptive composite balance closure exceeds requested maximum level"),
            )
            _aca_ensure!(refined, a)
            balanced = Base.checked_add(balanced, 1)
        end
    end
    mesh = CartesianAMRHierarchy(
        oldmesh.nbase,
        oldmesh.L,
        Tuple(sort!(collect(refined)));
        periodic = oldmesh.periodic,
    )
    selection = merge(
        selection,
        (
            max_level = maximum_level,
            buffer = margin,
            refine_tags = length(tags),
            coarsened = coarsened,
            balance_refined = balanced,
        ),
    )
    return mesh, selection
end
"Adapt the actual composite field with multilevel curvature, hysteresis and balance closure."
function adapt!(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    species::AbstractVector;
    refine_threshold::T,
    coarsen_threshold::T,
    max_level::Integer,
    buffer::Integer = 1,
) where {D,T}
    _ace_acquire(st)
    try
        _ace_ready(st, species)
        mesh, selection = _aca_proposal(st, refine_threshold, coarsen_threshold, max_level, buffer)
        return _ace_regrid!(st, species, mesh, selection)
    finally
        _ace_release(st)
    end
end

function _aca_periodic_buffer(index::Int, count::Int, margin::Int)
    margin >= count ÷ 2 && return 1:count
    # The short range has <=count entries. BigInt protects addition near Intmax;
    # all wrapped indices are converted only after reducing into1:count.
    return sort!(
        unique!([Int(mod(BigInt(index - 1) + offset, count)) + 1 for offset = -margin:margin]),
    )
end
