# Physical-pencil particle ownership and CIC data transport. No global field copy.
const _PencilParticle{T} = HybridPlasmaPIC._ElectrostaticHostParticle{3,T}
const _PencilNode = NTuple{3,Int}
const _PencilCharge = Rational{BigInt}

function _pencil_pic_population(ps::HybridPlasmaPIC.ParticleSet)
    HybridPlasmaPIC._check_particle_storage_axes(ps)
    HybridPlasmaPIC._require_particle_ordinary_stores(ps, "distributed PIC")
    HybridPlasmaPIC._particle_components_mightalias(ps) &&
        throw(ArgumentError("distributed PIC particle components must be independent"))
    HybridPlasmaPIC._validated_boris_factor(ps)
    HybridPlasmaPIC._validate_particle_positions(ps)
    HybridPlasmaPIC._require_finite_particle_velocities(ps, "distributed PIC")
    all(w -> isfinite(w) && HybridPlasmaPIC._stored_nonnegative(w), ps.weight) ||
        throw(ArgumentError("distributed PIC weights must be finite and nonnegative"))
    return nothing
end

function _pencil_pic_phase_values(
    initialized::Bool,
    primed::Bool,
    last_dt::T,
    time::T,
    step::Int,
) where {T}
    isfinite(time) && time >= zero(time) && step >= 0 ||
        throw(ArgumentError("distributed PIC clock is invalid"))
    isfinite(last_dt) && (
        primed ? HybridPlasmaPIC._stored_positive(last_dt) :
        HybridPlasmaPIC._particle_stored_zero(last_dt)
    ) || throw(ArgumentError("distributed PIC last_dt must be positive exactly when primed"))
    primed &&
        !initialized &&
        throw(ArgumentError("a primed distributed PIC model must be initialized"))
    if step == 0
        HybridPlasmaPIC._particle_stored_zero(time) && !primed || throw(
            ArgumentError(
                "an unstepped distributed PIC model must have zero time and unprimed velocity",
            ),
        )
    else
        time > zero(time) && primed && last_dt <= time || throw(
            ArgumentError(
                "a stepped distributed PIC model must have positive time and primed velocity",
            ),
        )
    end
    HybridPlasmaPIC._require_resolved_dt("distributed PIC history", last_dt)
    return nothing
end

function _pencil_pic_charges(model, species, context)
    comm = model.g.fft.comm
    local_charges =
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "physical charge preparation") do
            HybridPlasmaPIC._pic_collection_charges(model.g, zero(model.background_charge), species)
        end
    particle = _pencil_exact_scalar_sum(
        comm,
        context,
        "particle charge serialization",
        "particle charge accumulation",
    ) do
        local_charges.particle
    end
    absolute = _pencil_exact_scalar_sum(
        comm,
        context,
        "absolute charge serialization",
        "absolute charge accumulation",
    ) do
        local_charges.scale
    end
    return HybridPlasmaPIC._mpi_collective_stage(comm, context, "physical neutrality") do
        background =
            HybridPlasmaPIC._exact_rational(model.background_charge) *
            prod(HybridPlasmaPIC._exact_rational, model.g.L)
        charges = (
            particle = particle,
            background = background,
            net = particle + background,
            scale = absolute + abs(background),
        )
        HybridPlasmaPIC._espic_check_neutrality(model, charges)
        charges
    end
end

struct _PencilParticleLayout
    by_rank::Vector{NTuple{3,UnitRange{Int}}}
    rank::Int
    partitions::NTuple{3,Vector{UnitRange{Int}}}
    starts::NTuple{3,Vector{Int}}
    lookup::Array{Int,3}
end

function _pencil_partition_cover(partitions, starts, n)
    for d = 1:3
        length(starts[d]) == length(partitions[d]) ||
            throw(ArgumentError("particle partition starts have invalid length"))
        previous = 0
        for (i, interval) in enumerate(partitions[d])
            previous < n[d] &&
                first(interval) == previous + 1 &&
                first(interval) <= last(interval) <= n[d] &&
                starts[d][i] == first(interval) ||
                throw(ArgumentError("particle axis partitions do not form an exact ordered cover"))
            previous = last(interval)
        end
        previous == n[d] || throw(ArgumentError("particle axis partitions omit cells"))
    end
    return nothing
end

function _pencil_partition_coordinates(partitions, starts, owned)
    return ntuple(3) do d
        i = searchsortedfirst(starts[d], first(owned[d]))
        i <= length(partitions[d]) && partitions[d][i] == owned[d] ||
            throw(ArgumentError("physical pencil ranges do not form a Cartesian partition product"))
        i
    end
end

function _pencil_particle_layout(g)
    pencil = PencilArrays.pencil(g.physical)
    comm = PencilArrays.get_comm(pencil)
    MPI.Comm_compare(comm, g.fft.comm) in (MPI.IDENT, MPI.CONGRUENT) ||
        throw(ArgumentError("particle and Fourier communicators must have the same rank ordering"))
    topology = PencilArrays.topology(pencil)
    nranks = MPI.Comm_size(g.fft.comm)
    ranges = Vector{NTuple{3,UnitRange{Int}}}(undef, nranks)
    seen = falses(nranks)
    for I in CartesianIndices(size(topology))
        rank = Int(MPI.Cart_rank(comm, [i - 1 for i in Tuple(I)])) + 1
        1 <= rank <= nranks && !seen[rank] ||
            throw(ArgumentError("physical pencil topology does not identify each rank once"))
        remote = PencilArrays.range_remote(pencil, I, PencilArrays.LogicalOrder())
        ranges[rank] = ntuple(d -> Int(first(remote[d])):Int(last(remote[d])), 3)
        seen[rank] = true
    end
    all(seen) || throw(ArgumentError("physical pencil topology omits a rank"))
    partitions = ntuple(d -> sort!(unique([owned[d] for owned in ranges]); by = first), 3)
    starts = map(axis -> [first(interval) for interval in axis], partitions)
    _pencil_partition_cover(partitions, starts, g.n)
    prod(axis -> BigInt(length(axis)), partitions) == nranks ||
        throw(ArgumentError("physical pencil partition product does not match the communicator"))
    lookup = zeros(Int, map(length, partitions))
    for (rank, owned) in enumerate(ranges)
        coordinates = _pencil_partition_coordinates(partitions, starts, owned)
        iszero(lookup[coordinates...]) ||
            throw(ArgumentError("physical pencil ranks duplicate a partition"))
        lookup[coordinates...] = rank
    end
    all(>(0), lookup) || throw(ArgumentError("physical pencil partition product omits a rank"))
    rank = MPI.Comm_rank(g.fft.comm) + 1
    ranges[rank] == PencilArrays.range_local(pencil, PencilArrays.LogicalOrder()) ||
        throw(ArgumentError("physical pencil local ranges disagree with its rank table"))
    return _PencilParticleLayout(ranges, rank, partitions, starts, lookup)
end

# Recheck only immutable layout representation, not physical values or clocks.
function _pencil_particle_layout_storage(layout::_PencilParticleLayout, g)
    pencil = PencilArrays.pencil(g.physical)
    comm = PencilArrays.get_comm(pencil)
    MPI.Comm_compare(comm, g.fft.comm) in (MPI.IDENT, MPI.CONGRUENT) ||
        throw(ArgumentError("particle and Fourier communicator rank ordering changed"))
    length(layout.by_rank) == MPI.Comm_size(g.fft.comm) &&
        layout.rank == MPI.Comm_rank(g.fft.comm) + 1 ||
        throw(ArgumentError("particle ownership table does not match the communicator"))
    _pencil_partition_cover(layout.partitions, layout.starts, g.n)
    size(layout.lookup) == map(length, layout.partitions) &&
        prod(axis -> BigInt(length(axis)), layout.partitions) == length(layout.by_rank) ||
        throw(ArgumentError("particle lookup dimensions do not match the physical pencil"))
    topology = PencilArrays.topology(pencil)
    coordinates = Vector{Cint}(undef, length(size(topology)))
    for I in CartesianIndices(size(topology))
        for d in eachindex(coordinates)
            coordinates[d] = I[d] - 1
        end
        rank = Int(MPI.Cart_rank(comm, coordinates)) + 1
        layout.by_rank[rank] == PencilArrays.range_remote(pencil, I, PencilArrays.LogicalOrder()) ||
            throw(ArgumentError("particle ownership table differs from the physical pencil"))
        cell = _pencil_partition_coordinates(layout.partitions, layout.starts, layout.by_rank[rank])
        layout.lookup[cell...] == rank ||
            throw(ArgumentError("particle owner lookup differs from the physical pencil"))
    end
    return nothing
end

@inline _pencil_contains(ranges, idx) = all(d -> idx[d] in ranges[d], 1:3)

function _pencil_node_owner(layout::_PencilParticleLayout, idx::_PencilNode)
    _pencil_contains(layout.by_rank[layout.rank], idx) && return layout.rank
    all(d -> 1 <= idx[d] <= last(last(layout.partitions[d])), 1:3) ||
        throw(ArgumentError("CIC node is outside the physical pencil ownership table"))
    coordinates = ntuple(d -> searchsortedlast(layout.starts[d], idx[d]), 3)
    return layout.lookup[coordinates...]
end

Base.@propagate_inbounds @inline function HybridPlasmaPIC._particle_cell_position(
    ps::HybridPlasmaPIC.ParticleSet{3,T},
    g::DistributedFourierGrid{T},
    d::Int,
    p::Int,
) where {T}
    x = ps.x[d][p]
    isfinite(x) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
    return HybridPlasmaPIC._periodic_cell_position(x, g.L[d], g.dx[d])
end

@inline _pencil_cic_stencil(ps, g, p) = ntuple(
    d -> HybridPlasmaPIC._stencil1d(
        HybridPlasmaPIC.CIC(),
        HybridPlasmaPIC._particle_cell_position(ps, g, d, p),
    ),
    3,
)
@inline _pencil_cic_node(st, offset, n) =
    ntuple(d -> HybridPlasmaPIC._wrapped_cell(st[d][1] + offset[d] - 1, n[d]), 3)
@inline _pencil_particle_owner(ps, g, layout, p) =
    _pencil_node_owner(layout, _pencil_cic_node(_pencil_cic_stencil(ps, g, p), (1, 1, 1), g.n))

function _pencil_owned_particles(species, g, layout)
    for ps in species, p in eachindex(ps.weight)
        all(d -> zero(eltype(ps)) <= ps.x[d][p] < g.L[d], 1:3) || throw(
            ArgumentError(
                "particle positions must be in the canonical periodic interval; refresh the distributed field",
            ),
        )
        _pencil_particle_owner(ps, g, layout, p) == layout.rank || throw(
            ArgumentError("particle is not on its CIC anchor owner; refresh the distributed field"),
        )
    end
    return nothing
end

# Every kept subset is detached, even when no particle moves. The retained
# pre-push cache is never published and can safely be reused on the next call.
function _pencil_migrate(species::Vector{P}, g, layout, context) where {P}
    nranks = length(layout.by_rank)
    prepared, outgoing =
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "migration preparation") do
            result = P[]
            packets = [Tuple{Int,P}[] for _ = 1:nranks]
            for (s, ps) in enumerate(species)
                indices = [Int[] for _ = 1:nranks]
                for p in eachindex(ps.weight)
                    push!(indices[_pencil_particle_owner(ps, g, layout, p)], p)
                end
                push!(
                    result,
                    HybridPlasmaPIC._prepare_particle_replacement(ps, indices[layout.rank]),
                )
                for rank = 1:nranks
                    rank == layout.rank && continue
                    isempty(indices[rank]) && continue
                    push!(packets[rank], (s, HybridPlasmaPIC._subset_particles(ps, indices[rank])))
                end
            end
            bytes = [
                isempty(packet) ? UInt8[] : HybridPlasmaPIC._serialize_mpi_payload(packet) for
                packet in packets
            ]
            result, bytes
        end
    chunks = HybridPlasmaPIC._mpi_alltoallv_bytes(outgoing, g.fft.comm)
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "migration assembly") do
        for chunk in chunks
            isempty(chunk) && continue
            incoming = HybridPlasmaPIC._deserialize_mpi_payload(Vector{Tuple{Int,P}}, chunk)
            for (s, ps) in incoming
                1 <= s <= length(prepared) || throw(ArgumentError("invalid incoming species index"))
                destination = prepared[s]
                isequal((ps.q, ps.m), (destination.q, destination.m)) ||
                    throw(ArgumentError("incoming particle charge/mass differ from the species"))
                _pencil_pic_population(ps)
                HybridPlasmaPIC.append_particles!(destination, ps)
            end
        end
        foreach(_pencil_pic_population, prepared)
        _pencil_owned_particles(prepared, g, layout)
        prepared
    end
end

function _pencil_gather(components, species, g, layout, context)
    T = eltype(parent(first(components)))
    requests, outgoing =
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "CIC request preparation") do
            sets = [Set{_PencilNode}() for _ in layout.by_rank]
            stamp = CartesianIndices((2, 2, 2))
            for ps in species, p in eachindex(ps.weight)
                st = _pencil_cic_stencil(ps, g, p)
                for offset in stamp
                    node = _pencil_cic_node(st, Tuple(offset), g.n)
                    rank = _pencil_node_owner(layout, node)
                    rank == layout.rank || push!(sets[rank], node)
                end
            end
            nodes = [sort!(collect(set)) for set in sets]
            nodes,
            [isempty(row) ? UInt8[] : HybridPlasmaPIC._serialize_mpi_payload(row) for row in nodes]
        end
    chunks = HybridPlasmaPIC._mpi_alltoallv_bytes(outgoing, g.fft.comm)
    replies = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        context,
        "CIC field response preparation",
    ) do
        fields = map(PencilArrays.global_view, components)
        result = [UInt8[] for _ in layout.by_rank]
        for rank in eachindex(chunks)
            isempty(chunks[rank]) && continue
            nodes = HybridPlasmaPIC._deserialize_mpi_payload(Vector{_PencilNode}, chunks[rank])
            values = NTuple{3,T}[]
            for node in nodes
                _pencil_contains(layout.by_rank[layout.rank], node) ||
                    throw(ArgumentError("CIC field request was delivered to a nonowner"))
                push!(values, ntuple(c -> fields[c][node...], 3))
            end
            result[rank] =
                HybridPlasmaPIC._serialize_mpi_payload((nodes = nodes, values = values))
        end
        result
    end
    returned = HybridPlasmaPIC._mpi_alltoallv_bytes(replies, g.fft.comm)
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "CIC gather") do
        ghosts = Dict{_PencilNode,NTuple{3,T}}()
        for rank in eachindex(returned)
            if isempty(requests[rank])
                isempty(returned[rank]) || throw(ArgumentError("unexpected CIC field response"))
                continue
            end
            row = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, returned[rank])
            keys(row) == (:nodes, :values) &&
                row.nodes isa Vector{_PencilNode} &&
                row.values isa Vector{NTuple{3,T}} &&
                isequal(row.nodes, requests[rank]) &&
                length(row.values) == length(row.nodes) ||
                throw(ArgumentError("CIC field response does not match its request"))
            for (node, values) in zip(row.nodes, row.values)
                all(isfinite, values) || throw(ArgumentError("CIC field response must be finite"))
                ghosts[node] = values
            end
        end
        fields = map(PencilArrays.global_view, components)
        stamp = CartesianIndices((2, 2, 2))
        gathers = NTuple{3,Vector{T}}[]
        field_at(node) =
            _pencil_contains(layout.by_rank[layout.rank], node) ?
            ntuple(c -> fields[c][node...], 3) : ghosts[node]
        for ps in species
            output = ntuple(_ -> Vector{T}(undef, length(ps.weight)), 3)
            for p in eachindex(ps.weight)
                values = HybridPlasmaPIC._gather_vector_stencil(
                    field_at,
                    _pencil_cic_stencil(ps, g, p),
                    g.n,
                    stamp,
                    T,
                )
                all(isfinite, values) || throw(
                    ArgumentError(
                        "distributed PIC gather exceeds the finite stored-field boundary",
                    ),
                )
                for c = 1:3
                    output[c][p] = values[c]
                end
            end
            push!(gathers, output)
        end
        gathers
    end
end

# Exact contributions cross ranks before division/rounding, so no partial
# density or charge product needs to fit the working precision.
function _pencil_deposit_cic!(
    out,
    species,
    g,
    layout,
    background,
    component::Int,
    context;
    midpoints = nothing,
)
    T = eltype(parent(out))
    sums, outgoing =
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "CIC moment preparation") do
            0 <= component <= 3 || throw(ArgumentError("CIC moment component must be in 0:3"))
            midpoints === nothing ||
                length(midpoints) == length(species) ||
                throw(ArgumentError("CIC midpoint rows must match species"))
            exact = fill!(similar(out, _PencilCharge), zero(_PencilCharge))
            owned = PencilArrays.global_view(exact)
            packets = [Dict{_PencilNode,_PencilCharge}() for _ in layout.by_rank]
            for (s, ps) in enumerate(species)
                charge = HybridPlasmaPIC._exact_rational(ps.q)
                factors = component == 0 ? (ps.weight,) : (ps.weight, ps.v[component])
                midpoints === nothing ||
                    HybridPlasmaPIC._swap_positions_with_midpoints!(ps, midpoints[s])
                try
                    HybridPlasmaPIC._foreach_exact_moment_contribution(
                        ps,
                        factors,
                        g,
                        HybridPlasmaPIC.CIC(),
                    ) do node, weight
                        value = charge * weight
                        if !iszero(value)
                            rank = _pencil_node_owner(layout, node)
                            if rank == layout.rank
                                owned[node...] += value
                            else
                                packets[rank][node] =
                                    get(packets[rank], node, zero(_PencilCharge)) + value
                            end
                        end
                    end
                finally
                    midpoints === nothing ||
                        HybridPlasmaPIC._swap_positions_with_midpoints!(ps, midpoints[s])
                end
            end
            bytes = [
                isempty(packet) ? UInt8[] :
                HybridPlasmaPIC._serialize_mpi_payload(sort!(collect(packet); by = first)) for
                packet in packets
            ]
            exact, bytes
        end
    chunks = HybridPlasmaPIC._mpi_alltoallv_bytes(outgoing, g.fft.comm)
    return HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        context,
        "CIC moment assembly and rounding",
    ) do
        owned = PencilArrays.global_view(sums)
        for chunk in chunks
            isempty(chunk) && continue
            packet = HybridPlasmaPIC._deserialize_mpi_payload(
                Vector{Pair{_PencilNode,_PencilCharge}},
                chunk,
            )
            for (node, value) in packet
                _pencil_contains(layout.by_rank[layout.rank], node) ||
                    throw(ArgumentError("CIC moment contribution was delivered to a nonowner"))
                owned[node...] += value
            end
        end
        Q = HybridPlasmaPIC._exact_rational
        volume, offset = prod(Q, g.dx), Q(background)
        for I in eachindex(parent(out), parent(sums))
            value = HybridPlasmaPIC._from_exact(T, offset + parent(sums)[I] / volume)
            isfinite(value) || throw(ArgumentError("complete distributed CIC moment is not finite"))
            parent(out)[I] = value
        end
        nothing
    end
end
