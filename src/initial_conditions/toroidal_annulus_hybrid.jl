# Setup helpers use the same represented point/node geometry as the coupled RHS.
function _tah_initial_scalar(name, value, ::Type{T}; positive = false) where {T}
    value isa Real && isfinite(value) ||
        throw(ArgumentError("$name must be a finite real scalar"))
    positive && !(value > 0) &&
        throw(ArgumentError("$name must be positive"))
    result = _require_finite_real(name, value, T)
    !iszero(value) && iszero(result) &&
        throw(ArgumentError("$name must not disappear during conversion to $T"))
    positive && !(result > zero(T)) &&
        throw(ArgumentError("$name must be positive in $T"))
    return result
end

"""
    load_toroidal_annulus_cold!(population, mesh; number_density,
        velocity=(0,0,0), particles_per_cell=(1,1,1))

Prepare a deterministic cold load in the actual logical cell interiors, with
Cartesian velocity and physical number weight `n0*J(point)*cell_volume/prod(M)`.
Requires `(Nr-1)*Nu*Nv*prod(M)` existing ordinary 3D host particles. All input,
alias, coordinate and weight checks precede publication; IDs and tags are kept.
Retained charge and mass must be finite positive values.
Call before binding a simulation. Deposited density and equilibrium residuals
are measurements, not values assigned by this loader.
"""
function load_toroidal_annulus_cold!(
    p::ParticleSet{3,T},
    mesh::ToroidalAnnulusGrid{T};
    number_density,
    velocity = (0, 0, 0),
    particles_per_cell = (1, 1, 1),
) where {T}
    _tah_geometry(mesh)
    _tah_particle_storage(p, T)
    isfinite(p.q) && p.q > zero(T) && isfinite(p.m) && p.m > zero(T) ||
        throw(ArgumentError("annular loading requires finite positive retained charge and mass"))
    particles_per_cell isa Union{Tuple,AbstractVector} && length(particles_per_cell) == 3 ||
        throw(DimensionMismatch("annular particles_per_cell requires three counts"))
    particles_per_cell isa AbstractVector && axes(particles_per_cell, 1) != Base.OneTo(3) &&
        throw(ArgumentError("annular particles_per_cell must use one-based axes"))
    M = ntuple(d -> _require_positive_intlike("particles_per_cell[$d]", particles_per_cell[d]), 3)
    cells = (mesh.n[1] - 1, mesh.n[2], mesh.n[3])
    markers_per_cell = prod(BigInt, M)
    required = prod(BigInt, cells) * markers_per_cell
    required == nparticles(p) || throw(
        DimensionMismatch("annular cold load requires (Nr-1)*Nu*Nv*prod(particles_per_cell) existing particles"),
    )
    density = _tah_initial_scalar("annular number_density", number_density, T; positive = true)
    velocity isa Union{Tuple,AbstractVector} && length(velocity) == 3 ||
        throw(DimensionMismatch("annular Cartesian loading velocity requires three components"))
    velocity isa AbstractVector && axes(velocity, 1) != Base.OneTo(3) &&
        throw(ArgumentError("annular Cartesian loading velocity must use one-based axes"))
    drift = ntuple(c -> _tah_initial_scalar("annular loading velocity[$c]", velocity[c], T), 3)

    # Include unchanged provenance leaves: sharing one can corrupt preserved IDs/tags.
    arrays = _tah_particle_arrays(p)
    geometry = _tah_mesh_arrays(mesh)
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("annular loading particle leaves must be independent"))
        end
        any(a -> _particle_arrays_mightalias(arrays[i], a), geometry) &&
            throw(ArgumentError("annular loading leaves overlap geometry, operator or retained plan storage"))
    end

    Q = _exact_rational
    positions = ntuple(_ -> Vector{T}(undef, nparticles(p)), 3)
    weights = Vector{T}(undef, nparticles(p))
    marker = 0
    for cell in CartesianIndices(cells)
        lo = ntuple(d -> Q(mesh.op.coordinates[d][cell[d]]), 3)
        hi = ntuple(3) do d
            (d == 1 || cell[d] < mesh.n[d]) ?
                Q(mesh.op.coordinates[d][cell[d]+1]) : Q(mesh.L[d])
        end
        widths = ntuple(d -> hi[d] - lo[d], 3)
        all(>(zero(_CompositeQ)), widths) ||
            throw(ArgumentError("annular logical cell edges must remain distinct"))
        logical_volume = prod(widths)
        for offset in CartesianIndices(M)
            marker += 1
            point = ntuple(3) do d
                fraction = (BigInt(2) * offset[d] - 1) // (BigInt(2) * M[d])
                coordinate = _tah_round(T, lo[d] + widths[d] * fraction,
                    "annular quiet-load coordinate"; preserve_nonzero = true)
                lower = lo[d] + widths[d] * (BigInt(offset[d]) - 1) / M[d]
                upper = lo[d] + widths[d] * BigInt(offset[d]) / M[d]
                lower < Q(coordinate) < upper ||
                    throw(ArgumentError("annular quiet-load marker must remain strictly inside its logical subcell"))
                coordinate
            end
            geo = _tah_point_geometry(mesh, point)
            weight = _tah_round(T, Q(density) * geo.J * logical_volume / markers_per_cell,
                "annular physical number weight"; preserve_nonzero = true)
            weight > zero(T) || throw(ArgumentError("annular physical number weight must be positive"))
            weights[marker] = weight
            for d = 1:3
                positions[d][marker] = point[d]
            end
        end
    end

    # Every expected rejection and allocation precedes fixed-length ordinary stores.
    for i in eachindex(p.weight)
        @inbounds p.weight[i] = weights[i]
        for c = 1:3
            @inbounds p.x[c][i] = positions[c][i]
            @inbounds p.v[c][i] = drift[c]
        end
    end
    return p
end

"""
    toroidal_annulus_vacuum_flux(mesh; strength=one(T))

Return detached primary flux arrays for the analytic profile
`B=(strength/R)*e_phi`, using the actual represented nodal metric factors.
Its only nonzero component is `beta[3]=J*strength/(h[3]*R)`. This seed does not
assign measured zero current, force balance or a discrete equilibrium.
"""
function toroidal_annulus_vacuum_flux(mesh::ToroidalAnnulusGrid{T}; strength = one(T)) where {T}
    _tah_geometry(mesh)
    C = _tah_initial_scalar("annular vacuum strength", strength, T)
    Q = _exact_rational
    beta = ntuple(_ -> zeros(T, mesh.n), 3)
    for I in CartesianIndices(mesh.n)
        h, J = _tah_metric(mesh, I)
        beta[3][I] = _tah_round(T, J * Q(C) / (h[3] * Q(mesh.R[I[1], I[2]])),
            "annular vacuum primary flux"; preserve_nonzero = true)
    end
    return beta
end
