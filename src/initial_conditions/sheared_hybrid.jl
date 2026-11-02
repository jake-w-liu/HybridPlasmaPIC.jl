"""
    load_sheared_uniform!(population, mesh; number_density, velocity=(0,0,0), particles_per_axis=1)

Prepare a real cold quiet load in logical cells, with physical Cartesian velocity
and number weights n0*dV/M^3. Requires Ng*M^3 existing ordinary particles. Call
before binding a simulation; deposited uniformity is a measured numerical result.
"""
function load_sheared_uniform!(p::ParticleSet{3,T}, mesh::PeriodicShearGrid{T};
    number_density, velocity = (0, 0, 0), particles_per_axis = 1) where {T}
    _boris_accuracy_mode() || throw(ArgumentError("sheared loading requires nearest rounding and gradual underflow"))
    _shear_geometry(mesh)
    _sh_particle(p; values = false)
    g = mesh.logical
    M = _require_positive_intlike("particles_per_axis", particles_per_axis)
    required = BigInt(length(g.cbuf)) * BigInt(M)^3
    required == nparticles(p) || throw(DimensionMismatch("uniform shear load requires Ng*M^3 existing particles"))
    density = _require_finite_positive_real("number_density", number_density, T)
    length(velocity) == 3 || throw(DimensionMismatch("physical loading velocity requires three components"))
    drift = ntuple(c -> _require_finite_real("loading velocity", velocity[c], T), 3)
    arrays = _sh_particle_arrays(p)
    geometry = (mesh._s, g.cbuf, g.tbuf, g.abuf, g.sbuf, g.midx..., g.kfull..., g.kvec..., g.ik...)
    for i in eachindex(arrays)
        for j in 1:(i - 1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(ArgumentError("loading particle leaves must be independent"))
        end
        any(a -> _particle_arrays_mightalias(arrays[i], a), geometry) && throw(ArgumentError("loading leaves overlap geometry or Fourier work"))
    end
    _visit_checkpoint_objects((g.plan, g.iplan)) do value
        value isa AbstractArray && any(a -> _particle_arrays_mightalias(a, value), arrays) &&
            throw(ArgumentError("loading leaves overlap retained FFT plan storage"))
    end
    Q = _exact_rational
    weight = _shear_round(T, Q(density) * prod(Q, g.dx) / BigInt(M)^3,
        "uniform shear number weight"; preserve_nonzero = true)
    positions = ntuple(_ -> Vector{T}(undef, nparticles(p)), 3)
    marker = 0
    for cell in CartesianIndices(g.n), offset in CartesianIndices((M, M, M))
        marker += 1
        for d in 1:3
            fraction = (BigInt(2) * offset[d] - 1 - M) // (BigInt(2) * M)
            x = _shear_round(T, Q(g.dx[d]) * (cell[d] - 1 + fraction), "quiet-load coordinate")
            positions[d][marker] = _periodic_coordinate(x, zero(T), g.L[d], g.L[d])
        end
    end
    # All validation, allocation and conversions precede the first public store.
    for i in eachindex(p.weight)
        @inbounds p.weight[i] = weight
        for c in 1:3
            @inbounds p.x[c][i] = positions[c][i]
            @inbounds p.v[c][i] = drift[c]
        end
    end
    return p
end
