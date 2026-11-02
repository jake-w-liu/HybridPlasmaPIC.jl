# The inverse-gradient component gain is bounded by 1/kmin on every resolved
# mode. Round the bound upward; a rounded-down reciprocal is not a safety bound.
function _pencil_poisson_gain(g::DistributedFourierGrid{T}) where {T}
    minimum_wave = T(Inf)
    for axis in g.kvec, wave in axis
        iszero(wave) || (minimum_wave = min(minimum_wave, abs(wave)))
    end
    isinf(minimum_wave) && return zero(T)
    Q = HybridPlasmaPIC._exact_rational
    exact = inv(Q(minimum_wave))
    gain = HybridPlasmaPIC._from_exact(T, exact)
    isfinite(gain) && Q(gain) < exact && (gain = nextfloat(gain))
    isfinite(gain) ||
        throw(ArgumentError("distributed Poisson requires a finite inverse-wavenumber bound"))
    return gain
end

"""
    distributed_poisson_field!(E, rho, grid; neutrality=:require)

Collectively solve the periodic 3D, derivative-resolved Poisson field from finite
stored charge density on the grid's exact physical pencil. Uses the existing
serial Poisson mode owner and distributed transforms; derivative-null DC/Nyquist
modes have zero field. This is a static field prerequisite, not a particle solver.

`:require` demands an exactly zero sum of represented input densities across
ranks. `:neutralize` supplies the explicit uniform compensator. Return
`(E, source_charge, neutralizing_charge)` with exact extensive charges using stored
cell volume. Density changes only when explicitly used as an output; no source
is reconstructed from particles.

Output components must be independent; one may alias rho. All three outputs
publish together after collective success. Nearest rounding, gradual underflow
and the existing finite FFT/scale boundaries apply. No field is gathered. Runtime,
compiler and resource acceptance is deferred during this implementation pass.
"""
function HybridPlasmaPIC.distributed_poisson_field!(
    out,
    rho,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R};
    neutrality = :require,
) where {T,R}
    context = "Distributed Poisson field"
    _agree_configuration(g.fft.comm, context) do
        Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() || throw(
            ArgumentError("distributed Poisson requires nearest rounding and gradual underflow"),
        )
        (neutrality === :require || neutrality === :neutralize) ||
            throw(ArgumentError("neutrality must be :require or :neutralize"))
        T, g.n, g.L, g.dx, neutrality
    end
    magnitude = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "storage preflight") do
        _require_physical_field("charge density", rho, g)
        _require_vector_field("electric output", out, g)
        _require_distinct_components(out, "Poisson output components must be independent")
        _finite_magnitude(rho)
    end
    total = _pencil_exact_scalar_sum(
        g.fft.comm,
        context,
        "local charge preparation",
        "exact charge accumulation",
    ) do
        HybridPlasmaPIC._exact_float_sum(parent(rho))
    end
    gain, result = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "charge admission") do
        neutrality === :require &&
            !iszero(total) &&
            throw(
                ArgumentError(
                    "periodic Poisson requires exact neutral charge or neutrality=:neutralize",
                ),
            )
        charge = total * prod(HybridPlasmaPIC._exact_rational, g.dx)
        report = (
            E = out,
            source_charge = charge,
            neutralizing_charge = neutrality === :neutralize ? -charge : zero(charge),
        )
        _pencil_poisson_gain(g), report
    end
    staging = _staging_buffers(g)
    if iszero(gain)
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "derivative-null field") do
            foreach(field -> fill!(parent(field), zero(T)), staging)
        end
    else
        scale = _agreed_scale(g, magnitude, gain)
        _forward_field!(g.spectral, rho, g, scale, context)
        buffers = (g.spectral, g.spectral_work...)
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "Poisson mode solve") do
            modes = PencilArrays.global_view(g.spectral)
            destination = map(PencilArrays.global_view, buffers)
            for I in CartesianIndices(modes)
                wave = ntuple(d -> g.kvec[d][I[d]], 3)
                values = HybridPlasmaPIC._spectral_poisson_mode(wave, modes[I])
                for c = 1:3
                    destination[c][I] = values[c]
                end
            end
        end
        for c = 1:3
            _inverse_field!(staging[c], buffers[c], g, scale, context)
        end
    end
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "complete output validation") do
        _finite_vector_magnitude(staging)
    end
    _publish_components!(out::NTuple{3,R}, staging)
    return result
end
