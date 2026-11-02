struct DistributedMaxwellVacuum{T,G,F,R}
    g::G
    fields::F
    c::T
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    _candidate::F
    _curl::NTuple{3,R}
    _cfl::T
end

function _distributed_maxwell_mode(::Type{T}) where {T}
    Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() || throw(
        ArgumentError("distributed vacuum Maxwell requires nearest rounding and gradual underflow"),
    )
    return nothing
end

function _distributed_maxwell_cfl(g::DistributedFourierGrid{T}, c::T) where {T}
    Q = HybridPlasmaPIC._exact_rational
    frequency_squared = Q(c)^2 * sum(gain -> Q(gain)^2, g.gain)
    iszero(frequency_squared) && return T(Inf)
    return HybridPlasmaPIC._sqrt_from_exact(
        T,
        4 / frequency_squared,
        Int128(0),
        RoundDown,
        false,
    )
end

"""
    distributed_maxwell_vacuum(grid; c=1)

Collectively allocate synchronized, initially zero E/B fields on a periodic 3D
`distributed_fourier_grid`. The speed converts to the grid's Float32/Float64
precision and must agree across ranks. Populate owned cells through
`PencilArrays.global_view(model.fields.E[c])` and the corresponding B components.
The grid, metadata and caller-owned communicator retain their existing lifetime
and one-caller policies; no global field is replicated or gathered.

`step_maxwell!` uses three canonical distributed curls and magnetic/electric/
magnetic KDK stages. Each scalar affine stage retains the converted duration and
speed exactly before one final rounding. This intentionally defines a new map,
not ordinary EMPIC/Yee bit identity; exact scalar arithmetic has unmeasured cost.
The source is the represented distributed curl, whose FFT/range limits still
apply. There is no original-field derivative range recovery. Matching-layout
checkpoint overloads are described in `docs/DISTRIBUTED_MAXWELL_CHECKPOINT.md`.
Nearest rounding and gradual underflow are required. Concurrent use of the model
or its grid is unsupported. See `docs/DISTRIBUTED_MAXWELL_IMPLEMENTATION.md`.
"""
function HybridPlasmaPIC.distributed_maxwell_vacuum(g::DistributedFourierGrid{T}; c = 1) where {T}
    context = "Distributed vacuum Maxwell"
    _, _, speed, _, _ = _agree_configuration(g.fft.comm, context) do
        _distributed_maxwell_mode(T)
        c isa Real || throw(ArgumentError("c must be real"))
        g.n, g.L, HybridPlasmaPIC._require_finite_positive_real("c", c, T), T, g.dx
    end
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "model allocation") do
        vector() = ntuple(_ -> HybridPlasmaPIC.distributed_field(g), 3)
        fields = (E = vector(), B = vector())
        candidate = (E = vector(), B = vector())
        _distributed_maxwell_finite_fields(fields)
        DistributedMaxwellVacuum(
            g,
            fields,
            speed,
            Ref(zero(T)),
            Ref(0),
            candidate,
            vector(),
            _distributed_maxwell_cfl(g, speed),
        )
    end
end

"""
    maxwell_cfl(model::DistributedMaxwellVacuum)

The downward-rounded periodic spectral KDK bound `2/(c*sqrt(sum(gain.^2)))`,
using the grid's derivative-resolved, odd-Nyquist-zeroed axis gains. Every step
must be strictly smaller. A completely derivative-null grid returns Inf.
This constant-step bound is not a variable-step stability or accuracy guarantee.
"""
HybridPlasmaPIC.maxwell_cfl(model::DistributedMaxwellVacuum) = model._cfl

function _distributed_maxwell_storage(model::DistributedMaxwellVacuum)
    arrays = (
        model.fields.E...,
        model.fields.B...,
        model._candidate.E...,
        model._candidate.B...,
        model._curl...,
    )
    for i in eachindex(arrays)
        _require_physical_field("Maxwell storage", arrays[i], model.g, i)
        for j = 1:(i-1)
            Base.mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("distributed Maxwell physical and private fields must be independent"))
        end
    end
    return nothing
end

function _distributed_maxwell_finite_fields(fields)
    _finite_vector_magnitude(fields.E)
    _finite_vector_magnitude(fields.B)
    return nothing
end

function _distributed_maxwell_state_agreement(model::DistributedMaxwellVacuum{T}, context) where {T}
    return _agree_configuration(model.g.fft.comm, context) do
        _distributed_maxwell_mode(T)
        HybridPlasmaPIC._require_finite_positive_real("c", model.c, T)
        isfinite(model.time[]) && model.time[] >= zero(T) ||
            throw(ArgumentError("distributed Maxwell time must be finite and nonnegative"))
        model.step[] >= 0 ||
            throw(ArgumentError("distributed Maxwell step must be nonnegative"))
        model.g.n, model.g.L, model.g.dx, model.c, model.time[], model.step[], T
    end
end

# The existing complete Maxwell expression owns every local affine update.
# Its zero-current coefficients select either sign without materializing h*c²
# or a rounded half step. Zero curls leave the old signed zero unchanged.
function _distributed_maxwell_affine!(out, old, curl, coefficients)
    T = eltype(first(out))
    @inbounds for c = 1:3
        destination, source, slope = parent(out[c]), parent(old[c]), parent(curl[c])
        for I in eachindex(destination, source, slope)
            derivative = slope[I]
            isfinite(derivative) || throw(
                ArgumentError(
                    "distributed Maxwell requires a finite materialized curl; " *
                    "completion through a later timestep factor is not supported",
                ),
            )
            value = if HybridPlasmaPIC._particle_stored_zero(derivative)
                source[I]
            else
                HybridPlasmaPIC._em_exact_electric_update(
                    source[I],
                    derivative,
                    zero(T),
                    coefficients,
                )
            end
            isfinite(value) ||
                throw(ArgumentError("complete distributed Maxwell stage is not finite in $T"))
            destination[I] = value
        end
    end
    return nothing
end

"""
    step_maxwell!(model::DistributedMaxwellVacuum, dt)

Collectively advance synchronized periodic spectral vacuum fields by one positive
strictly sub-CFL step. Converted dt, speed and clocks must agree on all ranks.
Mode, storage, aliases and finite physical fields validate before staging writes.
The three curls and complete scalar stages write only private fields; all six
physical arrays and clocks commit after the final collective success. Recoverable
local failure preserves physical state on every rank; private scratch may change.
MPI process/transport loss and fatal resource failures retain the operator's
exclusions. No projection removes an initial divergence residual.
"""
function HybridPlasmaPIC.step_maxwell!(model::DistributedMaxwellVacuum{T}, dt) where {T}
    g = model.g
    context = "Distributed vacuum Maxwell"
    parameters = _agree_configuration(g.fft.comm, context) do
        _distributed_maxwell_mode(T)
        HybridPlasmaPIC._require_finite_positive_real("c", model.c, T)
        dt isa Real || throw(ArgumentError("dt must be real"))
        h = HybridPlasmaPIC._require_finite_positive_real("dt", dt, T)
        h < HybridPlasmaPIC.maxwell_cfl(model) ||
            throw(ArgumentError("dt must be strictly below the distributed spectral Maxwell CFL limit"))
        next_time, next_step = HybridPlasmaPIC._validated_step_clock(model.time[], model.step[], h)
        next_time > model.time[] ||
            throw(ArgumentError("distributed Maxwell time must advance finitely"))
        h, model.c, model.time[], model.step[], next_time, next_step, g.n, g.L, g.dx, T
    end
    coefficients = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "step preflight") do
        _distributed_maxwell_storage(model)
        _distributed_maxwell_finite_fields(model.fields)
        Q = HybridPlasmaPIC._exact_rational
        h, speed = Q(parameters[1]), Q(model.c)
        z = zero(h)
        (magnetic = (z, -h / 2), electric = (z, h * speed^2))
    end
    original, candidate, work = model.fields, model._candidate, model._curl
    SpectralOperators.curl!(work, original.E, g)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "first magnetic half stage") do
        _distributed_maxwell_affine!(candidate.B, original.B, work, coefficients.magnetic)
    end
    SpectralOperators.curl!(work, candidate.B, g)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "electric stage") do
        _distributed_maxwell_affine!(candidate.E, original.E, work, coefficients.electric)
    end
    SpectralOperators.curl!(work, candidate.E, g)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "final magnetic half stage") do
        _distributed_maxwell_affine!(candidate.B, candidate.B, work, coefficients.magnetic)
    end
    for c = 1:3
        copyto!(parent(original.E[c]), parent(candidate.E[c]))
        copyto!(parent(original.B[c]), parent(candidate.B[c]))
    end
    model.time[] = parameters[5]
    model.step[] = parameters[6]
    return model
end

# Only exact scalar energies cross ranks, never physical fields. Prepare the
# serialized payload on every rank before any broadcast; guard allocation and
# decode/add before the next collective, retaining local failures through the
# same owner as the FFT orchestration. This is not an object-reduction API.
function _distributed_maxwell_energy(model::DistributedMaxwellVacuum{T}) where {T}
    comm = model.g.fft.comm
    context = "Distributed Maxwell energy"
    total = _pencil_exact_scalar_sum(
        comm,
        context,
        "local energy preparation",
        "exact energy accumulation",
    ) do
        HybridPlasmaPIC._field_energy_rational(
            map(parent, model.fields.E),
            map(parent, model.fields.B),
            HybridPlasmaPIC._exact_rational(model.c)^2,
            model.g.dx,
        )
    end
    return HybridPlasmaPIC._mpi_collective_stage(comm, context, "energy rounding") do
        HybridPlasmaPIC._from_exact(T, total)
    end
end

"""
    em_field_energy(model::DistributedMaxwellVacuum)

Collectively return the normalized global energy `sum(E²+c²B²)*prod(dx)/2`.
Owned-cell energy and the sum of scalar rank contributions remain exact until
one final rounding to the field precision; a nonrepresentable total returns Inf.
Only scalar serialized contributions are exchanged. No field is gathered.
"""
function HybridPlasmaPIC.em_field_energy(model::DistributedMaxwellVacuum)
    context = "Distributed Maxwell energy"
    _distributed_maxwell_state_agreement(model, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "storage validation") do
        _distributed_maxwell_storage(model)
        _distributed_maxwell_finite_fields(model.fields)
    end
    return _distributed_maxwell_energy(model)
end

"""
    distributed_maxwell_diagnostics(model)

Collectively return synchronized time/step, global energy and maximum absolute
electric/magnetic spectral divergence. The two existing divergence operators use
private fields; no physical state is changed, projected or gathered. Every rank
must call collective model operations in the same order. Runtime, compiler and
resource acceptance of this source implementation remain deferred.
"""
function HybridPlasmaPIC.distributed_maxwell_diagnostics(model::DistributedMaxwellVacuum)
    energy = HybridPlasmaPIC.em_field_energy(model)
    g = model.g
    SpectralOperators.divergence!(model._curl[1], model.fields.E, g)
    SpectralOperators.divergence!(model._curl[2], model.fields.B, g)
    maxima = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        "Distributed Maxwell diagnostics",
        "divergence reduction preparation",
    ) do
        maximum(abs, parent(model._curl[1])), maximum(abs, parent(model._curl[2]))
    end
    divE = MPI.Allreduce(maxima[1], MPI.MAX, g.fft.comm)
    divB = MPI.Allreduce(maxima[2], MPI.MAX, g.fft.comm)
    return (
        time = model.time[],
        step = model.step[],
        energy = energy,
        max_divE = divE,
        max_divB = divB,
    )
end
