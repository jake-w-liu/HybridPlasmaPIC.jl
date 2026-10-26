const SpectralOperators = HybridPlasmaPIC.SpectralOperators

struct DistributedFourierGrid{T,F,A,B,R}
    n::NTuple{3,Int}
    L::NTuple{3,T}
    dx::NTuple{3,T}
    midx::NTuple{3,Vector{Int}}
    kfull::NTuple{3,Vector{T}}
    ik::NTuple{3,Vector{Complex{T}}}
    kvec::NTuple{3,Vector{T}}
    gain::NTuple{3,T}
    count::Int
    fft::F
    physical::A
    spectral::B
    staging::R
    # Vector operators keep up to three spectral components live and stage all
    # three real outputs before publication.
    spectral_work::NTuple{2,B}
    staging_work::NTuple{2,R}
end

"""
    distributed_fourier_grid(n, L; comm=MPI.COMM_WORLD)

Collectively construct a periodic 3D Fourier workspace with rank-owned fields.
Every rank must supply the same three grid sizes and positive finite lengths.
Lengths have one common Float32 or Float64 type, which selects field precision.
The workspace belongs to one caller at a time; its metadata vectors are read-only.
The communicator remains caller-owned. Each rank owns one complex physical,
three complex spectral and three real publication arrays over its cells.
This workspace provides derivatives, gradient, divergence, curl, Laplacian,
divergence-free projection and spectral filters; it does not provide particle
halos, distributed stepping or checkpoint storage.
"""
function HybridPlasmaPIC.distributed_fourier_grid(n, L; comm = MPI.COMM_WORLD)
    HybridPlasmaPIC.ensure_mpi_initialized!()
    nn, lengths = _agree_configuration(comm, "Distributed Fourier grid") do
        nn = _checked_fft_shape(n)
        (L isa NTuple{3,Float32} || L isa NTuple{3,Float64}) || throw(
            ArgumentError(
                "distributed Fourier lengths must be a three-tuple of Float32 or Float64",
            ),
        )
        all(x -> isfinite(x) && x > zero(x), L) ||
            throw(ArgumentError("distributed Fourier lengths must be positive and finite"))
        nn, L
    end
    metadata, count = HybridPlasmaPIC._mpi_collective_stage(
        comm,
        "Distributed Fourier grid",
        "metadata construction",
    ) do
        count = Base.checked_mul(Base.checked_mul(nn[1], nn[2]), nn[3])
        SpectralOperators._fourier_metadata(nn, lengths), count
    end
    T = eltype(lengths)
    fft = HybridPlasmaPIC.distributed_fft_plan(nn; comm, real_type = T)
    return HybridPlasmaPIC._mpi_collective_stage(
        comm,
        "Distributed Fourier grid",
        "workspace allocation",
    ) do
        physical = HybridPlasmaPIC.distributed_fft_input(fft)
        spectral = HybridPlasmaPIC.distributed_fft_output(fft)
        staging = similar(physical, T)
        spectral_work = (
            HybridPlasmaPIC.distributed_fft_output(fft),
            HybridPlasmaPIC.distributed_fft_output(fft),
        )
        staging_work = (similar(physical, T), similar(physical, T))
        gain = ntuple(d -> maximum(abs, metadata.kvec[d]), 3)
        DistributedFourierGrid(
            metadata.n,
            metadata.L,
            metadata.dx,
            metadata.midx,
            metadata.kfull,
            metadata.ik,
            metadata.kvec,
            gain,
            count,
            fft,
            physical,
            spectral,
            staging,
            spectral_work,
            staging_work,
        )
    end
end

"""
    distributed_field(grid)

Allocate a zeroed real CPU PencilArray on `grid`'s exact physical pencil.
Use `PencilArrays.global_view(field)` to access global logical cell indices.
Allocation is local; every collective derivative still requires all grid ranks.
"""
function HybridPlasmaPIC.distributed_field(g::DistributedFourierGrid)
    return fill!(similar(g.staging), zero(eltype(g.staging)))
end

# Component labels are formatted only when an error is reported.
_operand_label(name, component::Int) = component == 0 ? name : "$name[$component]"

# Homogeneous loops keep this check allocation-free.
function _aliases_workspace(a, g::DistributedFourierGrid)
    Base.mightalias(a, g.physical) && return true
    for workspace in (g.spectral, g.spectral_work...)
        Base.mightalias(a, workspace) && return true
    end
    for workspace in (g.staging, g.staging_work...)
        Base.mightalias(a, workspace) && return true
    end
    # Thin-pencil fields can share read-only metadata storage, including
    # detectable Array views with a different element type.
    for metadata in (g.kfull..., g.kvec...)
        Base.mightalias(a, metadata) && return true
    end
    for metadata in g.ik
        Base.mightalias(a, metadata) && return true
    end
    for metadata in g.midx
        Base.mightalias(a, metadata) && return true
    end
    return false
end

function _require_physical_field(
    name,
    a,
    g::DistributedFourierGrid{T},
    component::Int = 0,
) where {T}
    a isa PencilArrays.PencilArray{T,3} || throw(
        ArgumentError(
            "$(_operand_label(name, component)) must be a three-dimensional PencilArray with element type $T",
        ),
    )
    parent(a) isa Array{T,3} ||
        throw(ArgumentError("$(_operand_label(name, component)) must use CPU Array storage"))
    PencilArrays.pencil(a) === PencilArrays.pencil(g.physical) || throw(
        ArgumentError("$(_operand_label(name, component)) must use this grid's physical pencil"),
    )
    axes(a) == axes(g.physical) || throw(
        DimensionMismatch(
            "$(_operand_label(name, component)) axes do not match the grid's owned cells",
        ),
    )
    _aliases_workspace(a, g) && throw(
        ArgumentError(
            "$(_operand_label(name, component)) must not alias private Fourier workspace storage",
        ),
    )
    return nothing
end

# Vector fields are three-tuples of rank-owned fields on the grid's pencil.
function _require_vector_field(name, v, g::DistributedFourierGrid)
    v isa Tuple && length(v) == 3 ||
        throw(ArgumentError("$name must be a three-tuple of distributed fields"))
    for c = 1:3
        _require_physical_field(name, v[c], g, c)
    end
    return nothing
end

function _require_distinct_components(v, message)
    for c = 1:2, d = (c+1):3
        Base.mightalias(v[c], v[d]) && throw(ArgumentError(message))
    end
    return nothing
end

function _finite_magnitude(field)
    magnitude = maximum(abs, parent(field))
    isfinite(magnitude) || throw(ArgumentError("spectral operator input values must be finite"))
    return magnitude
end

# Like the serial vector operators, check each used component before combining.
function _finite_vector_magnitude(v)
    magnitude = _finite_magnitude(v[1])
    magnitude = max(magnitude, _finite_magnitude(v[2]))
    return max(magnitude, _finite_magnitude(v[3]))
end

"""
    deriv!(out, field, grid::DistributedFourierGrid, axis)

Collectively differentiate a real rank-owned field along the agreed global axis.
Input and output may alias each other. Both must use the grid's physical pencil,
precision and CPU Array storage, and must not alias its private workspaces.
Recoverable local failures are reported on every rank before publishing output.
MPI process/transport failures and fatal resource failures are excluded.
"""
function SpectralOperators.deriv!(out, field, g::DistributedFourierGrid{T}, axis) where {T}
    comm = g.fft.comm
    j::Int, local_magnitude::Union{Float32,Float64} =
        HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "input validation") do
            axis isa Integer && 1 <= axis <= 3 ||
                throw(ArgumentError("distributed derivative axis must be an integer in 1:3"))
            _require_physical_field("input", field, g)
            _require_physical_field("output", out, g)
            magnitude = maximum(abs, parent(field))
            isfinite(magnitude) ||
                throw(ArgumentError("spectral operator input values must be finite"))
            Int(axis), magnitude
        end
    agreed_axis = MPI.Bcast(j, 0, comm)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "axis agreement") do
        j == agreed_axis ||
            throw(ArgumentError("distributed derivative axis differs between MPI ranks"))
    end
    magnitude = MPI.Allreduce(local_magnitude, MPI.MAX, comm)
    safe = SpectralOperators._fft_safe_magnitude(T, g.count, g.gain[j])
    scale = SpectralOperators._fft_input_scale(magnitude, safe)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "input scaling") do
        SpectralOperators._copy_fft_input!(parent(g.physical), parent(field), scale)
    end
    HybridPlasmaPIC.distributed_fft_forward!(g.spectral, g.fft, g.physical)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "mode multiplication") do
        modes = PencilArrays.global_view(g.spectral)
        SpectralOperators._apply_ik_store!(modes, modes, g.ik[j], j)
    end
    HybridPlasmaPIC.distributed_fft_inverse!(g.physical, g.fft, g.spectral)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "output preparation") do
        SpectralOperators._store_fft_output!(parent(g.staging), parent(g.physical), scale)
    end
    copyto!(parent(out), parent(g.staging))
    return out
end

function SpectralOperators.deriv(field, g::DistributedFourierGrid, axis)
    out = HybridPlasmaPIC.distributed_field(g)
    return SpectralOperators.deriv!(out, field, g, axis)
end

# Operators below follow deriv!'s contract: validate every operand in one
# collective stage, reduce the input maximum over ranks (and components), use
# the shared scale with the global transform length, agree each local step
# before the next transpose, and publish only after the final agreement.
# Public methods accept operands of any type for that validation, which
# establishes the grid's exact field type `R`; type assertions state that
# invariant before each typed numerical core, so rejected calls never compile one.

# Common FFT input scale from the communicator-wide maximum magnitude.
function _agreed_scale(
    g::DistributedFourierGrid{T},
    local_magnitude::Union{Float32,Float64},
    gain::Union{Float32,Float64},
) where {T}
    magnitude = MPI.Allreduce(local_magnitude, MPI.MAX, g.fft.comm)
    safe = SpectralOperators._fft_safe_magnitude(T, g.count, gain)
    return SpectralOperators._fft_input_scale(magnitude, safe)
end

# Scale one validated field into the physical transform input, agree any
# lost-nonzero rescaling error, then transform it into `spectral`.
function _forward_field!(spectral, field, g::DistributedFourierGrid, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "input scaling") do
        SpectralOperators._copy_fft_input!(parent(g.physical), parent(field), scale)
    end
    return HybridPlasmaPIC.distributed_fft_forward!(spectral, g.fft, g.physical)
end

# Normalized inverse transform, then stage the rescaled real part privately.
function _inverse_field!(staging, spectral, g::DistributedFourierGrid, scale, context)
    HybridPlasmaPIC.distributed_fft_inverse!(g.physical, g.fft, spectral)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "output preparation") do
        SpectralOperators._store_fft_output!(parent(staging), parent(g.physical), scale)
    end
    return staging
end

_staging_buffers(g::DistributedFourierGrid) = (g.staging, g.staging_work[1], g.staging_work[2])

function _publish_components!(out::NTuple{3,R}, staging::NTuple{3,R}) where {R}
    for c = 1:3
        copyto!(parent(out[c]), parent(staging[c]))
    end
    return out
end

"""
    gradient!(out::NTuple{3}, f, grid::DistributedFourierGrid)

Collectively store `∂_j f` in `out[j]` for a real rank-owned field: one forward
and three inverse transforms. Output components must not alias each other; one
may alias `f`. Operand, precision, alias and failure rules follow `deriv!`; all
three components are published together after the last collective success.
"""
function SpectralOperators.gradient!(
    out,
    f,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    local_magnitude::Union{Float32,Float64} = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        "Distributed gradient",
        "input validation",
    ) do
        _require_physical_field("input", f, g)
        _require_vector_field("output", out, g)
        _require_distinct_components(out, "gradient! output components must not alias each other")
        _finite_magnitude(f)
    end
    return _gradient!(out::NTuple{3,R}, f::R, g, local_magnitude::T)
end

function _gradient!(
    out::NTuple{3,R},
    f::R,
    g::DistributedFourierGrid{T,F,A,B,R},
    local_magnitude::T,
) where {T,F,A,B,R}
    context = "Distributed gradient"
    # The serial gradient bounds every component by the largest axis gain.
    scale = _agreed_scale(g, local_magnitude, max(g.gain[1], g.gain[2], g.gain[3]))
    _forward_field!(g.spectral, f, g, scale, context)
    modes = PencilArrays.global_view(g.spectral)
    product = PencilArrays.global_view(g.spectral_work[1])
    staging = _staging_buffers(g)
    for j = 1:3
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode multiplication") do
            SpectralOperators._apply_ik_store!(product, modes, g.ik[j], j)
        end
        _inverse_field!(staging[j], g.spectral_work[1], g, scale, context)
    end
    return _publish_components!(out, staging)
end

"""
    divergence!(out, v, grid::DistributedFourierGrid)

Collectively store `Σ_j ∂_j v[j]` for a tuple of at least three rank-owned fields;
as in the serial method only the first three components are used. `out` may
alias any component. Rules otherwise follow `deriv!`.
"""
function SpectralOperators.divergence!(
    out,
    v,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    # The serial divergence bounds the summed multiplier by the axis gains.
    gain = g.gain[1] + g.gain[2] + g.gain[3]
    local_magnitude::Union{Float32,Float64} = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        "Distributed divergence",
        "input validation",
    ) do
        v isa Tuple && length(v) >= 3 ||
            throw(ArgumentError("divergence! needs at least 3 vector components"))
        _require_physical_field("output", out, g)
        for c = 1:3
            _require_physical_field("input", v[c], g, c)
        end
        isfinite(gain) || throw(
            ArgumentError("divergence spectral multiplier bound is not representable as $T"),
        )
        _finite_vector_magnitude(v)
    end
    components = (v[1], v[2], v[3])
    return _divergence!(out::R, components::NTuple{3,R}, g, local_magnitude::T, gain)
end

function _divergence!(
    out::R,
    v::NTuple{3,R},
    g::DistributedFourierGrid{T,F,A,B,R},
    local_magnitude::T,
    gain::T,
) where {T,F,A,B,R}
    context = "Distributed divergence"
    scale = _agreed_scale(g, local_magnitude, gain)
    accumulator = PencilArrays.global_view(g.spectral)
    modes = PencilArrays.global_view(g.spectral_work[1])
    for j = 1:3
        _forward_field!(g.spectral_work[1], v[j], g, scale, context)
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode accumulation") do
            j == 1 && fill!(parent(g.spectral), zero(Complex{T}))
            SpectralOperators._apply_ik_accum!(accumulator, modes, g.ik[j], j, one(T))
        end
    end
    _inverse_field!(g.staging, g.spectral, g, scale, context)
    copyto!(parent(out), parent(g.staging))
    return out
end

"""
    curl!(out::NTuple{3}, A::NTuple{3}, grid::DistributedFourierGrid)

Collectively store the periodic curl of three rank-owned components. Output
components must be distinct from each other and from every input component.
Three forward and three inverse transforms reuse three spectral buffers; all
outputs are published together. Rules otherwise follow `deriv!`.
"""
function SpectralOperators.curl!(
    out,
    A,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    gain = g.gain[1] + g.gain[2] + g.gain[3]
    local_magnitude::Union{Float32,Float64} =
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, "Distributed curl", "input validation") do
            _require_vector_field("output", out, g)
            _require_vector_field("input", A, g)
            _require_distinct_components(out, "curl! output components must not alias each other")
            for c = 1:3, d = 1:3
                Base.mightalias(out[c], A[d]) &&
                    throw(ArgumentError("curl! output components must not alias input components"))
            end
            isfinite(gain) ||
                throw(ArgumentError("curl spectral multiplier bound is not representable as $T"))
            _finite_vector_magnitude(A)
        end
    return _curl!(out::NTuple{3,R}, A::NTuple{3,R}, g, local_magnitude::T, gain)
end

# dst = plus * ik[plus_axis] - minus * ik[minus_axis] on global logical modes.
function _curl_component!(dst, plus, plus_axis, minus, minus_axis, g, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode multiplication") do
        SpectralOperators._curl_component_store!(
            dst,
            plus,
            g.ik[plus_axis],
            plus_axis,
            minus,
            g.ik[minus_axis],
            minus_axis,
        )
    end
end

function _curl!(
    out::NTuple{3,R},
    A::NTuple{3,R},
    g::DistributedFourierGrid{T,F,A_,B,R},
    local_magnitude::T,
    gain::T,
) where {T,F,A_,B,R}
    context = "Distributed curl"
    scale = _agreed_scale(g, local_magnitude, gain)
    s0, s1, s2 = g.spectral, g.spectral_work[1], g.spectral_work[2]
    hat0, hat1, hat2 = map(PencilArrays.global_view, (s0, s1, s2))
    staging = _staging_buffers(g)
    _forward_field!(s1, A[2], g, scale, context)
    _forward_field!(s2, A[3], g, scale, context)
    # (curl A)_x = ∂_y A_z - ∂_z A_y
    _curl_component!(hat0, hat2, 2, hat1, 3, g, context)
    _inverse_field!(staging[1], s0, g, scale, context)
    _forward_field!(s0, A[1], g, scale, context)
    # (curl A)_y = ∂_z A_x - ∂_x A_z, stored over A_z's modes (no longer needed).
    _curl_component!(hat2, hat0, 3, hat2, 1, g, context)
    _inverse_field!(staging[2], s2, g, scale, context)
    # (curl A)_z = ∂_x A_y - ∂_y A_x, stored over A_y's modes.
    _curl_component!(hat1, hat1, 1, hat0, 2, g, context)
    _inverse_field!(staging[3], s1, g, scale, context)
    return _publish_components!(out, staging)
end

"""
    laplacian!(out, f, grid::DistributedFourierGrid)

Collectively store the periodic Laplacian `-|k|² f̂`, retaining Nyquist modes.
Throws `ArgumentError` on every rank before mutating `out` when the squared
wavenumber bound overflows. `out` may alias `f`; rules otherwise follow `deriv!`.
"""
function SpectralOperators.laplacian!(
    out,
    f,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    gain::Union{Float32,Float64}, local_magnitude::Union{Float32,Float64} =
        HybridPlasmaPIC._mpi_collective_stage(
            g.fft.comm,
            "Distributed Laplacian",
            "input validation",
        ) do
            _require_physical_field("input", f, g)
            _require_physical_field("output", out, g)
            SpectralOperators._laplacian_gain(g.kfull), _finite_magnitude(f)
        end
    return _laplacian!(out::R, f::R, g, local_magnitude::T, gain::T)
end

function _laplacian!(
    out::R,
    f::R,
    g::DistributedFourierGrid{T,F,A,B,R},
    local_magnitude::T,
    gain::T,
) where {T,F,A,B,R}
    context = "Distributed Laplacian"
    scale = _agreed_scale(g, local_magnitude, gain)
    _forward_field!(g.spectral, f, g, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode multiplication") do
        SpectralOperators._apply_laplacian_multiplier!(
            PencilArrays.global_view(g.spectral),
            g.kfull,
        )
    end
    _inverse_field!(g.staging, g.spectral, g, scale, context)
    copyto!(parent(out), parent(g.staging))
    return out
end

"""
    project_divfree!(B::NTuple{3}, grid::DistributedFourierGrid)

Collectively replace three distinct rank-owned components by their
divergence-free part `(I - k kᵀ/|k|²) B̂`, using the derivative-resolved
(Nyquist-zeroed) wavenumbers; modes whose such wavenumber is zero, including DC
and mixed DC/Nyquist modes, are unchanged. Rules otherwise follow `deriv!`.
"""
function SpectralOperators.project_divfree!(
    B,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    local_magnitude::Union{Float32,Float64} = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        "Distributed divergence-free projection",
        "input validation",
    ) do
        _require_vector_field("component", B, g)
        _require_distinct_components(B, "project_divfree! components must not alias each other")
        _finite_vector_magnitude(B)
    end
    return _project_divfree!(B::NTuple{3,R}, g, local_magnitude::T)
end

function _project_divfree!(
    B::NTuple{3,R},
    g::DistributedFourierGrid{T,F,A,B_,R},
    local_magnitude::T,
) where {T,F,A,B_,R}
    context = "Distributed divergence-free projection"
    # The projector is dimensionless; the serial method uses a unit gain.
    scale = _agreed_scale(g, local_magnitude, one(T))
    spectral = (g.spectral, g.spectral_work[1], g.spectral_work[2])
    for c = 1:3
        _forward_field!(spectral[c], B[c], g, scale, context)
    end
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode projection") do
        SpectralOperators._project_fourier_modes!(map(PencilArrays.global_view, spectral), g.kvec)
    end
    staging = _staging_buffers(g)
    for c = 1:3
        _inverse_field!(staging[c], spectral[c], g, scale, context)
    end
    return _publish_components!(B, staging)
end

"""
    exp_filter!(field, grid::DistributedFourierGrid; α=36, p=8)

Collectively apply the serial mode-ratio filter
`exp(-α max_d(|m_d|/(n_d÷2))^(2p))` in place; DC is retained exactly. `α` must
be real and nonnegative (converted to the grid precision) and `p` an integer in
`1:typemax(Int)`; both must agree on every rank. Rules otherwise follow `deriv!`.
"""
function SpectralOperators.exp_filter!(
    field,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R};
    α = 36,
    p = 8,
) where {T,R}
    alpha::Union{Float32,Float64}, power::Int, local_magnitude::Union{Float32,Float64} =
        HybridPlasmaPIC._mpi_collective_stage(
            g.fft.comm,
            "Distributed exponential filter",
            "input validation",
        ) do
            _require_physical_field("field", field, g)
            p isa Integer || throw(ArgumentError("p must be an integer"))
            p >= 1 || throw(ArgumentError("p must be ≥ 1"))
            p <= typemax(Int) || throw(ArgumentError("p must not exceed typemax(Int)"))
            α isa Real || throw(ArgumentError("α must be a real number"))
            α >= 0 || throw(ArgumentError("α must be ≥ 0"))
            T(α), Int(p), _finite_magnitude(field)
        end
    return _exp_filter!(field::R, g, local_magnitude::T, alpha::T, power)
end

function _exp_filter!(
    field::R,
    g::DistributedFourierGrid{T,F,A,B,R},
    local_magnitude::T,
    alpha::T,
    power::Int,
) where {T,F,A,B,R}
    comm = g.fft.comm
    context = "Distributed exponential filter"
    # Float32 values widen exactly, so one Float64 broadcast covers both precisions.
    agreed_alpha = MPI.Bcast(Float64(alpha), 0, comm)
    agreed_power = MPI.Bcast(power, 0, comm)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "parameter agreement") do
        isequal(Float64(alpha), agreed_alpha) && power == agreed_power || throw(
            ArgumentError("distributed exponential filter parameters differ between MPI ranks"),
        )
    end
    scale = _agreed_scale(g, local_magnitude, one(T))
    _forward_field!(g.spectral, field, g, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "mode multiplication") do
        SpectralOperators._apply_exp_filter_multiplier!(
            PencilArrays.global_view(g.spectral),
            g.midx,
            g.n,
            alpha,
            power,
        )
    end
    _inverse_field!(g.staging, g.spectral, g, scale, context)
    copyto!(parent(field), parent(g.staging))
    return field
end

"""
    dealias_two_thirds!(field, grid::DistributedFourierGrid)

Collectively retain exactly the modes with `3abs(m_d) < n_d` on every axis, in
place; DC is retained. Rules otherwise follow `deriv!`.
"""
function SpectralOperators.dealias_two_thirds!(
    field,
    g::DistributedFourierGrid{T,<:Any,<:Any,<:Any,R},
) where {T,R}
    local_magnitude::Union{Float32,Float64} = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        "Distributed two-thirds dealiasing",
        "input validation",
    ) do
        _require_physical_field("field", field, g)
        _finite_magnitude(field)
    end
    return _dealias_two_thirds!(field::R, g, local_magnitude::T)
end

function _dealias_two_thirds!(
    field::R,
    g::DistributedFourierGrid{T,F,A,B,R},
    local_magnitude::T,
) where {T,F,A,B,R}
    context = "Distributed two-thirds dealiasing"
    scale = _agreed_scale(g, local_magnitude, one(T))
    _forward_field!(g.spectral, field, g, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "mode multiplication") do
        SpectralOperators._apply_two_thirds_mask!(PencilArrays.global_view(g.spectral), g.midx, g.n)
    end
    _inverse_field!(g.staging, g.spectral, g, scale, context)
    copyto!(parent(field), parent(g.staging))
    return field
end
