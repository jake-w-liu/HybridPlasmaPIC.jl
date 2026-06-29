#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

struct ValidationOffsetVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    first_index::Int
end

Base.size(v::ValidationOffsetVector) = size(v.data)
Base.axes(v::ValidationOffsetVector) = (v.first_index:(v.first_index+length(v.data)-1),)
Base.IndexStyle(::Type{<:ValidationOffsetVector}) = IndexLinear()
Base.getindex(v::ValidationOffsetVector, i::Int) = v.data[i-v.first_index+1]

function _throws_expected(f)
    try
        f()
        return false
    catch
        return true
    end
end

_shape_width(::NGP) = 1
_shape_width(::CIC) = 2
_shape_width(::TSC) = 3

function _particle_roundtrip_error(a::ParticleSet, b::ParticleSet)
    err = Float64(nparticles(a) == nparticles(b) && a.q == b.q && a.m == b.m ? 0.0 : 1.0)
    for d in eachindex(a.x)
        err = max(err, maximum(abs, collect(a.x[d]) .- collect(b.x[d])))
    end
    for c = 1:3
        err = max(err, maximum(abs, collect(a.v[c]) .- collect(b.v[c])))
    end
    err = max(err, maximum(abs, collect(a.weight) .- collect(b.weight)))
    err = max(err, a.id == b.id ? 0.0 : 1.0)
    err = max(err, a.tag == b.tag ? 0.0 : 1.0)
    return err
end

function case_13_threaded_backend_api_validation(artifact_dir::AbstractString)
    id = "13_threaded_backend_api_validation"
    max_deposit_error = 0.0
    max_density_error = 0.0
    max_conservation_error = 0.0
    offset_error = 0.0

    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(d -> 7 + 2d, D)
        g = FourierGrid(n, ntuple(d -> 1.25 + 0.25d, D))
        np = 307 + 17D
        ps = ParticleSet{D,Float64}(np)
        load_uniform!(
            ps,
            MersenneTwister(1000 + 17D + _shape_width(shape)),
            ntuple(_ -> 0.0, D),
            g.L,
        )
        ps.weight .= 0.5 .+ rand(MersenneTwister(2000 + 19D + _shape_width(shape)), np)
        vals = randn(MersenneTwister(3000 + 23D + _shape_width(shape)), np)

        serial = zeros(Float64, n)
        threaded = similar(serial)
        deposit_scalar!(serial, ps, vals, g, shape)
        deposit_scalar_threaded!(threaded, ps, vals, g, shape)
        max_deposit_error = max(max_deposit_error, maximum(abs, threaded .- serial))
        max_conservation_error = max(max_conservation_error, abs(sum(threaded) - sum(vals)))

        nserial = zeros(Float64, n)
        nthreaded = similar(nserial)
        density!(nserial, ps, g, shape)
        density_threaded!(nthreaded, ps, g, shape)
        max_density_error = max(max_density_error, maximum(abs, nthreaded .- nserial))

        if D == 2 && shape isa CIC
            offset_vals = ValidationOffsetVector(vals, -7)
            offset_threaded = similar(serial)
            deposit_scalar_threaded!(offset_threaded, ps, offset_vals, g, shape)
            offset_error = maximum(abs, offset_threaded .- serial)
        end
    end

    g0 = FourierGrid((8, 8), (1.0, 1.0))
    ps0 = ParticleSet{2,Float64}(0)
    z = ones(Float64, g0.n)
    deposit_scalar_threaded!(z, ps0, Float64[], g0, TSC())
    zero_particle_error = maximum(abs, z)

    gs = FourierGrid((32,), (2π,))
    xs = [(i - 1) * gs.dx[1] for i = 1:gs.n[1]]
    f_workspace = cos.(3 .* xs)
    f_default = copy(f_workspace)
    work = BinomialSmoothWorkspace(gs)
    binomial_smooth!(f_workspace, gs, work; passes = 2)
    binomial_smooth!(f_default, gs; passes = 2)
    workspace_error = maximum(abs, f_workspace .- f_default)

    ps = ParticleSet{2,Float64}(4; q = 2.0, m = 3.0)
    ps.x[1] .= [0.1, 0.2, 0.3, 0.4]
    ps.x[2] .= [1.1, 1.2, 1.3, 1.4]
    ps.v[1] .= [2.1, 2.2, 2.3, 2.4]
    ps.v[2] .= [3.1, 3.2, 3.3, 3.4]
    ps.v[3] .= [4.1, 4.2, 4.3, 4.4]
    ps.weight .= [0.5, 0.6, 0.7, 0.8]
    ps.id .= UInt64[11, 12, 13, 14]
    ps.tag .= UInt32[21, 22, 23, 24]
    host = copy_particles_to_backend(Val(:cpu), ps)
    roundtrip = copy_particles_to_host(host)
    backend_error =
        particle_storage_backend(ps) == :cpu &&
        all(
            particle_array_backend(A) == :cpu for A in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
        ) &&
        particle_storage_backend(host) == :cpu ? 0.0 : 1.0
    host_copy_error =
        max(_particle_roundtrip_error(host, ps), _particle_roundtrip_error(roundtrip, host))

    mem = BackendMemoryStatus(:cpu, true, false; total_bytes = 128, free_bytes = 32)
    memory_pressure_error = abs(memory_pressure(mem) - 0.75)
    cpu_status = backend_memory_status(Val(:cpu))
    cpu_status_error =
        cpu_status.backend == :cpu &&
        cpu_status.device_available &&
        !cpu_status.pool_supported &&
        backend_memory_status(ps).backend == :cpu &&
        reclaim_backend_memory!(Val(:cpu)) == false ? 0.0 : 1.0

    extension_error =
        supported_extensions() == (:cuda, :metal, :io, :pencilfft) &&
        extension_name(Val(:cuda)) == :HybridPlasmaPICCUDAExt &&
        extension_dependency_name(Val(:io)) == :HDF5 &&
        extension_loaded(Val(:cuda)) isa Bool &&
        loaded_extensions() isa Tuple ? 0.0 : 1.0
    missing_extension_error =
        extension_loaded(Val(:pencilfft)) ? 0.0 :
        (
            _throws_expected(() -> require_extension(Val(:pencilfft))) &&
            _throws_expected(() -> extension_dependency_module(Val(:pencilfft))) &&
            _throws_expected(() -> extension_device_array_type(Val(:pencilfft))) &&
            _throws_expected(() -> disallow_scalar_indexing!(Val(:pencilfft))) &&
            _throws_expected(() -> 17_distributed_fft_roundtrip_error()) ? 0.0 : 1.0
        )
    prepare_cpu_error = prepare_gpu_backend!(Val(:cpu)) === nothing ? 0.0 : 1.0
    abstract_contract_error =
        CIC() isa ShapeFunction && IsothermalElectrons(1.0) isa ElectronClosure ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "13_threaded_backend_api_validation.csv")
    rows = (
        (
            "threaded_deposit_max_abs_error",
            max_deposit_error,
            0.0,
            "absolute",
            max_deposit_error,
            1e-10,
        ),
        (
            "threaded_density_max_abs_error",
            max_density_error,
            0.0,
            "absolute",
            max_density_error,
            1e-10,
        ),
        (
            "threaded_deposit_conservation_error",
            max_conservation_error,
            0.0,
            "absolute",
            max_conservation_error,
            1e-10,
        ),
        (
            "threaded_offset_values_max_abs_error",
            offset_error,
            0.0,
            "absolute",
            offset_error,
            1e-12,
        ),
        (
            "threaded_zero_particle_error",
            zero_particle_error,
            0.0,
            "absolute",
            zero_particle_error,
            0.0,
        ),
        (
            "binomial_workspace_equivalence_error",
            workspace_error,
            0.0,
            "absolute",
            workspace_error,
            0.0,
        ),
        ("cpu_backend_contract_error", backend_error, 0.0, "absolute", backend_error, 0.0),
        (
            "cpu_backend_copy_roundtrip_error",
            host_copy_error,
            0.0,
            "absolute",
            host_copy_error,
            0.0,
        ),
        (
            "backend_memory_pressure_error",
            memory_pressure(mem),
            0.75,
            "absolute",
            memory_pressure_error,
            0.0,
        ),
        (
            "cpu_backend_status_contract_error",
            cpu_status_error,
            0.0,
            "absolute",
            cpu_status_error,
            0.0,
        ),
        (
            "extension_registry_contract_error",
            extension_error,
            0.0,
            "absolute",
            extension_error,
            0.0,
        ),
        (
            "missing_extension_fallback_contract_error",
            missing_extension_error,
            0.0,
            "absolute",
            missing_extension_error,
            0.0,
        ),
        (
            "prepare_cpu_backend_contract_error",
            prepare_cpu_error,
            0.0,
            "absolute",
            prepare_cpu_error,
            0.0,
        ),
        (
            "abstract_api_type_contract_error",
            abstract_contract_error,
            0.0,
            "absolute",
            abstract_contract_error,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "threaded_backend_api",
        reference_kind = "analytic_or_contract",
        reference = "threaded-vs-serial deposition, CPU backend roundtrip, extension fallback contracts",
        rows = rows,
        artifact = artifact,
    )
end

VALIDATION_CASE = ValidationCase(
    id = "13_threaded_backend_api_validation",
    default = true,
    description = "Threaded deposition/density, CPU backend storage, memory telemetry, and extension API contracts.",
    runner = case_13_threaded_backend_api_validation,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
