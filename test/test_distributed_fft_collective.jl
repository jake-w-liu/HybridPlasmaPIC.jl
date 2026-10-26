module DistributedFFTCollectiveTests
using HybridPlasmaPIC, Test
import MPI, PencilArrays, PencilFFTs, FFTW, Serialization
const H = HybridPlasmaPIC
const F = PencilFFTs.Transforms
H.ensure_mpi_initialized!()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const SIZE = MPI.Comm_size(COMM)

# Fault injection at the dependency boundary, with the normal FFT's metadata.
# These types are used only to verify collective failure propagation.
struct SerializationFailure <: F.AbstractTransform end
struct PlanningFailure <: F.AbstractTransform end
const InjectedTransform = Union{SerializationFailure,PlanningFailure}
F.is_inplace(::InjectedTransform) = false
F.length_output(::InjectedTransform, n::Integer) = n
F.eltype_input(::InjectedTransform, ::Type{T}) where {T} = Complex{T}
F.eltype_output(::InjectedTransform, ::Type{Complex{T}}) where {T} = Complex{T}
F.binv(::InjectedTransform, d::Integer) = F.BFFT()
function Serialization.serialize(::Serialization.AbstractSerializer, ::SerializationFailure)
    throw(ArgumentError("injected transform serialization failure"))
end
function F.plan(::PlanningFailure, array, dims; kwargs...)
    RANK == SIZE - 1 && throw(ArgumentError("injected transform planning failure"))
    return F.plan(F.FFT(), array, dims; kwargs...)
end

function collective_failure(f, local_failure)
    error = try
        f()
        nothing
    catch exception
        exception
    end
    @test error isa Exception
    @test MPI.Allreduce(Int(error !== nothing), +, COMM) == SIZE
    @test local_failure ? error isa ArgumentError : error isa ErrorException
    @test MPI.Allreduce(1, +, COMM) == SIZE
    return error
end

@testset "Distributed FFT collective validation and agreement" begin
    defaults =
        (n = (8, 6, 4), transform = F.FFT(), periodic = (true, true, true), real_type = Float64)
    invalid = (
        (n = (8, 6),),
        (n = [8, 6, 4],),
        (n = (8.0, 6, 4),),
        (n = (8, 0, 4),),
        (n = (8, -1, 4),),
        (n = (BigInt(typemax(Int)) + 1, 6, 4),),
        (periodic = (true, false, true),),
        (periodic = (true, true),),
        (periodic = (1, 1, 1),),
        (transform = :fft,),
        (transform = (F.FFT(), F.FFT()),),
        (transform = (F.FFT(), F.RFFT(), F.FFT()),),
        (transform = (F.FFT!(), F.FFT(), F.FFT()),),
        (transform = SerializationFailure(),),
        (real_type = Float16,),
        (real_type = BigFloat,),
        (real_type = Int,),
        (real_type = ComplexF64,),
        (real_type = 32,),
        (real_type = nothing,),
        (real_type = "Float32",),
    )
    for origin = 0:SIZE-1, bad in invalid
        settings = RANK == origin ? merge(defaults, bad) : defaults
        collective_failure(RANK == origin) do
            distributed_fft_plan(
                settings.n;
                comm = COMM,
                transform = settings.transform,
                periodic = settings.periodic,
                real_type = settings.real_type,
            )
        end
    end
    if SIZE > 1
        for scenario in (:size, :direction, :real_parity, :precision)
            n = scenario === :size && RANK == SIZE - 1 ? (10, 6, 4) : defaults.n
            transform = if scenario === :real_parity
                F.BRFFT(RANK == SIZE - 1 ? 7 : 6)
            elseif scenario === :direction && RANK == SIZE - 1
                F.BFFT()
            else
                F.FFT()
            end
            collective_failure(RANK == SIZE - 1) do
                real_type = scenario === :precision && RANK == SIZE - 1 ? Float32 : Float64
                distributed_fft_plan(n; comm = COMM, transform, real_type)
            end
        end
        collective_failure(true) do
            distributed_fft_plan((1, 1, 1); comm = COMM)
        end
    end
    collective_failure(RANK == SIZE - 1) do
        distributed_fft_plan(defaults.n; comm = COMM, transform = PlanningFailure())
    end
end

function transform_reference(kind, array, n)
    kind in (:fft, :fft_inplace) && return FFTW.fft(array)
    kind === :bfft && return FFTW.bfft(array)
    kind in (:rfft, :rfft_inplace) && return FFTW.rfft(array)
    kind in (:r2r, :r2r_inplace) && return FFTW.r2r(array, FFTW.REDFT10)
    kind === :identity && return copy(array)
    kind === :mixed && return FFTW.rfft(array, (2, 3))
    kind === :brfft_even && return FFTW.brfft(array, 2n[3] - 2, (3, 2, 1))
    kind === :brfft_odd && return FFTW.brfft(array, 2n[3] - 1, (3, 2, 1))
    error("unknown transform fixture")
end

@testset "Distributed FFT supported transforms retain serial mathematics" begin
    configurations = (
        (:fft, F.FFT()),
        (:bfft, F.BFFT()),
        (:rfft, F.RFFT()),
        (:r2r, F.R2R(FFTW.REDFT10)),
        (:identity, F.NoTransform()),
        (:mixed, (F.NoTransform(), F.RFFT(), F.FFT())),
        (:brfft_even, F.BRFFT(6)),
        (:brfft_odd, F.BRFFT(7)),
        (:fft_inplace, F.FFT!()),
        (:rfft_inplace, F.RFFT!()),
        (:r2r_inplace, F.R2R!(FFTW.REDFT10)),
    )
    for T in (Float32, Float64), n in ((8, 6, 4), (9, 7, 5)), (kind, requested) in configurations
        # Single-transform and expanded forms must describe the same plan.
        transforms = requested isa Tuple ? requested : F.expand_dims(requested, Val(3))
        argument = isodd(RANK) ? transforms : requested
        shape = isodd(RANK) ? (BigInt(n[1]), Int32(n[2]), UInt8(n[3])) : n
        # Explicit/default Float64 requests must also agree across ranks.
        plan = if T === Float64 && isodd(RANK)
            distributed_fft_plan(shape; comm = COMM, transform = argument)
        else
            distributed_fft_plan(shape; comm = COMM, transform = argument, real_type = T)
        end
        @test plan.n == n
        inplace = F.is_inplace(transforms...)
        input = distributed_fft_input(plan)
        output = inplace ? input : distributed_fft_output(plan)
        input_view = inplace ? first(input) : input
        output_view = inplace ? last(output) : output
        @test real(eltype(input_view)) === T
        @test real(eltype(output_view)) === T
        global_input = PencilArrays.global_view(input_view)
        for I in CartesianIndices(global_input)
            value = if kind in (:brfft_even, :brfft_odd)
                # A real constant spectrum is Hermitian for either final parity.
                1.0
            else
                cospi(2(I[1] - 1) / n[1]) + 2sinpi(2(I[2] - 1) / n[2]) + 0.3cospi(4(I[3] - 1) / n[3])
            end
            global_input[I] =
                eltype(global_input) <: Complex && !(kind in (:brfft_even, :brfft_odd)) ?
                complex(value, 0.2sinpi(2(I[3] - 1) / n[3])) : value
        end
        original_local = copy(parent(input_view))
        original = PencilArrays.gather(input_view)
        @test MPI.Allreduce(length(parent(input_view)), +, COMM) == prod(n)
        @test SIZE == 1 || length(parent(input_view)) < prod(n)
        @test distributed_fft_forward!(output, plan, input) === output
        transformed = PencilArrays.gather(output_view)
        if RANK == 0
            expected = transform_reference(kind, original, n)
            scale = max(1.0, maximum(abs, original))
            @test size(transformed) == size(expected)
            @test isapprox(
                transformed,
                expected;
                rtol = 256eps(T),
                atol = 256eps(T) * length(original) * scale,
            )
        end
        @test distributed_fft_inverse!(input, plan, output) === input
        @test isapprox(
            parent(input_view),
            original_local;
            rtol = 256eps(T),
            atol = 256eps(T) * max(1.0, maximum(abs, original_local)),
        )
        @test MPI.Allreduce(1, +, COMM) == SIZE
    end
end
end
