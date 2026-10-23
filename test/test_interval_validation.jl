module IntervalValidationTests
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

function population(::Val{D}, ::Type{T}) where {D,T}
    ps = ParticleSet{D,T}(2^D)
    for d = 1:D
        fill!(ps.x[d], T(0.25))
    end
    for c = 1:3
        fill!(ps.v[c], T(c))
    end
    ps.weight .= T(2)
    ps.id .= UInt64.(1:nparticles(ps))
    ps.tag .= UInt32.(11:(10+nparticles(ps)))
    return ps
end
state(ps) = (map(copy, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)), ps.q, ps.m)
function caught(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end
function operation!(op, ps::ParticleSet{D}, lo, hi, rng) where {D}
    op === :uniform && return load_uniform!(ps, rng, lo, hi)
    op === :lattice && return load_lattice!(ps, lo, hi, ntuple(_ -> 2, D))
    op === :periodic && return apply_periodic!(ps, lo, hi)
    op === :reflecting && return apply_reflecting!(ps, lo, hi)
    op === :absorbing && return apply_absorbing!(ps, lo, hi)
    error("unknown test operation")
end

@testset "Interval errors preserve indexed messages, order and state" begin
    for T in (Float32, Float64), D = 1:3, axis = 1:D
        good_lo, good_hi = ntuple(_ -> zero(T), D), ntuple(_ -> one(T), D)
        invalid = Tuple[]
        for bad in (T(NaN), T(Inf), T(-Inf))
            push!(invalid, (Base.setindex(good_lo, bad, axis), good_hi, "lo[$axis] must be finite"))
            push!(invalid, (good_lo, Base.setindex(good_hi, bad, axis), "hi[$axis] must be finite"))
        end
        push!(
            invalid,
            (
                good_lo,
                Base.setindex(good_hi, zero(T), axis),
                "hi[$axis] must be greater than lo[$axis]",
            ),
        )
        push!(
            invalid,
            (
                Base.setindex(good_lo, T(Inf), axis),
                Base.setindex(good_hi, T(NaN), 1),
                "lo[$axis] must be finite",
            ),
        )
        for (lo, hi, message) in invalid,
            op in (:uniform, :lattice, :periodic, :reflecting, :absorbing)

            ps = population(Val(D), T)
            before = state(ps)
            refs = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            rng = MersenneTwister(0x5566)
            reference_rng = copy(rng)
            err = caught(() -> operation!(op, ps, lo, hi, rng))
            @test err isa ArgumentError
            if err isa ArgumentError
                @test err.msg == message
            end
            @test isequal(state(ps), before)
            @test all(
                (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)[i] === refs[i] for i in eachindex(refs)
            )
            @test rand(rng, UInt64, 8) == rand(reference_rng, UInt64, 8)
        end
    end
    # Finite source values may overflow when converted to the destination type.
    err = caught(() -> H._validated_open_interval((0.0,), (floatmax(Float64),), Float32))
    @test err isa ArgumentError
    @test err.msg == "hi[1] must be finite"
    @test H._require_finite_real("control", 2, Float32) === 2.0f0
    err = caught(() -> H._require_finite_real("control", Inf, Float32))
    @test err isa ArgumentError
    @test err.msg == "control must be finite"
end

function allocation_samples(f::F, args::A) where {F,A<:Tuple}
    for _ = 1:3
        f(args...)
    end
    return [@allocated(f(args...)) for _ = 1:3]
end
function interval_allocation_samples(lo, hi, ::Type{T}) where {T}
    for _ = 1:3
        H._validated_open_interval(lo, hi, T)
    end
    return [@allocated(H._validated_open_interval(lo, hi, T)) for _ = 1:3]
end
@testset "Finite interval and public boundary allocation" begin
    for T in (Float32, Float64), D = 1:3
        lo, hi = ntuple(_ -> zero(T), D), ntuple(_ -> one(T), D)
        @test H._validated_open_interval(lo, hi, T) === (lo, hi)
        @test interval_allocation_samples(lo, hi, T) == zeros(Int, 3)
        ps = population(Val(D), T)
        rng = MersenneTwister(0x7788)
        for op in (apply_periodic!, apply_reflecting!, apply_absorbing!)
            @test allocation_samples(op, (ps, lo, hi)) == zeros(Int, 3)
        end
        @test allocation_samples(load_uniform!, (ps, rng, lo, hi)) == zeros(Int, 3)
        @test allocation_samples(load_lattice!, (ps, lo, hi, ntuple(_ -> 2, D))) == zeros(Int, 3)
    end
    setprecision(BigFloat, 113) do
        lo, hi = (BigFloat(-1),), (BigFloat(2),)
        actual = H._validated_open_interval(lo, hi, BigFloat)
        @test actual == (lo, hi)
        @test all(precision(x) == 113 for interval in actual for x in interval)
        err = caught(() -> H._validated_open_interval(lo, (BigFloat(Inf),), BigFloat))
        @test err isa ArgumentError
        @test err.msg == "hi[1] must be finite"
    end
end
end
