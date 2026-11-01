module TestYeeParticleCharge
using HybridPlasmaPIC, Test
import JET
include("yee_charge_oracle.jl")
const O = YeeChargeOracle

function population(::Type{T}, x; q = one(T), weight = one(T)) where {T}
    ps = ParticleSet{1,T}(length(x); q)
    ps.x[1] .= x
    ps.weight .= weight
    return ps
end

function abstract_storage(ps::ParticleSet{1,T}) where {T}
    P = ParticleSet{
        1,T,AbstractVector{T},AbstractVector{T},AbstractVector{T},
        AbstractVector{UInt64},AbstractVector{UInt32},
    }
    return P(ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
end

state(ps) = (
    map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m,
)

function check_values(actual, exact, magnitude, count)
    T = eltype(actual)
    for i in eachindex(actual)
        @test isfinite(actual[i])
        @test abs(O.rat(actual[i]) - exact[i]) <=
              O.error_bound(T, exact[i], magnitude[i], count)
    end
end

@noinline function in_mode(operation::F, hardware, flush) where {F}
    precompile(operation, ())
    saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
    requested, active = Base.Rounding.to_fenv(hardware), nothing
    try
        Base.Rounding.setrounding_raw(Float64, requested)
        Base.set_zero_subnormals(flush)
        active = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
        return operation()
    finally
        Base.Rounding.setrounding_raw(Float64, saved[1])
        Base.set_zero_subnormals(saved[2])
        @test active == (requested, flush)
    end
end

@testset "Closed CIC charge and path impulse" begin
    for T in (Float32, Float64), length in (one(T), T(0.1))
        g = YeeGrid((8,), (length,))
        dense = population(T, length .* T[0, 0.125, 0.375, 0.625, 1])
        dense.weight .= T[1, 0, 0.5, 2, 1]
        second = population(T, reverse(dense.x[1]); q = T(-0.5))
        viewed = ParticleSet{1,T}(
            (view(second.x[1], :),),
            second.v,
            second.weight,
            second.id,
            second.tag,
            second.q,
            second.m,
        )
        # Exercise abstract component declarations as well as broad collections.
        broad = abstract_storage(viewed)
        broad.m = T(NaN)
        broad.v[1] .= T(NaN) # Unused values do not enter charge/current physics.
        species = (dense, broad)
        old = AbstractVector{T}[copy(ps.x[1]) for ps in species]
        new = AbstractVector{T}[view(reverse(copy(x)), :) for x in old]
        new_before = map(copy, new)
        expected_charge = O.charge(species, old, 8, length)
        expected_impulse = O.impulse(species, old, new, 8, length)
        original = map(state, species)
        for collection in (species, ParticleSet{1,T}[species...], Any[species...])
            work = YeeChargeWorkspace(g)
            charge, impulse = fill(T(7), 9), fill(T(8), 8)
            @test yee_particle_charge!(charge, collection, g, work) === charge
            check_values(charge, expected_charge..., 10)
            @test yee_longitudinal_impulse!(impulse, collection, old, new, g, work) === impulse
            check_values(impulse, expected_impulse..., 10)
            reversed = similar(impulse)
            yee_longitudinal_impulse!(reversed, collection, new, old, g, work)
            @test reversed == -impulse
            @test isequal(map(state, species), original)
            @test old == [collect(ps.x[1]) for ps in species]
            @test new == new_before
        end
        work = YeeChargeWorkspace(g)
        charge, impulse = ones(T, 9), ones(T, 8)
        yee_particle_charge!(charge, ParticleSet{1,T}[], g, work)
        yee_longitudinal_impulse!(impulse, Any[], (), (), g, work)
        @test all(iszero, charge) && all(iszero, impulse)
    end
end

@testset "Closed-wall continuity, exact range and publication" begin
    for T in (Float32, Float64)
        g, ps = YeeGrid((8,), (one(T),)), population(T, T[0])
        work = YeeChargeWorkspace(g)
        q0, q1, impulse = zeros(T, 9), zeros(T, 9), zeros(T, 8)
        yee_particle_charge!(q0, ps, g, work)
        yee_longitudinal_impulse!(impulse, ps, T[0], T[1], g, work)
        ps.x[1][1] = one(T)
        yee_particle_charge!(q1, ps, g, work)
        @test O.rat.(q1) - O.rat.(q0) + O.incidence(impulse) == zeros(O.Q, 9)
        yee_longitudinal_impulse!(impulse, ps, T[0], T[nextfloat(zero(T))], g, work)
        @test impulse[1] == 8nextfloat(zero(T))
        @test all(iszero, impulse[2:end])

        large = population(T, T[0]; q = floatmax(T), weight = floatmax(T))
        negative = population(T, T[0]; q = -floatmax(T), weight = floatmax(T))
        unit = population(T, T[0])
        yee_particle_charge!(q0, Any[large, negative, unit], g, work)
        @test q0 == T[1, 0, 0, 0, 0, 0, 0, 0, 0]
        q0 .= T(7)
        @test_throws ArgumentError yee_particle_charge!(q0, large, g, work)
        @test q0 == fill(T(7), 9)
        yee_particle_charge!(q0, unit, g, work) # Scratch is reusable after rejection.
        @test q0[1] == one(T)

        bad = population(T, T[-1])
        bad.q = zero(T) # Neutral particles retain the same geometric admission.
        q0 .= T(7)
        work.sum .= 11.0
        work.candidate .= T(13)
        before = (copy(q0), copy(work.sum), copy(work.candidate))
        @test_throws ArgumentError yee_particle_charge!(q0, Any[unit, bad], g, work)
        @test (q0, work.sum, work.candidate) == before
        original_impulse = copy(impulse)
        @test_throws ArgumentError yee_longitudinal_impulse!(
            impulse, Any[unit, unit], Any[T[0], T[0]], Any[T[1], T[NaN]], g, work,
        )
        @test impulse == original_impulse
        @test (q0, work.sum, work.candidate) == before
        @test_throws ArgumentError yee_particle_charge!(work.candidate, unit, g, work)
        @test (q0, work.sum, work.candidate) == before
    end
end

@testset "Stored operands in nondefault IEEE contexts" begin
    saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
    for T in (Float32, Float64), (hardware, flush) in ((RoundUp, false), (RoundNearest, true))
        g, sub = YeeGrid((8,), (one(T),)), nextfloat(zero(T))
        work = YeeChargeWorkspace(g)
        out, impulse = zeros(T, 9), zeros(T, 8)
        ps = abstract_storage(population(T, T[0]; q = sub))
        species = Any[ps]
        exact = O.rat(sub)
        in_mode(() -> yee_particle_charge!(out, species, g, work), hardware, flush)
        @test O.rounded(out[1], exact) && all(iszero, out[2:end])
        ps.q = one(T)
        for (old, new, expected) in ((T[0], T[sub], 8exact), (T[sub], T[0], -8exact))
            in_mode(
                () -> yee_longitudinal_impulse!(impulse, ps, old, new, g, work),
                hardware,
                flush,
            )
            @test O.rounded(impulse[1], expected) && all(iszero, impulse[2:end])
        end
        ps.q = -zero(T)
        in_mode(() -> yee_particle_charge!(out, ps, g, work), hardware, flush)
        @test all(iszero, out)
        ps.q = one(T)
        for bad_weight in (false, true)
            ps.x[1][1] = bad_weight ? zero(T) : -sub
            ps.weight[1] = bad_weight ? -sub : one(T)
            out .= T(7)
            work.sum .= 11.0
            work.candidate .= T(13)
            before = (copy(out), copy(work.sum), copy(work.candidate))
            @test_throws ArgumentError in_mode(
                () -> yee_particle_charge!(out, ps, g, work),
                hardware,
                flush,
            )
            @test (out, work.sum, work.candidate) == before
        end
        # This admitted grid has a subnormal stored L, not a rounded zero.
        small = YeeGrid((2,), (prevfloat(floatmin(T)),))
        wall = population(T, T[small.L[1]])
        wall_charge, wall_work = zeros(T, 3), YeeChargeWorkspace(small)
        in_mode(() -> yee_particle_charge!(wall_charge, wall, small, wall_work), hardware, flush)
        @test wall_charge == T[0, 0, 1]
        wrong = YeeChargeWorkspace(YeeGrid((2,), (prevfloat(small.L[1]),)))
        wrong.sum .= 11.0
        wrong.candidate .= T(13)
        before = (copy(wall_charge), copy(wrong.sum), copy(wrong.candidate))
        @test_throws ArgumentError in_mode(
            () -> yee_particle_charge!(wall_charge, wall, small, wrong),
            hardware,
            flush,
        )
        @test (wall_charge, wrong.sum, wrong.candidate) == before
    end
    @test (Base.Rounding.rounding_raw(Float64), get_zero_subnormals()) == saved
end

@testset "Yee scalar domain through broad public collections" begin
    for T in (Float32, Float64)
        G, W = YeeGrid{1,T}, YeeChargeWorkspace{T}
        for P in (Vector{ParticleSet{1,T}}, Vector{Any})
            for (operation, types) in (
                (yee_particle_charge!, Tuple{Vector{T},P,G,W}),
                (
                    yee_longitudinal_impulse!,
                    Tuple{Vector{T},P,Vector{AbstractVector{T}},Vector{AbstractVector{T}},G,W},
                ),
            )
                report = JET.report_call(
                    operation, types; sourceinfo = :full, stacktrace_types_limit = 0,
                )
                if !isempty(JET.get_reports(report))
                    show(IOContext(stdout, :limit => false), MIME"text/plain"(), report)
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end
end
