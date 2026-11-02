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
        1,
        T,
        AbstractVector{T},
        AbstractVector{T},
        AbstractVector{T},
        AbstractVector{UInt64},
        AbstractVector{UInt32},
    }
    return P(ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
end

state(ps) =
    (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m)

function check_values(actual, exact, magnitude, count)
    T = eltype(actual)
    for i in eachindex(actual)
        @test isfinite(actual[i])
        @test abs(O.rat(actual[i]) - exact[i]) <= O.error_bound(T, exact[i], magnitude[i], count)
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

@testset "Canonical unrounded Yee charge and current values" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64), shape in (CIC(), YeeTriangularCloud(T(0.25)))
        g = YeeGrid((2,), (one(T),))
        ps = population(T, T[0, 0.5])
        ps.weight .= T[1, 0.5]
        species, old, new = (ps,), (copy(ps.x[1]),), (T[0.5, 0],)
        original = (state(ps), map(copy, old), map(copy, new))
        work = YeeChargeWorkspace(g)
        for mode in (Val(:charge), Val(:impulse))
            charging = mode isa Val{:charge}
            before, after = charging ? (nothing, nothing) : (old, new)
            entries = charging ? 3 : 2
            expected = if shape isa CIC
                first(
                    charging ? O.charge(species, old, 2, one(T)) :
                    O.impulse(species, old, new, 2, one(T)),
                )
            else
                charging ? O.cloud_charge(species, old, 2, one(T), shape.half_width) :
                O.cloud_impulse(species, old, new, 2, one(T), shape.half_width)
            end
            raw = @inferred H._yee_charge_exact_sums(species, before, after, g, mode, shape)
            @test length(raw) == 3
            @test raw[1:entries] == expected
            if shape isa YeeTriangularCloud
                # These physical-volume fractions cannot be recovered from rounded current.
                @test any(value -> O.rat(T(value)) != value, expected)
            end
            saved = copy(raw)
            @test H._yee_charge_store_exact_candidate!(work, raw, entries) === nothing
            @test isequal(work.candidate[1:entries], T.(expected))
            @test raw == saved
            @test H._yee_charge_exact_candidate!(
                work,
                species,
                before,
                after,
                g,
                mode,
                entries,
                shape,
            ) === nothing
            @test isequal(work.candidate[1:entries], T.(expected))
            bad = copy(raw)
            bad[entries] = 2 * O.rat(floatmax(T))
            @test_throws ArgumentError H._yee_charge_store_exact_candidate!(work, bad, entries)
            @test isequal((state(ps), old, new), original)
        end
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
            impulse,
            Any[unit, unit],
            Any[T[0], T[0]],
            Any[T[1], T[NaN]],
            g,
            work,
        )
        @test impulse == original_impulse
        @test (q0, work.sum, work.candidate) == before
        @test_throws ArgumentError yee_particle_charge!(work.candidate, unit, g, work)
        @test (q0, work.sum, work.candidate) == before
    end
end

@testset "Longitudinal path voltage and current work" begin
    for T in (Float32, Float64)
        g = YeeGrid((3,), (T(0.1),))
        field = T[2, -1, 3]
        original = copy(field)
        L = g.L[1]
        for (a, b) in (
            (zero(T), L),
            (L, zero(T)),
            (L / 16, L / 8),
            (L / 8, 3L / 8),
            (L / 8, 7L / 8),
            (zero(T), zero(T)),
            (L / 2, L / 2),
            (L, L),
            (zero(T), nextfloat(zero(T))),
        )
            exact = O.voltage(field, a, b, 3, L, g.dx[1])
            voltage = yee_longitudinal_voltage(field, a, b, g)
            @test O.rounded(voltage, exact)
            @test O.rounded(yee_longitudinal_voltage(view(field, :), b, a, g), -exact)
        end
        @test isequal(field, original)
        # The stored Maxwell metric is not the exact physical cell length.
        # A midpoint product makes omitting this factor change the rounding cell.
        metric = 3 * O.rat(g.dx[1]) / O.rat(L)
        @test metric != 1
        constant = metric > 1 ? one(T) + 3eps(T) : nextfloat(one(T))
        displacement = T(3) / T(64)
        dc = fill(constant, 3)
        exact = metric * O.rat(constant) * O.rat(displacement)
        @test O.rounded(yee_longitudinal_voltage(dc, zero(T), displacement, g), exact)
        @test !O.rounded(constant * displacement, exact)

        # A full path has exactly unit normalized overlap on every edge, so
        # both independently rounded APIs evaluate the same dyadic field work.
        ps = population(T, T[0])
        impulse = zeros(T, 3)
        yee_longitudinal_impulse!(impulse, ps, T[0], T[L], g, YeeChargeWorkspace(g))
        electric = (field, zeros(T, 4), zeros(T, 4))
        source = (impulse, zeros(T, 4), zeros(T, 4))
        @test isequal(
            yee_longitudinal_voltage(field, zero(T), L, g),
            yee_source_work(electric, electric, source, g),
        )

        wide = YeeGrid((4,), (T(8),))
        extremes = T[floatmax(T), -floatmax(T), floatmax(T), -floatmax(T)]
        @test iszero(yee_longitudinal_voltage(extremes, zero(T), T(8), wide))
        @test_throws ArgumentError yee_longitudinal_voltage(
            fill(floatmax(T), 4),
            zero(T),
            T(8),
            wide,
        )
        @test_throws ArgumentError yee_longitudinal_voltage(field, -nextfloat(zero(T)), L, g)
        @test_throws ArgumentError yee_longitudinal_voltage(field, zero(T), nextfloat(L), g)
        @test_throws ArgumentError yee_longitudinal_voltage(field, T(NaN), L, g)
        @test_throws ArgumentError yee_longitudinal_voltage(field, 0, L, g)
        @test_throws ArgumentError yee_longitudinal_voltage([1, 2, 3], zero(T), L, g)
        @test_throws ArgumentError yee_longitudinal_voltage(T[NaN, 0, 0], zero(T), zero(T), g)
        @test_throws DimensionMismatch yee_longitudinal_voltage(field[1:2], zero(T), L, g)
        @test_throws ArgumentError yee_longitudinal_voltage(
            field,
            zero(T),
            L,
            YeeGrid((3,), (L,); boundaries = (:periodic,)),
        )
        @test_throws ArgumentError yee_longitudinal_voltage(
            field,
            zero(T),
            L,
            YeeGrid((3, 3), (L, L)),
        )
        @test isequal(field, original)
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
        voltage_field = ones(T, 8)
        voltage =
            in_mode(() -> yee_longitudinal_voltage(voltage_field, zero(T), sub, g), hardware, flush)
        @test O.rounded(voltage, exact)
        @test_throws ArgumentError in_mode(
            () -> yee_longitudinal_voltage(voltage_field, -sub, zero(T), g),
            hardware,
            flush,
        )
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

@testset "Public Yee impulse with strided reinterpret endpoints" begin
    for T in (Float32, Float64)
        g = YeeGrid((4,), (one(T),))
        ps = population(T, T[0.5, 0.5])
        target = T[0.25, 0.75]
        U = T === Float32 ? UInt32 : UInt64
        same_size = reinterpret(T, U[reinterpret(U, target[1]), reinterpret(U, target[2])])
        parent = Complex{T}[Complex{T}(target[1], target[2])]
        scalar_view = view(parent, 1)
        @test ndims(scalar_view) == 0
        reshaped = reinterpret(reshape, T, scalar_view)
        for after in (same_size, reshaped), collection in (ParticleSet{1,T}[ps], Any[ps])
            @test after isa StridedVector{T}
            @test collect(after) == target
            old, new = AbstractVector{T}[T[0, 1]], AbstractVector{T}[after]
            expected = O.impulse(collection, old, (target,), g.n[1], g.L[1])
            work = YeeChargeWorkspace(g)
            out, original = fill(T(7), g.n[1]), copy(after)
            @test yee_longitudinal_impulse!(out, collection, old, new, g, work) === out
            check_values(out, expected..., nparticles(ps))
            @test after == original
        end

        # A scalar view drops its dimension; a range view retains it. The
        # reshape-reinterpret constructor consequently makes a matrix here.
        ranged = reinterpret(reshape, T, view(parent, 1:1))
        @test size(ranged) == (2, 1)
        for collection in (ParticleSet{1,T}[ps], Any[ps])
            work = YeeChargeWorkspace(g)
            out = fill(T(7), g.n[1])
            work.sum .= 11.0
            work.candidate .= T(13)
            before = (copy(out), copy(work.sum), copy(work.candidate))
            @test_throws ArgumentError yee_longitudinal_impulse!(
                out,
                collection,
                (T[0, 1],),
                (ranged,),
                g,
                work,
            )
            @test (out, work.sum, work.candidate) == before
        end
    end
end

@testset "Yee scalar domain through broad public collections" begin
    for T in (Float32, Float64)
        G, W = YeeGrid{1,T}, YeeChargeWorkspace{T}
        report = JET.report_call(
            yee_longitudinal_voltage,
            Tuple{Vector{T},T,T,G};
            sourceinfo = :full,
            stacktrace_types_limit = 0,
        )
        if !isempty(JET.get_reports(report))
            show(IOContext(stdout, :limit => false), MIME"text/plain"(), report)
        end
        @test isempty(JET.get_reports(report))
        for P in (Vector{ParticleSet{1,T}}, Vector{Any})
            for (operation, types) in (
                (yee_particle_charge!, Tuple{Vector{T},P,G,W}),
                (
                    yee_longitudinal_impulse!,
                    Tuple{Vector{T},P,Vector{AbstractVector{T}},Vector{AbstractVector{T}},G,W},
                ),
            )
                report = JET.report_call(
                    operation,
                    types;
                    sourceinfo = :full,
                    stacktrace_types_limit = 0,
                )
                if !isempty(JET.get_reports(report))
                    show(IOContext(stdout, :limit => false), MIME"text/plain"(), report)
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end

@testset "Optional physical-width triangular Yee cloud" begin
    for T in (Float32, Float64)
        g = YeeGrid((4,), (one(T),))
        # Width precision is an original operand, independent of grid precision.
        a = T === Float32 ? Float64(0.3) : Float32(0.3)
        cloud = YeeTriangularCloud(a)
        ps = population(T, T[0, 0.25, 0.5, 0.5, 1]; q = T(2))
        second = abstract_storage(population(T, reverse(ps.x[1]); q = -one(T)))
        second.x = (view(second.x[1], :),)
        species = Any[ps, second]
        original = map(state, species)
        old = [copy(p.x[1]) for p in species]
        new = [T[0.125, 0.375, 0.5, 0.75, 0.875], T[1, 0.75, 0.625, 0.25, 0]]
        q0 = O.cloud_charge(species, old, 4, g.L[1], a)
        q1 = O.cloud_charge(species, new, 4, g.L[1], a)
        flux = O.cloud_impulse(species, old, new, 4, g.L[1], a)
        @test any(!iszero, flux)
        @test sum(q0) == sum(O.rat(p.q) * sum(O.rat, p.weight) for p in species)
        @test q1 - q0 + O.incidence(flux) == zeros(O.Q, 5)
        work = YeeChargeWorkspace(g)
        charge, impulse = zeros(T, 5), zeros(T, 4)
        yee_particle_charge!(charge, species, g, work; shape = cloud)
        yee_longitudinal_impulse!(impulse, species, old, new, g, work; shape = cloud)
        @test all(O.rounded(v, q) for (v, q) in zip(charge, q0))
        @test all(O.rounded(v, q) for (v, q) in zip(impulse, flux))
        reversed = similar(impulse)
        yee_longitudinal_impulse!(reversed, species, new, old, g, work; shape = cloud)
        @test reversed == -impulse
        field = T[1, -2, 3, 4]
        particle_work, work_error = zero(O.Q), zero(O.Q)
        for (pop, x0, x1) in zip(species, old, new), p in eachindex(pop.weight)
            reference =
                O.rat(g.dx[1]) *
                sum(O.rat(field[e+1]) * O.cloud_path(x0[p], x1[p], e, 4, g.L[1], a) for e = 0:3)
            voltage = yee_longitudinal_voltage(field, x0[p], x1[p], g; shape = cloud)
            @test O.rounded(voltage, reference)
            source = O.rat(pop.q) * O.rat(pop.weight[p])
            particle_work += source * O.rat(voltage)
            work_error += abs(source) * O.rounding_error(T, reference)
        end
        mesh_work = O.rat(g.dx[1]) * sum(O.rat(v) * O.rat(i) for (v, i) in zip(field, impulse))
        mesh_error =
            O.rat(g.dx[1]) * sum(abs(O.rat(field[e])) * O.rounding_error(T, flux[e]) for e = 1:4)
        @test abs(particle_work - mesh_work) <= work_error + mesh_error
        for x in T[0, prevfloat(T(0.5)), 0.5, nextfloat(T(0.5)), 1]
            expected =
                O.rat(g.dx[1]) *
                sum(O.rat(field[e+1]) * O.cloud_edge(x, e, 4, g.L[1], a) for e = 0:3)
            @test O.rounded(yee_longitudinal_field(field, x, g; shape = cloud), expected)
        end
        for (x, y) in ((zero(T), one(T)), (T(0.25), T(0.5)), (T(0.5), T(0.5)))
            expected =
                O.rat(g.dx[1]) *
                sum(O.rat(field[e+1]) * O.cloud_path(x, y, e, 4, g.L[1], a) for e = 0:3)
            @test O.rounded(yee_longitudinal_voltage(field, x, y, g; shape = cloud), expected)
        end
        @test isequal(map(state, species), original)
        # Compact interior support keeps its final charge at the following
        # node. Also exercise exact support endpoints on grid nodes.
        for (center, width) in ((T(0.375), T(0.0625)), (T(0.5), T(0.25)))
            interior = population(T, T[center])
            aligned = YeeTriangularCloud(width)
            yee_particle_charge!(charge, interior, g, work; shape = aligned)
            expected = O.cloud_charge((interior,), (interior.x[1],), 4, g.L[1], width)
            @test all(O.rounded(v, q) for (v, q) in zip(charge, expected))
        end
        # Global preflight still precedes scratch/output writes in the new mode.
        second.weight[end] = -one(T)
        before = (copy(charge), copy(work.sum), copy(work.candidate))
        @test_throws ArgumentError yee_particle_charge!(charge, species, g, work; shape = cloud)
        @test (charge, work.sum, work.candidate) == before
        second.weight[end] = one(T)
        @test_throws ArgumentError yee_particle_charge!(work.candidate, ps, g, work; shape = cloud)
        bad = copy(field)
        bad[end] = T(NaN)
        @test_throws ArgumentError yee_longitudinal_voltage(bad, T(0.5), T(0.5), g; shape = cloud)
        @test_throws ArgumentError yee_longitudinal_field(field, -one(T), g; shape = cloud)
        for width in T[0, -1, Inf, NaN]
            @test_throws ArgumentError YeeTriangularCloud(width)
        end
        # Exact whole-entry completion retains overflowing q*weight cancellation.
        positive = population(T, T[0.5]; q = floatmax(T), weight = T(2))
        negative = population(T, T[0.5]; q = -floatmax(T), weight = T(2))
        yee_particle_charge!(charge, (positive, negative), g, work; shape = cloud)
        @test all(iszero, charge)
        before = copy(charge)
        @test_throws ArgumentError yee_particle_charge!(charge, positive, g, work; shape = cloud)
        @test charge == before
        yee_particle_charge!(charge, ParticleSet{1,T}[], g, work; shape = cloud)
        @test all(iszero, charge)
        # Use the existing context helper; width admission itself must retain
        # positive stored subnormal bits under FTZ, before any field operation.
        tiny = nextfloat(zero(T))
        result = in_mode(RoundUp, true) do
            narrow = YeeTriangularCloud(tiny)
            yee_longitudinal_field(ones(T, 4), T(0.5), g; shape = narrow)
        end
        @test result === one(T)
    end
end

@testset "Cloud work metric, wall response and stationary Gram pairing" begin
    for T in (Float32, Float64)
        g = YeeGrid((11,), (T(0.1),))
        cloud = YeeTriangularCloud(T(0.01))
        constant = fill(T(3), 11)
        eta = 11 * O.rat(g.dx[1]) / O.rat(g.L[1])
        @test O.rounded(yee_longitudinal_field(constant, T(0.05), g; shape = cloud), 3eta)
        @test O.rounded(yee_longitudinal_field(constant, zero(T), g; shape = cloud), 3eta / 2)
        @test O.rounded(yee_longitudinal_field(constant, g.L[1], g; shape = cloud), 3eta / 2)
        # A midpoint fixture makes omission of the nonunit stored-dx metric
        # observable after final rounding, rather than merely checking its formula.
        metric_grid = YeeGrid((10,), (T(0.1),))
        metric = 10 * O.rat(metric_grid.dx[1]) / O.rat(metric_grid.L[1])
        @test metric != 1
        left = metric > 1 ? one(T) : nextfloat(one(T))
        right = nextfloat(left)
        sampled = zeros(T, 10)
        sampled[5], sampled[6] = left, right
        midpoint = (O.rat(left) + O.rat(right)) / 2
        U = T === Float32 ? UInt32 : UInt64
        omitted = iseven(reinterpret(U, left)) ? left : right
        expected = metric * midpoint
        @test !O.rounded(omitted, expected)
        @test O.rounded(
            yee_longitudinal_field(
                sampled,
                metric_grid.L[1] / 2,
                metric_grid;
                shape = YeeTriangularCloud(metric_grid.dx[1] / 4),
            ),
            expected,
        )
        # Reuse the retained unstable geometry, now with the matching cloud pair.
        n, length, width = 2, O.Q(1), O.Q(1 // 4)
        basis =
            [[O.cloud_edge(x, e, n, length, width) for e = 0:1] for x in (O.Q(1 // 8), O.Q(5 // 8))]
        matrix = [
            sum(w * b[e] * b[f] / 2 for (w, b) in zip((O.Q(1), O.Q(5 // 3)), basis)) for
            e = 1:2, f = 1:2
        ]
        @test matrix[1, 2] == matrix[2, 1]
        @test matrix[1, 1] >= 0 && matrix[2, 2] >= 0
        @test matrix[1, 1] * matrix[2, 2] >= matrix[1, 2]^2
        symmetric = T[-1, 1]
        two = YeeGrid((2,), (one(T),))
        for (x, b) in zip(T[0.125, 0.625], basis), e = 1:2
            unit = zeros(T, 2)
            unit[e] = one(T)
            @test O.rounded(
                yee_longitudinal_field(unit, x, two; shape = YeeTriangularCloud(T(0.25))),
                b[e] / 2,
            )
        end
        @test yee_longitudinal_field(
            symmetric,
            T(0.5),
            two;
            shape = YeeTriangularCloud(T(0.25)),
        ) === zero(T)
        # An exact volume quadrature checks the declared finite cloud wall law.
        a = O.Q(1 // 4)
        @test O.polynomial_integral(r -> r * O.density(r, O.Q(0), a), [O.Q(0), a]) == a / 6
        # Widths larger than L retain endpoint tails; the continuum a/6 wall
        # formula no longer applies. Check the actual sampled-field contraction.
        wide = YeeTriangularCloud(T(2))
        sampled = T[0.25, 0.75] # exact cell averages of E=r on this grid
        expected =
            O.rat(two.dx[1]) * sum(O.rat(sampled[e+1]) * O.cloud_edge(0, e, 2, 1, 2) for e = 0:1)
        @test O.rounded(yee_longitudinal_field(sampled, zero(T), two; shape = wide), expected)
        @test O.polynomial_integral(r -> r * O.density(r, O.Q(0), O.Q(2)), O.Q[0, 1]) == 1 // 6
    end
end

@testset "Fixed physical cloud width force bound" begin
    # Independent volume quadrature, not a mesh-shaped point-sheet substitute.
    x, a = O.rat(0.3), O.rat(0.13)
    physical(r) = -r + O.cloud_mass(r, x, a)
    force = O.polynomial_integral(r -> O.density(r, x, a) * physical(r), [x - a, x, x + a])
    for n in (32, 64, 128)
        g = YeeGrid((n,), (1.0,))
        h = O.Q(1) / n
        exact_field = [
            begin
                lo, hi = e * h, (e + 1) * h
                knots = [lo, hi]
                append!(knots, filter(r -> lo < r < hi, [x - a, x, x + a]))
                O.polynomial_integral(physical, knots) / h
            end for e = 0:(n-1)
        ]
        field = Float64.(exact_field)
        result =
            yee_longitudinal_field(field, Float64(x), g; shape = YeeTriangularCloud(Float64(a)))
        represented =
            O.rat(g.dx[1]) * sum(O.rat(field[e+1]) * O.cloud_edge(x, e, n, 1, a) for e = 0:(n-1))
        @test O.rounded(result, represented)
        eta = O.rat(g.dx[1]) / h
        spatial_bound = eta * (1 / a^2) * (1 + 1 / a) * h^2
        if n == 128
            @test spatial_bound < abs(force)
        end
        # Separately enclose initial cell-average publication and final rounding.
        field_error = maximum(O.rounding_error(Float64, v) for v in exact_field)
        @test abs(O.rat(result) - eta * force) <=
              spatial_bound + eta * field_error + O.rounding_error(Float64, represented)
    end
end
end
