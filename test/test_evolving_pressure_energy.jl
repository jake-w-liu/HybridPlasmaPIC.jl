using HybridPlasmaPIC, Test, Random

function ep_energy_exact(pressure, gamma, grid)
    Q = Rational{BigInt}
    return sum(Q, pressure; init = Q(0)) * prod(Q.(grid.dx)) / (Q(gamma) - 1)
end

# Check the independent exact value against the output's rounding cell. This
# does not call the production rational-to-float conversion or use BigFloat.
function ep_energy_rounding_cell(value::T, exact) where {T}
    Q = Rational{BigInt}
    @test isfinite(value) || value == T(Inf)
    if iszero(value)
        @test 0 <= exact <= Q(nextfloat(zero(T))) / 2
    elseif isinf(value)
        edge = Q(floatmax(T)) + (Q(floatmax(T)) - Q(prevfloat(floatmax(T)))) / 2
        @test exact >= edge
    elseif isfinite(value)
        lower = (Q(prevfloat(value)) + Q(value)) / 2
        upper =
            isfinite(nextfloat(value)) ? (Q(value) + Q(nextfloat(value))) / 2 :
            Q(value) + (Q(value) - Q(prevfloat(value))) / 2
        @test lower <= exact <= upper
        if exact == lower || exact == upper
            bits = T === Float32 ? reinterpret(UInt32, value) : reinterpret(UInt64, value)
            @test iseven(bits)
        end
    end
end

function ep_energy_check_public(pressure::Array{T,D}, gamma::T, grid) where {T,D}
    density = ones(T, grid.n)
    field = ntuple(_ -> zeros(T, grid.n), 3)
    closure = EvolvingElectrons(gamma)
    exact = ep_energy_exact(pressure, gamma, grid)
    before = deepcopy((pressure, density, field))
    result = electron_internal_energy(density, closure, grid; pe = pressure)
    ep_energy_rounding_cell(result, exact)
    @test isequal(electron_internal_energy(density, field, closure, grid; pe = pressure), result)
    pressure_view = view(pressure, ntuple(_ -> Colon(), D)...)
    @test isequal(electron_internal_energy(density, closure, grid; pe = pressure_view), result)
    ions = ParticleSet{D,T}(0)
    for populations in (ions, [ions])
        budget = energy_budget(populations, field, density, closure, grid; pe = pressure)
        @test isequal(budget.electron_internal, result)
    end
    @test isequal((pressure, density, field), before)
    return result
end

@testset "Stored pressure complete energy and rounding boundaries" begin
    for T in (Float32, Float64)
        for D = 1:3
            shape = ntuple(_ -> 2, D)
            grid = FourierGrid(shape, ntuple(_ -> one(T), D))
            ep_energy_check_public(fill(floatmax(T), shape), T(3), grid)
            ep_energy_check_public(fill(floatmax(T), shape), nextfloat(one(T)), grid)
        end

        # Normal additions can each round back to floatmax while their exact
        # total reaches or exceeds the nearest-even overflow midpoint.
        quarter = (floatmax(T) - prevfloat(floatmax(T))) / T(4)
        grid = FourierGrid((3,), (T(3),))
        for tail in (prevfloat(quarter), quarter, nextfloat(quarter))
            ep_energy_check_public(T[floatmax(T), quarter, tail], T(2), grid)
        end

        pressure = T[nextfloat(zero(T)), 0]
        # Below, at, and above the exact halfway point of the zero/qmin cell.
        for spacing in (prevfloat(one(T)), one(T), nextfloat(one(T)), T(4))
            grid = FourierGrid((2,), (2spacing,))
            ep_energy_check_public(pressure, T(3), grid)
        end

        for (length, level) in (
            (T === Float32 ? T(1e14) : T(2e103), T === Float32 ? T(1e-42) : T(1e-310)),
            (T === Float32 ? T(1e-16) : T(1e-110), T === Float32 ? T(1e30) : T(1e300)),
        )
            grid = FourierGrid((2, 2, 2), (length, length, length))
            ep_energy_check_public(fill(level, grid.n), T(2), grid)
            ep_energy_check_public(zeros(T, grid.n), T(2), grid)
        end

        # One-cell dimensions isolate volume association without a nonzero
        # Fourier mode; all three spacings are positive and representable.
        small_exponent = T === Float32 ? -72 : -540
        large_exponent = T === Float32 ? 100 : 1000
        small = ldexp(nextfloat(one(T)), small_exponent)
        large = ldexp(one(T), large_exponent)
        for lengths in ((small, small, large), (large, small, small), (small, large, small))
            grid = FourierGrid((1, 1, 1), lengths)
            ep_energy_check_public(fill(T(3), grid.n), T(2), grid)
        end

        grid = FourierGrid((8,), (T(2),))
        mixed = T[floatmax(T), floatmax(T), nextfloat(zero(T)), 0, 1, 2, 3, 4]
        ep_energy_check_public(mixed, T(5), grid)
    end
end

@testset "Stored pressure overflow upper bound vs exact integrals" begin
    random = MersenneTwister(0xB0A1D)
    Q = Rational{BigInt}
    for T in (Float32, Float64), D = 1:3, _ = 1:32
        count = rand(random, (1, 3, 8, 31, 128))
        # These ranges keep the rounded positive sum finite; spacing products
        # may overflow/underflow and are evaluated only by the rational oracle.
        exponents = T === Float32 ? (-120:110) : (-1000:1000)
        pressure = [ldexp(one(T) + rand(random, T), rand(random, exponents)) for _ = 1:count]
        rounded = sum(pressure)
        gamma = rand(random, (nextfloat(one(T)), T(5 / 3), T(2), T(100), floatmax(T)))
        spacings = ntuple(_ -> ldexp(one(T) + rand(random, T), rand(random, exponents)), D)
        upper = HybridPlasmaPIC._pressure_energy_upper_exponent(rounded, gamma, spacings, count)
        bound =
            upper >= 0 ? (BigInt(1) << Int(upper)) // BigInt(1) :
            BigInt(1) // (BigInt(1) << Int(-upper))
        exact = sum(Q, pressure; init = Q(0)) * prod(Q.(spacings)) / (Q(gamma) - 1)
        @test all(p -> p <= rounded, pressure)
        @test exact <= bound
        if upper <= exponent(floatmax(T))
            edge = Q(floatmax(T)) + (Q(floatmax(T)) - Q(prevfloat(floatmax(T)))) / 2
            @test exact < edge
        end
    end
end

function ep_energy_allocation(density, closure, grid, pressure)
    return @allocated electron_internal_energy(density, closure, grid; pe = pressure)
end

@testset "Stored pressure ordinary arithmetic, validation and allocation" begin
    random = MersenneTwister(0xE1EC7)
    for T in (Float32, Float64), D = 1:3
        shape = ntuple(_ -> 8, D)
        grid = FourierGrid(shape, ntuple(d -> T(2 + d / 4), D))
        density = ones(T, shape)
        pressure = rand(random, T, shape) .+ T(0.2)
        field = ntuple(_ -> zeros(T, shape), 3)
        ions = ParticleSet{D,T}(0)
        for gamma in (nextfloat(one(T)), T(5 / 3), T(3), T(100))
            closure = EvolvingElectrons(gamma)
            original = sum(pressure) / (gamma - one(T)) * prod(grid.dx)
            actual = electron_internal_energy(density, closure, grid; pe = pressure)
            @test isequal(actual, original)
            ep_energy_allocation(density, closure, grid, pressure)
            @test ep_energy_allocation(density, closure, grid, pressure) == 0
            @test energy_budget(ions, field, density, closure, grid; pe = pressure).total == actual
        end
        closure = EvolvingElectrons(T(5 / 3))
        for value in (T(-1), T(Inf), T(NaN))
            invalid = copy(pressure)
            invalid[1] = value
            before = copy(invalid)
            @test_throws ArgumentError electron_internal_energy(
                density,
                closure,
                grid;
                pe = invalid,
            )
            @test isequal(invalid, before)
        end
    end
end
