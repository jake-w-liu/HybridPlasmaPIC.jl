using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

# Exponent limits are process-wide MPFR state. This fixture is only run in a
# dedicated process with one Julia worker thread and its own restored settings.
@assert Threads.nthreads() == 1
@assert length(ARGS) == 1
@assert realpath(pkgdir(H)) == realpath(only(ARGS))

dyadic(e::Integer) = e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << -e)

function stored_rational(value::BigFloat)
    significand, power, sign = Base.decompose(value)
    return (BigInt(sign) * significand) * dyadic(power)
end

get_emin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
get_emax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
set_emin(value) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), value)
set_emax(value) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), value)

mutable struct EndpointRNG <: AbstractRNG
    uniforms::Int
    normals::Int
end
Random.rand(rng::EndpointRNG, ::Type{BigFloat}) = (rng.uniforms += 1; BigFloat(0.5))
Random.randn(rng::EndpointRNG, ::Type{BigFloat}) = (rng.normals += 1; BigFloat(0))
population_state(p) = (map(copy, (p.x..., p.v..., p.weight, p.id, p.tag)), p.q, p.m)

function check_public_threshold(speed, speed_exact, sign, bits)
    electrons = ParticleSet{1,BigFloat}(1; q = BigFloat(-1), m = BigFloat(1))
    ions = ParticleSet{1,BigFloat}(1; q = BigFloat(1), m = BigFloat(100))
    electrons.x[1][1], ions.x[1][1] = BigFloat(0.25), BigFloat(0.75)
    electrons.v[1][1] = sign * speed
    electrons.v[2][1], electrons.v[3][1] = -BigFloat(0), BigFloat(0)
    electrons.weight[1] = BigFloat(2)
    electrons.id[1], ions.id[1] = UInt64(7), UInt64(11)
    electrons.tag[1], ions.tag[1] = 3, 4
    before = population_state(electrons), population_state(ions)
    ec, ic = Ref(UInt64(20)), Ref(UInt64(30))
    rng = EndpointRNG(0, 0)
    # m=1 and E_iz=1/2 give strict eligibility speed^2>1. For these
    # neighboring speeds, rate=1 and dt=16 put probability well above 1/2.
    selected = speed_exact > 1
    count = ionize_mcc!(
        electrons,
        ions,
        BigFloat(16);
        nσ_iz = BigFloat(1),
        E_iz = BigFloat(0.5),
        e_nextid = ec,
        i_nextid = ic,
        rng,
    )
    @test count == Int(selected)
    @test (nparticles(electrons), nparticles(ions)) == (1 + Int(selected), 1 + Int(selected))
    @test (rng.uniforms, rng.normals) == (Int(selected), 6 * Int(selected))
    @test (ec[], ic[]) == (UInt64(20 + Int(selected)), UInt64(30 + Int(selected)))
    if selected && count == 1
        actual = electrons.v[1][1]
        @test signbit(actual) == (sign < 0)
        @test isfinite(actual) && !iszero(actual)
        @test precision(actual) == bits
        # The complete cooled component square simplifies exactly to v^2-1.
        squared = speed_exact^2 - 1
        magnitude = abs(actual)
        value = stored_rational(magnitude)
        previous, following =
            stored_rational(prevfloat(magnitude)), stored_rational(nextfloat(magnitude))
        lower, upper = (previous + value) / 2, (value + following) / 2
        significand = value / dyadic(exponent(magnitude) - bits + 1)
        @test denominator(significand) == 1
        even = iseven(numerator(significand))
        @test lower^2 < squared < upper^2 || (even && (squared == lower^2 || squared == upper^2))
        @test isequal(electrons.v[2][1], -BigFloat(0))
        @test isequal(electrons.v[3][1], BigFloat(0))
        @test electrons.id == UInt64[7, 20]
        @test ions.id == UInt64[11, 30]
        @test electrons.weight == BigFloat[2, 2]
        @test ions.weight == BigFloat[1, 2]
        @test electrons.x[1] == BigFloat[0.25, 0.25]
        @test ions.x[1] == BigFloat[0.75, 0.25]
        @test all(iszero(electrons.v[c][2]) && iszero(ions.v[c][2]) for c = 1:3)
        @test all(isequal(ions.v[c][1], before[2][1][1+c][1]) for c = 1:3)
        @test electrons.tag[1] == 3 && ions.tag[1] == 4
    elseif !selected
        @test isequal((population_state(electrons), population_state(ions)), before)
    end
end

@testset "Isolated BigFloat square-root endpoints" begin
    saved_emin, saved_emax = get_emin(), get_emax()
    saved_bits, saved_mode = precision(BigFloat), rounding(BigFloat)
    stored_wide = setprecision(BigFloat, 256) do
        setrounding(BigFloat, RoundNearest) do
            value = BigFloat(1) + ldexp(BigFloat(1), -192)
            (value, -value)
        end
    end
    try
        @assert set_emin(-64) == 0
        @assert set_emax(64) == 0
        @test (get_emin(), get_emax()) == (-64, 64)
        for bits in (24, 53, 96, 128)
            # MPFR uses normalized significands only: the minimum is 2^(emin-1).
            # The largest p-bit value is (2^p-1)*2^(emax-p).
            minimum = dyadic(-65)
            quantum_at_maximum = dyadic(64 - bits)
            maximum = ((BigInt(1) << bits) - 1) * quantum_at_maximum
            half_minimum_squared = (minimum / 2)^2
            underflow_delta = half_minimum_squared * dyadic(-bits - 10)
            overflow = maximum + quantum_at_maximum / 2
            overflow_delta = quantum_at_maximum * dyadic(-bits - 10)
            zero_exact, two_exact = BigInt(0) // BigInt(1), BigInt(2) // BigInt(1)
            cases = [
                ("zero", zero_exact, zero_exact),
                ("below half minimum", half_minimum_squared - underflow_delta, zero_exact),
                ("half minimum tie", half_minimum_squared, zero_exact),
                ("above half minimum", half_minimum_squared + underflow_delta, minimum),
                ("minimum", minimum^2, minimum),
                ("ordinary square", BigInt(4) // BigInt(1), two_exact),
                ("maximum", maximum^2, maximum),
                ("below overflow", (overflow - overflow_delta)^2, maximum),
                ("overflow tie", overflow^2, nothing),
                ("above overflow", (overflow + overflow_delta)^2, nothing),
            ]
            setprecision(BigFloat, bits) do
                for mode in (RoundNearest, RoundUp, RoundDown)
                    setrounding(BigFloat, mode) do
                        @testset "bits=$bits mode=$mode case=$name" for (name, exact, expected) in
                                                                        cases
                            actual = H._sqrt_from_exact(BigFloat, exact)
                            @test actual isa BigFloat
                            @test precision(actual) == bits
                            @test !signbit(actual)
                            if expected === nothing
                                @test isinf(actual)
                            else
                                @test isfinite(actual)
                                @test stored_rational(actual) == expected
                            end
                            @test precision(BigFloat) == bits
                            @test rounding(BigFloat) == mode
                            @test (get_emin(), get_emax()) == (-64, 64)
                        end
                        @testset "Stored conversion bits=$bits mode=$mode" begin
                            # Expected rationals are constructed from the test
                            # values, without calling the decomposition oracle.
                            pairs = (
                                (BigFloat(0), zero_exact),
                                (-BigFloat(0), zero_exact),
                                (BigFloat(2), two_exact),
                                (BigFloat(-2), -two_exact),
                                (BigFloat(0.5), dyadic(-1)),
                                (BigFloat(-0.5), -dyadic(-1)),
                                (ldexp(BigFloat(1), -65), minimum),
                                (floatmax(BigFloat), maximum),
                                (stored_wide[1], 1 + dyadic(-192)),
                                (stored_wide[2], -1 - dyadic(-192)),
                            )
                            for (value, expected) in pairs
                                @test H._exact_rational(value) == expected
                                @test precision(BigFloat) == bits
                                @test rounding(BigFloat) == mode
                                @test (get_emin(), get_emax()) == (-64, 64)
                            end
                        end
                        @testset "Public threshold bits=$bits mode=$mode" begin
                            speeds = (
                                (prevfloat(BigFloat(1)), 1 - dyadic(-bits)),
                                (BigFloat(1), BigInt(1) // BigInt(1)),
                                (nextfloat(BigFloat(1)), 1 + dyadic(1 - bits)),
                            )
                            for (speed, speed_exact) in speeds, sign in (-1, 1)
                                check_public_threshold(speed, speed_exact, sign, bits)
                                @test precision(BigFloat) == bits
                                @test rounding(BigFloat) == mode
                                @test (get_emin(), get_emax()) == (-64, 64)
                            end
                        end
                    end
                end
            end
        end
    finally
        @test set_emin(saved_emin) == 0
        @test set_emax(saved_emax) == 0
    end
    @test (get_emin(), get_emax()) == (saved_emin, saved_emax)
    @test precision(BigFloat) == saved_bits
    @test rounding(BigFloat) == saved_mode
end
println("BigFloat endpoint settings restored")
