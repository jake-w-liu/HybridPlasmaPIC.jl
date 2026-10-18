module IsothermalFirstPassChecks
using HybridPlasmaPIC, Test, Random
include("isothermal_energy/oracles.jl")
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
@testset "Degree-eight first pass and unresolved continuation" begin
    rng = MersenneTwister(920733)
    selected = Dict(true => 0, false => 0)
    for repetition = 1:32768
        d = 2.0^-8
        n = [1 + d, 1 - d]
        temperature = 0.5 + 1.5rand(rng)
        dx = (1.0,)
        state = H._isothermal_statistics(n)
        lower, upper, exponent =
            H._isothermal_complete_dyadic(n, temperature, dx, state, 69, Val(8))
        lo, hi = H._from_exact(Float64, lower, exponent), H._from_exact(Float64, upper, exponent)
        certified = isequal(lo, hi)
        selected[certified] >= 8 && continue
        selected[certified] += 1
        reflo, refhi = Oracle.oracle_integral(n, temperature, dx)
        @test Q(lower) * Q(2)^exponent <= reflo <= refhi <= Q(upper) * Q(2)^exponent
        result = H._isothermal_energy(n, temperature, dx)
        @test result !== nothing
        @test Oracle.rounded(result, reflo, refhi)
        @test isequal(result, H._isothermal_exact_energy(n, temperature, dx))
        @test !certified || isequal(result, lo)
        @test isequal(result, H._isothermal_energy(view(n, :), temperature, dx))
        all(==(8), values(selected)) && break
    end
    println((selected,))
    flush(stdout)
    @test selected[true] == 8
    @test selected[false] == 8
end
end
