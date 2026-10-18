module ShockFluxRangeTests
using HybridPlasmaPIC, Test, Random
include("legacy.jl")

const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function independent_round(::Type{T}, q::Q) where {T}
    setprecision(BigFloat, 4096) do
        lower = setrounding(BigFloat, RoundDown) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        upper = setrounding(BigFloat, RoundUp) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        @test Q(lower) <= q <= Q(upper)
        # Destination conversion is nearest-even regardless of the endpoint
        # construction mode. Julia's default BigFloat conversion inherits it.
        lo, hi = T(lower, RoundNearest), T(upper, RoundNearest)
        @test isequal(lo, hi)
        return lo
    end
end

function population(T, velocities, weights)
    p = ParticleSet{1,T}(length(weights))
    p.x[1] .= zero(T)
    foreach(a -> fill!(a, zero(T)), p.v)
    p.v[1] .= velocities
    p.weight .= weights
    return p
end
function exact_nodes(sh, ps)
    out = fill(zero(Q), sh.s.n)
    for p in eachindex(ps.weight)
        # Independent stencil for dyadic witness positions. It agrees with
        # the declared folded linear stencil without calling the production helper.
        x = Q(ps.x[1][p]) / Q(sh.s.dx)
        left = floor(Int, x)
        fraction = x - left
        a = clamp(left + 1, 1, length(out))
        b = clamp(left + 2, 1, length(out))
        term =
            Q(ps.weight[p]) *
            (Q(ps.v[1][p])^2 + Q(ps.v[2][p])^2 + Q(ps.v[3][p])^2) *
            Q(ps.v[1][p]) / 2
        out[a] += (1 - fraction) * term
        out[b] += fraction * term
    end
    return [q / Q(h) for (q, h) in zip(out, sh.s.H)]
end
function allocations(sh::S, ps::P) where {S,P}
    Legacy._ion_energy_flux(sh, ps)
    H._ion_energy_flux(sh, ps)
    return (@allocated(Legacy._ion_energy_flux(sh, ps)), @allocated(H._ion_energy_flux(sh, ps)))
end
println((VERSION,))
flush(stdout)
@testset "Complete scalar and nodal shock flux correction" begin
    for T in (Float32, Float64)
        sh = PerpShock(8, T(7); γe = T(1.5), η = zero(T))
        ranges = T === Float32 ? ((70, -100), (-80, 120)) : ((600, -1000), (-600, 1000))
        for (ve, ne) in ranges, sign in (-1, 1)
            v = sign * ldexp(one(T), ve)
            n = ldexp(one(T), ne)
            expected = independent_round(T, Q(n) * Q(v)^3 / 2)
            @test isequal(H._shock_bulk_energy_flux(n, v, zero(T)), expected)
            ps = population(T, T[v], T[n])
            before = deepcopy((ps.x, ps.v, ps.weight))
            exact = exact_nodes(sh, ps)
            actual = H._ion_energy_flux(sh, ps)
            @test all(isequal(actual[i], independent_round(T, exact[i])) for i in eachindex(actual))
            @test isequal((ps.x, ps.v, ps.weight), before)
        end
        for sign in (-1, 1)
            gamma = T(1.5)
            pressure = floatmax(T) / 2
            u = T(sign * 0.25)
            expected = independent_round(T, Q(gamma) * Q(pressure) * Q(u) / (Q(gamma) - 1))
            @test isequal(H._shock_electron_enthalpy_flux(gamma, pressure, u), expected)
        end
        q = nextfloat(zero(T))
        v = ldexp(one(T), T === Float32 ? 70 : 600)
        for ps in (
            population(T, ones(T, 2), fill(q, 2)),
            population(T, T[v, -v, 1], ones(T, 3)),
            population(T, T[], T[]),
        )
            exact = exact_nodes(sh, ps)
            actual = H._ion_energy_flux(sh, ps)
            @test all(isequal(actual[i], independent_round(T, exact[i])) for i in eachindex(actual))
        end
        scale = ldexp(one(T), T === Float32 ? 100 : 800)
        large = PerpShock(8, T(7) * scale)
        v = ldexp(one(T), T === Float32 ? 50 : 400)
        ps = population(T, T[v], ones(T, 1))
        actual = H._ion_energy_flux(large, ps)
        exact = exact_nodes(large, ps)
        @test all(isequal(actual[i], independent_round(T, exact[i])) for i in eachindex(actual))
    end
end
@testset "Ordinary scalar compatibility and public nodal resources" begin
    rng = MersenneTwister(7433)
    for T in (Float32, Float64)
        for _ = 1:100
            n = T(rand(rng))
            ux = T(randn(rng))
            uy = T(randn(rng))
            gamma = T(5 / 3)
            pressure = T(rand(rng))
            @test isequal(H._shock_bulk_energy_flux(n, ux, uy), T(0.5) * n * (ux^2 + uy^2) * ux)
            @test isequal(
                H._shock_electron_enthalpy_flux(gamma, pressure, ux),
                (gamma / (gamma - one(T))) * pressure * ux,
            )
        end
        for count in (0, 1, 8, 64), nodes in (8, 64)
            sh = PerpShock(nodes, T(nodes - 1))
            ps = population(T, T.(randn(rng, count)), T.(rand(rng, count)))
            ps.x[1] .= T.(rand(rng, count) .* (nodes - 1))
            expected = Legacy._ion_energy_flux(sh, ps)
            actual = H._ion_energy_flux(sh, ps)
            @test isequal(actual, expected)
            # Concrete helpers prevent heterogeneous harness-loop allocation.
            old_bytes, new_bytes = allocations(sh, ps)
            @test new_bytes == old_bytes
            println((T, count, nodes, old_bytes, new_bytes))
            flush(stdout)
        end
    end
end
@testset "Public boundary flux components and lattice certificate" begin
    for T in (Float32, Float64), sign in (-1, 1)
        sh = PerpShock(8, T(7); γe = T(1.5))
        ve, ne = T === Float32 ? (70, -100) : (600, -1000)
        velocity, density = sign * ldexp(one(T), ve), ldexp(one(T), ne)
        sh.ux .= velocity
        sh.n .= density
        expected = sign * ldexp(one(T), 3ve + ne - 1)
        @test boundary_energy_flux(sh).kinetic == (expected, expected)
        ps = population(T, T[velocity], T[density])
        @test boundary_energy_flux(sh; ps = ps).kinetic == (zero(T), 2expected)
        sh.ux .= T(sign * 0.25)
        sh.pe .= floatmax(T) / 2
        expected = independent_round(T, Q(floatmax(T)) * 3sign / 8)
        @test boundary_energy_flux(sh).enthalpy == (expected, expected)
        @test boundary_energy_flux(sh; ps = ps).enthalpy == (expected, expected)
    end
    for T in (Float32, Float64)
        sh = PerpShock(8, T(7))
        large = ldexp(one(T), T === Float32 ? 20 : 200)
        ps = population(T, T[1, large], ones(T, 2))
        ps.x[1][2] = T(7)
        # Both nonzero stencils and all input factors are integer powers of
        # two. The half-flux factor sets the minimum term quantum to 2^-1.
        @test H._shock_flux_quantum(sh, ps) == -1
        @test isequal(H._ion_energy_flux(sh, ps), Legacy._ion_energy_flux(sh, ps))
    end
end
println("SHOCK_FLUX_INITIAL_CHECKS_COMPLETE")
flush(stdout)

end
