module ShockFluxPropertiesTests
using HybridPlasmaPIC, Test, Random


const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function rounded(::Type{T}, q::Q) where {T}
    setprecision(BigFloat, 8192) do
        lo = setrounding(BigFloat, RoundDown) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        hi = setrounding(BigFloat, RoundUp) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        @test Q(lo) <= q <= Q(hi)
        a, b = T(lo, RoundNearest), T(hi, RoundNearest)
        @test isequal(a, b)
        a
    end
end

function finite_bits(rng, ::Type{T}) where {T}
    U = T === Float32 ? UInt32 : UInt64
    while true
        x = reinterpret(T, rand(rng, U))
        isfinite(x) && return x
    end
end

function exact_nodes(sh, ps)
    out = fill(zero(Q), sh.s.n)
    for p in eachindex(ps.weight)
        # Fixtures use dyadic positions and unit spacing, making the declared
        # folded linear stencil exact without calling the implementation.
        position = Q(ps.x[1][p]) / Q(sh.s.dx)
        cell = floor(Int, position)
        fraction = position - cell
        a, b = clamp(cell + 1, 1, sh.s.n), clamp(cell + 2, 1, sh.s.n)
        flux = Q(ps.weight[p]) * sum(Q(v[p])^2 for v in ps.v) * Q(ps.v[1][p]) / 2
        out[a] += (1 - fraction) * flux
        out[b] += fraction * flux
    end
    [value / Q(norm) for (value, norm) in zip(out, sh.s.H)]
end

println((VERSION, Threads.nthreads()))
flush(stdout)
@testset "Dyadic quantum is exact and minimal" begin
    rng = MersenneTwister(324813)
    for T in (Float32, Float64), k = 1:600
        x = k == 1 ? nextfloat(zero(T)) : k == 2 ? floatmax(T) : finite_bits(rng, T)
        iszero(x) && continue
        exponent = H._shock_flux_factor_quantum(x)
        quantum = exponent >= 0 ? Q(BigInt(1) << exponent) : Q(1, BigInt(1) << -exponent)
        integer = Q(x) / quantum
        @test denominator(integer) == 1
        @test isodd(numerator(integer))
    end
end

@testset "Complete cold nodal expressions, folds, ties and overflow" begin
    rng = MersenneTwister(774231)
    for T in (Float32, Float64)
        huge = ldexp(one(T), T === Float32 ? 70 : 600)
        tiny = nextfloat(zero(T))
        for k = 1:120
            sh = PerpShock(8, T(7))
            count = 3 + mod(k, 10)
            ps = ParticleSet{1,T}(count)
            ps.x[1] .= T.(rand(rng, -8:64, count)) ./ 8
            ps.weight .= T.(rand(rng, -16:16, count)) ./ 8
            for v in ps.v
                v .= T.(rand(rng, -16:16, count)) ./ 8
            end
            # One huge square forces complete arithmetic; pair cancellation,
            # folded a==b nodes and zero weights exercise replay ownership.
            ps.v[1][1] = huge
            if iseven(k)
                ps.v[1][2] = -huge
                ps.v[2][2] = ps.v[2][1]
                ps.v[3][2] = ps.v[3][1]
                ps.weight[2] = ps.weight[1]
                ps.x[1][2] = ps.x[1][1]
            end
            before = deepcopy((ps.x, ps.v, ps.weight, sh.s.H))
            exact = exact_nodes(sh, ps)
            actual = H._ion_energy_flux(sh, ps)
            @test all(isequal(actual[i], rounded(T, exact[i])) for i in eachindex(actual))
            @test isequal((ps.x, ps.v, ps.weight, sh.s.H), before)
        end
        for weight in (tiny, 3tiny, 5tiny),
            position in (T(0), T(0.5), T(7.5)),
            norm in (T(0.5), one(T), T(2))

            sh = PerpShock(8, T(7))
            sh.s.H .= norm
            ps = ParticleSet{1,T}(2)
            ps.x[1] .= position
            ps.weight .= weight
            ps.v[1] .= one(T)
            expected = exact_nodes(sh, ps)
            @test all(
                isequal(a, rounded(T, q)) for (a, q) in zip(H._ion_energy_flux(sh, ps), expected)
            )
        end
        for sign in (-1, 1),
            n in (tiny, 3tiny, 5tiny),
            u in (prevfloat(one(T)), one(T), nextfloat(one(T)))

            q = Q(n) * Q(sign * u)^3 / 2
            @test isequal(H._shock_bulk_energy_flux(n, sign * u, zero(T)), rounded(T, q))
            qh = Q(2) * Q(n) * Q(sign * T(0.25))
            @test isequal(H._shock_electron_enthalpy_flux(T(2), n, sign * T(0.25)), rounded(T, qh))
        end
    end
end
println("SHOCK_FLUX_PROPERTIES_COMPLETE")
flush(stdout)

end
