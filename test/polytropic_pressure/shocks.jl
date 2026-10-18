module ShockChecks
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function shock(::Type{T}, D, temperature, gamma; energy = false) where {T}
    parameters = (Te = temperature, γe = gamma, η = zero(T), τ = zero(T), B0 = zero(T))
    D == 1 && return PerpShock(3, one(T); parameters..., closure = energy ? :energy : :polytropic)
    D == 2 && return PerpShock2D(3, 4, one(T), one(T); parameters...)
    return PerpShock3D(3, 4, 4, one(T), one(T), one(T); parameters...)
end
function electric!(sh, D)
    D == 1 && return H.compute_E!(sh)
    D == 2 && return H.compute_E2d!(sh)
    return H.compute_E3d!(sh)
end
fields(sh, D) = D == 1 || D == 2 ? (sh.Ex, sh.Ey) : sh.E
@testset "Shock force owners preserve complete pointwise pressure range" begin
    for T in (Float32, Float64), D = 1:3, kind in (:overflow, :underflow)
        k = T === Float32 ? 100 : 600
        density = ldexp(one(T), kind === :overflow ? k : -k)
        temperature = inv(density)
        expected = H._from_exact(T, Q(temperature) * Q(density)^2)
        @test expected == density && isfinite(expected) && expected > 0
        sh = shock(T, D, temperature, T(2))
        fill!(sh.n, density)
        before = copy(sh.n)
        @test electric!(sh, D) === sh
        @test all(==(expected), sh.pe)
        @test sh.n == before
        @test all(a -> all(iszero, a), fields(sh, D))
    end
end
@testset "Shock initialization and energy boundaries retain complete pressure" begin
    for T in (Float32, Float64), kind in (:overflow, :underflow)
        k = T === Float32 ? 100 : 600
        density = ldexp(one(T), kind === :overflow ? k : -k)
        sh = shock(T, 1, inv(density), T(2); energy = true)
        particles = ParticleSet{1,T}(3)
        particles.x[1] .= sh.x
        particles.weight .= density .* sh.s.H
        @test H.init_shock!(sh, particles) === sh
        @test all(==(density), sh.n)
        @test all(==(density), sh.pe)
        @test H._rk4_pe_energy!(sh, T(0.001)) === sh
        @test first(sh.pe) == density && last(sh.pe) == density
        # The energy equation's interior pressure floor remains its own policy.
        kind === :overflow && @test sh.pe[2] == density
    end
end
@testset "Signed shock densities retain their distinct power contract" begin
    for T in (Float32, Float64), D = 1:3
        even = shock(T, D, T(0.7), T(2))
        fill!(even.n, -one(T))
        electric!(even, D)
        @test all(==(T(0.7)), even.pe)
        fractional = shock(T, D, T(0.7), T(1.5))
        fill!(fractional.n, -one(T))
        @test_throws DomainError electric!(fractional, D)
    end
end
end
