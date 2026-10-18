module ShockExactReuseResourcesTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
state_cost(sh::S, ps::P) where {S,P} = @timed H._shock_ion_flux_state(sh, ps)
public_cost(sh::S, ps::P) where {S,P} = @timed boundary_energy_flux(sh; ps = ps)
@testset "Complete totals reuse already accumulated exact boundary nodes" begin
    for T in (Float32, Float64),
        count in (256, 4096),
        profile in (:overflow, :subnormal, :normalization)

        sh = PerpShock(8, T(7); γe = T(2))
        ps = ParticleSet{1,T}(count)
        ps.x[1] .= T[isodd(i) ? 0 : 7 for i = 1:count]
        sh.ux .= T(2)
        sh.pe .= T(0.25)
        if profile === :overflow
            n = ldexp(one(T), exponent(floatmax(T)) - 1)
            sh.n .= n
            sh.Bz .= T(4)
            sh.Ey .= -n
            ps.v[1] .= T(2)
            ps.weight .= n / count
            expected = one(T)
        elseif profile === :subnormal
            q = nextfloat(zero(T))
            sh.ux .= T(0.25)
            sh.pe .= q
            sh.Ey .= q
            sh.Bz .= T(0.5)
            ps.weight .= 32q
            ps.v[1] .= T(0.25)
            expected = (T(count) / 4 + 1) * q
        else
            ve, he = T === Float32 ? (50, 100) : (400, 800)
            sh.s.H .= ldexp(one(T), he)
            ps.v[1] .= ldexp(one(T), ve)
            sh.Bz .= one(T)
            sh.Ey .= -ldexp(T(count) / 2, 3ve - 1 - he)
            expected = one(T)
        end
        state_cost(sh, ps)
        public_cost(sh, ps)
        states = [state_cost(sh, ps) for _ = 1:3]
        totals = [public_cost(sh, ps) for _ = 1:3]
        for (state, total) in zip(states, totals)
            @test state.value.exact !== nothing
            @test all(haskey(state.value.exact, i) for i in (1, 8))
            @test total.value.total == (expected, expected)
        end
        # The enclosing total may add bounded field arithmetic, but replaying
        # both boundary populations nearly doubles the measured state work.
        # A 10% allowance retains that distinction at both tested sizes.
        @test maximum(s.bytes for s in totals) <= 1.1minimum(s.bytes for s in states)
    end
end
end
