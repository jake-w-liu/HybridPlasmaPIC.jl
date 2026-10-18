module ShockFluxResourcesTests
using HybridPlasmaPIC, Test, Random
include("legacy.jl")

const H = HybridPlasmaPIC
measure_old(sh::S, ps::P) where {S,P} = @timed Legacy._ion_energy_flux(sh, ps)
measure_new(sh::S, ps::P) where {S,P} = @timed H._ion_energy_flux(sh, ps)
println((VERSION,))
flush(stdout)
@testset "Larger ordinary shock flux compatibility and resources" begin
    rng = MersenneTwister(183992)
    for T in (Float32, Float64),
        (count, nodes) in ((512, 64), (8192, 64), (8192, 512)),
        profile in (:positive, :signed)

        sh = PerpShock(nodes, T(nodes - 1))
        ps = ParticleSet{1,T}(count)
        ps.x[1] .= T.(rand(rng, count) .* (nodes - 1))
        ps.weight .= T.(rand(rng, count))
        for v in ps.v
            v .= T.(randn(rng, count))
        end
        profile === :positive && (ps.v[1] .= abs.(ps.v[1]))
        measure_old(sh, ps)
        measure_new(sh, ps)
        old = measure_old(sh, ps)
        new = measure_new(sh, ps)
        mismatches = Base.count(!isequal(a, b) for (a, b) in zip(old.value, new.value))
        println((
            T,
            particles = length(ps.weight),
            nodes,
            profile,
            old_bytes = old.bytes,
            new_bytes = new.bytes,
            old_seconds = old.time,
            new_seconds = new.time,
            mismatches,
        ))
        flush(stdout)
        @test isequal(old.value, new.value)
        @test new.bytes == old.bytes
    end
end

end
