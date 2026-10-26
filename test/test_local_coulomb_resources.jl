module LocalCoulombResourceTests
using Test, Random, HybridPlasmaPIC
@noinline function call!(species::S, g::G, workspace::W, rng::R) where {S,G,W,R}
    collide_coulomb!(
        species,
        g,
        0.001;
        coulomb_log = 1,
        epsilon0 = 1,
        weight_to_number = 1,
        workspace,
        rng,
    )
end
@noinline function allocation_call(species::S, g::G, w::W, r::R) where {S,G,W,R}
    return @allocated call!(species, g, w, r)
end
function profile(::Type{T}, D, n, species_count, collection) where {T}
    grid = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> T(8), D))
    parts = ntuple(s -> ParticleSet{D,T}(n; m = T(s), q = T(s)), species_count)
    for (s, ps) in enumerate(parts)
        for x in ps.x
            fill!(x, T(0.25))
        end
        for c = 1:3, p = 1:n
            ps.v[c][p] = T(s + c + mod(p, 11)) / 4
        end
    end
    species = collection == :tuple ? parts : collect(parts)
    workspace = CoulombCollisionWorkspace(species, grid)
    rng = MersenneTwister(5)
    for _ = 1:3
        call!(species, grid, workspace, rng)
    end
    samples = ntuple(_ -> allocation_call(species, grid, workspace, rng), 3)
    @info "local allocation profile" T D n species_count collection samples
    @test all(iszero, samples)
end
@testset "Warmed reusable local collision allocations" begin
    for T in (Float32, Float64),
        D = 1:3,
        n in (0, 16, 4096),
        s in (1, 2, 4),
        collection in (:tuple, :vector)

        profile(T, D, n, s, collection)
    end
end

end
