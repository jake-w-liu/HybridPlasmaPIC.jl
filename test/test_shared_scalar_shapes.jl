using Test
using HybridPlasmaPIC

# A public custom shape may return weights wider than the particle/grid type.
# Dyadic coordinates keep the independent CIC moments exactly representable.
struct _SharedWideCIC <: ShapeFunction end
HybridPlasmaPIC.width(::_SharedWideCIC) = 2
function HybridPlasmaPIC._stencil1d(::_SharedWideCIC, s)
    base, weights = HybridPlasmaPIC._stencil1d(CIC(), s)
    return base, map(Float64, weights)
end

@testset "Shared scalar mixed-precision shape contract" begin
    grid = FourierGrid((4,), (4.0f0,))
    particle = ParticleSet{1,Float32}(1; q = 2.0f0)
    particle.x[1][1] = 0.25f0
    particle.v[1][1], particle.v[2][1] = 2.0f0, -4.0f0
    moments = ntuple(_ -> zeros(Float32, 4), 3)
    momentum!(moments, particle, grid, _SharedWideCIC())
    @test moments == (Float32[1.5, 0.5, 0, 0], Float32[-3, -1, 0, 0], zeros(Float32, 4))
    fields = HybridFields{1,Float32}(grid.n)
    compute_moments_multi!(fields, [particle], grid, _SharedWideCIC(), 0.125f0)
    @test fields.n == Float32[1.5, 0.5, 0, 0]
    @test fields.ui == (Float32[2, 2, 0, 0], Float32[-4, -4, 0, 0], zeros(Float32, 4))
end
