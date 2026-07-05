# Float-edge regressions for the particle boundaries and the cell_index ↔
# deposition convention: the half-open [lo,hi) contract at the periodic seam,
# non-finite rejection at the open (absorbing) boundary, and cell/deposit
# agreement for positions on or just outside the box edges.

using HybridPlasmaPIC, Test

@testset "boundary float edges" begin

    @testset "apply_periodic! keeps the seam inside [lo,hi)" begin
        # mod(-1e-20, 1.0) rounds to exactly 1.0: the wrap must fold hi back to lo
        ps = ParticleSet{1,Float64}(2)
        ps.x[1] .= [-1e-20, 0.5]
        apply_periodic!(ps, (0.0,), (1.0,))
        @test all(0.0 .<= ps.x[1] .< 1.0)
        @test ps.x[1][1] == 0.0              # was exactly 1.0 == hi before the fix
        @test ps.x[1][2] == 0.5              # interior particle untouched

        # nonzero lo: l + mod(x−l, L) can round up to hi even when mod < L
        ps = ParticleSet{1,Float64}(1)
        ps.x[1][1] = prevfloat(1.0)
        apply_periodic!(ps, (1.0,), (1.1,))
        @test 1.0 <= ps.x[1][1] < 1.1        # was exactly 1.1 == hi before the fix
        @test ps.x[1][1] == 1.0
    end

    @testset "apply_absorbing! rejects non-finite positions" begin
        for bad in (NaN, Inf, -Inf)
            ps = ParticleSet{1,Float64}(2)
            ps.x[1] .= [bad, 0.5]
            @test_throws ArgumentError apply_absorbing!(ps, (0.0,), (1.0,))
            @test nparticles(ps) == 2        # nothing silently removed (was nremoved == 1)
        end
        # non-finite on a later axis is caught too
        ps2 = ParticleSet{2,Float64}(1)
        ps2.x[1] .= 0.5
        ps2.x[2] .= NaN
        @test_throws ArgumentError apply_absorbing!(ps2, (0.0, 0.0), (1.0, 1.0))
        # finite out-of-box particles are still the sanctioned loss
        ps = ParticleSet{1,Float64}(3)
        ps.x[1] .= [-0.1, 0.5, 1.0]
        @test apply_absorbing!(ps, (0.0,), (1.0,)) == 2
        @test ps.x[1] == [0.5]
    end

    @testset "cell_index matches the deposition cell at the box edges" begin
        n = 16
        g = FourierGrid((n,), (1.0,))
        # positions with an integral stencil coordinate s = mod(x,L)/dx: CIC puts
        # the full weight on the base node, so argmax(ρ) IS the deposit's home cell
        cases = (
            (0.0, 1),
            (0.5, 9),
            (1.0, 1),                        # x == L is the periodic image of 0 (was cell 16)
            (-1e-20, 1),                     # mod rounds to L: deposition wraps to node 1
            (-0.25, 13),                     # out-of-box folds to 0.75 (clamp gave cell 1)
        )
        for (x, expected) in cases
            ps = ParticleSet{1,Float64}(1)
            ps.x[1][1] = x
            ps.weight[1] = 1.0
            ρ = zeros(n)
            deposit_scalar!(ρ, ps, ps.weight, g, CIC())
            ci = cell_index(ps, g)[1]
            @test ci == argmax(ρ)            # the invariant: same cell as the deposit
            @test ci == expected
        end
        # x = prevfloat(L): truly inside the last cell [15/16, 1), so cell 16; the
        # CIC stencil anchors at the same base cell (weight split over its two
        # nodes, node 16 and the wrapped node 1 — nowhere else)
        ps = ParticleSet{1,Float64}(1)
        ps.x[1][1] = prevfloat(1.0)
        ps.weight[1] = 1.0
        ρ = zeros(n)
        deposit_scalar!(ρ, ps, ps.weight, g, CIC())
        @test cell_index(ps, g)[1] == 16
        @test ρ[16] + ρ[1] ≈ 1.0
        @test all(iszero, ρ[2:15])
    end

    @testset "sort/histogram smoke across the wrapped seam" begin
        g = FourierGrid((16,), (1.0,))
        ps = ParticleSet{1,Float64}(3)
        ps.x[1] .= [-1e-20, 0.97, 0.5]
        ps.id .= UInt64[1, 2, 3]
        apply_periodic!(ps, (0.0,), (1.0,))
        ppc = particles_per_cell(ps, g)
        @test sum(ppc) == 3
        @test ppc[1] == 1                    # seam particle counted with its deposit cell
        @test ppc[16] == 1 && ppc[9] == 1
        sort_particles!(ps, g)
        @test issorted(cell_index(ps, g))
        @test sort(ps.id) == UInt64[1, 2, 3]
    end
end
