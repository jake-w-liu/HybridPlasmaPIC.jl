using HybridPlasmaPIC, Test

function _weighted_particles(weights, vx, vy, vz)
    ps = ParticleSet{1,Float64}(length(weights))
    ps.weight .= weights
    ps.v[1] .= vx
    ps.v[2] .= vy
    ps.v[3] .= vz
    return ps
end

@testset "diagnostic scalar and array reductions" begin
    @test sum_diagnostics([1, 2, 3]) == 6
    @test min_diagnostics([3.0, -2.0, 5.0]) == -2.0
    @test max_diagnostics([3.0, -2.0, 5.0]) == 5.0

    arrays = [[1.0, 2.0, 3.0], [0.5, -1.0, 4.0], [2.0, 3.0, -2.0]]
    @test sum_diagnostics(arrays) == [3.5, 4.0, 5.0]
    @test min_diagnostics(arrays) == [0.5, -1.0, -2.0]
    @test max_diagnostics(arrays) == [2.0, 3.0, 4.0]

    @test_throws ArgumentError sum_diagnostics(Any[])
    @test_throws ArgumentError reduce_diagnostics([1, 2]; op = :mean)
    @test_throws DimensionMismatch sum_diagnostics([[1.0, 2.0], [1.0 2.0]])
    @test_throws ArgumentError sum_diagnostics(Any[1.0, "bad"])
end

@testset "nested diagnostic reductions preserve structure" begin
    locals = [
        (
            energy = (kinetic = 1.0, magnetic = 2.0, total = 3.0),
            momentum = (1.0, 2.0, 3.0),
            hist = [1, 0, 2],
        ),
        (
            energy = (kinetic = 4.0, magnetic = 5.0, total = 9.0),
            momentum = (-1.0, 1.0, 0.0),
            hist = [0, 3, 4],
        ),
        (
            energy = (kinetic = 2.5, magnetic = 1.5, total = 4.0),
            momentum = (0.5, -3.0, 1.0),
            hist = [2, 2, 0],
        ),
    ]

    reduced = sum_diagnostics(locals)
    @test keys(reduced) == (:energy, :momentum, :hist)
    @test reduced.energy == (kinetic = 7.5, magnetic = 8.5, total = 16.0)
    @test reduced.momentum == (0.5, 0.0, 4.0)
    @test reduced.hist == [3, 5, 6]

    @test_throws ArgumentError sum_diagnostics([
        (energy = 1.0, momentum = 2.0),
        (energy = 3.0, current = 4.0),
    ])
    @test_throws ArgumentError sum_diagnostics([(1.0, 2.0), (3.0,)])
end

@testset "split particle diagnostics reduce to serial values" begin
    weights = [1.0, 0.5, 2.0, 1.5, 0.25]
    vx = [1.0, -2.0, 0.5, 3.0, -4.0]
    vy = [0.0, 1.0, -1.5, 2.5, 3.0]
    vz = [2.0, 0.0, -1.0, 1.0, -0.5]

    serial = _weighted_particles(weights, vx, vy, vz)
    left = _weighted_particles(weights[1:2], vx[1:2], vy[1:2], vz[1:2])
    mid = _weighted_particles(weights[3:4], vx[3:4], vy[3:4], vz[3:4])
    right = _weighted_particles(weights[5:5], vx[5:5], vy[5:5], vz[5:5])

    local_diagnostics = [
        (number = sum(left.weight), momentum = total_momentum(left)),
        (number = sum(mid.weight), momentum = total_momentum(mid)),
        (number = sum(right.weight), momentum = total_momentum(right)),
    ]
    reduced = sum_diagnostics(local_diagnostics)

    @test reduced.number ≈ sum(serial.weight)
    @test reduced.momentum[1] ≈ total_momentum(serial)[1]
    @test reduced.momentum[2] ≈ total_momentum(serial)[2]
    @test reduced.momentum[3] ≈ total_momentum(serial)[3]
end
