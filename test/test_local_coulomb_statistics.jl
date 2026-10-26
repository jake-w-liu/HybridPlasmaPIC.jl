module LocalCoulombStatisticsTests

using Test, Random, HybridPlasmaPIC, Statistics

# Frozen acceptance before first candidate measurements:
# 256 independent cells, 128 seeds, dt=1e-4; all expectation tests use seven
# measured standard errors plus 2.5% (drift) / 5% (diffusion) finite-dt allowance.
# Statistics are accumulated per cell, preserving within-cell correlations.
function cross_generator(::Type{T}, na, nb, weighted, reverse_order, cell_width) where {T}
    cells, trials = 256, 128
    grid = FourierGrid((cells,), (T(cells * cell_width),))
    a = ParticleSet{1,T}(cells * na; m = 1, q = -1)
    b = ParticleSet{1,T}(cells * nb; m = 4, q = 1)
    for (ps, count) in ((a, na), (b, nb)), cell = 1:cells, k = 1:count
        p = (cell - 1) * count + k
        ps.x[1][p] = (T(cell) - T(0.75)) * T(cell_width)
        ps.weight[p] = weighted ? T(2.0^(mod(k, 5) - 2)) : one(T)
    end
    species = reverse_order ? (b, a) : (a, b)
    workspace = CoulombCollisionWorkspace(species, grid)
    dt = T(1e-4)
    logarithm = T[0 8pi; 8pi 0]
    totals = zeros(Float64, 2)
    squares = zeros(Float64, 2)
    for trial = 1:trials
        fill!(a.v[1], one(T))
        fill!(b.v[1], -one(T))
        for ps in (a, b), c = 2:3
            fill!(ps.v[c], zero(T))
        end
        collide_coulomb!(
            species,
            grid,
            dt;
            coulomb_log = logarithm,
            epsilon0 = 1,
            weight_to_number = 1,
            workspace,
            rng = MersenneTwister(0x521 + trial),
        )
        for cell = 1:cells
            dx, dy2 = 0.0, 0.0
            for k = 1:na
                p = (cell - 1) * na + k
                dx += Float64(a.v[1][p]) - 1
                dy2 += Float64(a.v[2][p])^2
            end
            values = (dx / na, dy2 / na)
            for i = 1:2
                totals[i] += values[i]
                squares[i] += values[i]^2
            end
        end
    end
    mu = Float64(a.m) * Float64(b.m) / (Float64(a.m) + Float64(b.m))
    K = (Float64(a.q) * Float64(b.q))^2 * Float64(logarithm[1, 2]) / (8pi * mu^2)
    fraction = Float64(b.m) / (Float64(a.m) + Float64(b.m))
    wb = sum(Float64, b.weight[1:nb]) / Float64(grid.dx[1])
    n = cells * trials
    expected = (-fraction * K * wb * Float64(dt) / 2, fraction^2 * K * wb * Float64(dt))
    for i = 1:2
        mean = totals[i] / n
        se = sqrt(max(0, squares[i] / n - mean^2) / (n - 1))
        allowance = 7se + (i == 1 ? 0.025 : 0.05) * abs(expected[i])
        @info "cross Landau generator" T na nb weighted reverse_order cell_width metric = i mean expected =
            expected[i] se allowance
        @test abs(mean - expected[i]) <= allowance
    end
end
function self_generator(::Type{T}, count, weighted) where {T}
    cells, trials = 256, 128
    grid = FourierGrid((cells,), (T(cells),))
    ps = ParticleSet{1,T}(cells * count; m = 1, q = 1)
    positive = count ÷ 2
    for cell = 1:cells, k = 1:count
        p = (cell - 1) * count + k
        ps.x[1][p] = T(cell) - T(0.75)
        ps.weight[p] = weighted ? T(2.0^(mod(k, 5) - 2)) : one(T)
    end
    workspace = CoulombCollisionWorkspace((ps,), grid)
    dt = T(1e-4)
    logarithm = T(8pi)
    totals = zeros(2)
    squares = zeros(2)
    for trial = 1:trials
        for cell = 1:cells, k = 1:count
            ps.v[1][(cell-1)*count+k] = k <= positive ? one(T) : -one(T)
        end
        fill!(ps.v[2], zero(T))
        fill!(ps.v[3], zero(T))
        collide_coulomb!(
            (ps,),
            grid,
            dt;
            coulomb_log = logarithm,
            epsilon0 = 1,
            weight_to_number = 1,
            workspace,
            rng = MersenneTwister(0x681 + trial),
        )
        for cell = 1:cells
            dx, dy2 = 0.0, 0.0
            for k = 1:positive
                p = (cell - 1) * count + k
                dx += Float64(ps.v[1][p]) - 1
                dy2 += Float64(ps.v[2][p])^2
            end
            values = (dx / positive, dy2 / positive)
            for i = 1:2
                totals[i] += values[i]
                squares[i] += values[i]^2
            end
        end
    end
    K = Float64(logarithm) / (8pi * 0.5^2)
    wb = sum(Float64, ps.weight[(positive+1):count])
    n = cells * trials
    expected = (-K * wb * Float64(dt) / 4, K * wb * Float64(dt) / 4)
    for i = 1:2
        mean = totals[i] / n
        se = sqrt(max(0, squares[i] / n - mean^2) / (n - 1))
        allowance = 7se + (i == 1 ? 0.025 : 0.05) * abs(expected[i])
        @info "self Landau generator" T count weighted metric = i mean expected = expected[i] se allowance
        @test abs(mean - expected[i]) <= allowance
    end
end

# A coupled four-species ring checks that the public operation sums multiple
# enabled physical channels and honors zero matrix entries. The oracle sums
# the independent first-order Landau generators at the initial velocities;
# it does not use the production pair scheduler or scattering implementation.
function multispecies_generator(::Type{T}, reverse_order) where {T}
    cells, trials = 256, 128
    grid = FourierGrid((cells,), (T(cells),))
    masses, charges = (1, 2, 4, 8), (1, -1, 2, 1)
    weights, initial = (0.5, 2.0, 1.0, 4.0), (-3.0, -1.0, 1.0, 4.0)
    populations = ntuple(s -> ParticleSet{1,T}(cells; m = masses[s], q = charges[s]), 4)
    logarithms = zeros(T, 4, 4)
    for (a, b, scale) in ((1, 2, 1.0), (2, 3, 0.5), (3, 4, 2.0), (1, 4, 0.25))
        logarithms[a, b] = logarithms[b, a] = T(8pi * scale)
    end
    for s = 1:4
        fill!(populations[s].weight, T(weights[s]))
        for cell = 1:cells
            populations[s].x[1][cell] = T(cell) - T(0.75)
        end
    end
    order = reverse_order ? (4, 3, 2, 1) : (1, 2, 3, 4)
    species = map(s -> populations[s], order)
    channels = logarithms[collect(order), collect(order)]
    workspace = CoulombCollisionWorkspace(species, grid)
    dt = T(1e-4)
    expected = zeros(Float64, 4, 2)
    for a = 1:4, b = 1:4
        iszero(logarithms[a, b]) && continue
        mu = masses[a] * masses[b] / (masses[a] + masses[b])
        K = (charges[a] * charges[b])^2 * Float64(logarithms[a, b]) / (8pi * mu^2)
        fraction = masses[b] / (masses[a] + masses[b])
        relative = initial[a] - initial[b]
        density = weights[b] / Float64(grid.dx[1])
        expected[a, 1] -= 2fraction * K * density * Float64(dt) * relative / abs(relative)^3
        expected[a, 2] += 2fraction^2 * K * density * Float64(dt) / abs(relative)
    end
    totals, squares = zeros(4, 2), zeros(4, 2)
    for trial = 1:trials
        for s = 1:4
            fill!(populations[s].v[1], T(initial[s]))
            fill!(populations[s].v[2], zero(T))
            fill!(populations[s].v[3], zero(T))
        end
        collide_coulomb!(
            species,
            grid,
            dt;
            coulomb_log = channels,
            epsilon0 = 1,
            weight_to_number = 1,
            workspace,
            rng = MersenneTwister(0x781 + trial),
        )
        for s = 1:4, cell = 1:cells
            values = (
                Float64(populations[s].v[1][cell]) - initial[s],
                Float64(populations[s].v[2][cell])^2,
            )
            for metric = 1:2
                totals[s, metric] += values[metric]
                squares[s, metric] += values[metric]^2
            end
        end
    end
    n = cells * trials
    for s = 1:4, metric = 1:2
        mean = totals[s, metric] / n
        se = sqrt(max(0, squares[s, metric] / n - mean^2) / (n - 1))
        allowance = 7se + (metric == 1 ? 0.025 : 0.05) * abs(expected[s, metric])
        @info "coupled multispecies Landau generator" T reverse_order species = s metric mean expected =
            expected[s, metric] se allowance
        @test abs(mean - expected[s, metric]) <= allowance
    end
end

@testset "Independent empirical Landau generator statistics" begin
    for T in (Float32, Float64),
        (na, nb) in ((1, 1), (2, 3), (3, 7), (7, 3)),
        weighted in (false, true),
        reverse_order in (false, true),
        cell_width in (1.0, 0.3)

        cross_generator(T, na, nb, weighted, reverse_order, cell_width)
    end
    for T in (Float32, Float64), n in (2, 3, 4, 7), weighted in (false, true)
        self_generator(T, n, weighted)
    end
    for T in (Float32, Float64), reverse_order in (false, true)
        multispecies_generator(T, reverse_order)
    end
end

end
