module LocalCoulombThermalTests

using Test, Random, HybridPlasmaPIC, Statistics

mutable struct AntitheticNormal <: AbstractRNG
    inner::MersenneTwister
    sign::Float64
end
Random.rand(r::AntitheticNormal, range::UnitRange{Int}) = rand(r.inner, range)
Random.rand(r::AntitheticNormal, T::Union{Type{Float32},Type{Float64}}) = rand(r.inner, T)
Random.randn(r::AntitheticNormal, T::Union{Type{Float32},Type{Float64}}) =
    T(r.sign) * randn(r.inner, T)
function make_species(::Type{T}, n) where {T}
    species = (ParticleSet{1,T}(n; m = 1), ParticleSet{1,T}(n; m = 4))
    for ps in species
        fill!(ps.x[1], T(0.25))
        fill!(ps.weight, inv(T(n)))
    end
    return species
end
energy(ps) = sum(
    p -> Float64(ps.m) * Float64(ps.weight[p]) * sum(c -> Float64(ps.v[c][p])^2, 1:3) / 2,
    1:nparticles(ps),
)
function thermal_rate(::Type{T}, dt0) where {T}
    n, trials = 512, 256
    grid = FourierGrid((1,), (one(T),))
    plus, minus = make_species(T, n), make_species(T, n)
    pw, mw = CoulombCollisionWorkspace(plus, grid), CoulombCollisionWorkspace(minus, grid)
    samples = zeros(trials)
    dt = T(dt0)
    logarithm = T(8pi)
    for trial = 1:trials
        rng = MersenneTwister(0x934 + trial)
        for (s, temperature) in enumerate((0.25, 2.0)), c = 1:3
            randn!(rng, plus[s].v[c])
            plus[s].v[c] .*= T(sqrt(temperature / Float64(plus[s].m)))
            copyto!(minus[s].v[c], plus[s].v[c])
        end
        initial = energy(plus[1])
        collide_coulomb!(
            plus,
            grid,
            dt;
            coulomb_log = T[0 logarithm; logarithm 0],
            epsilon0 = 1,
            weight_to_number = 1,
            u_floor = T(0.001),
            workspace = pw,
            rng = AntitheticNormal(MersenneTwister(0x714 + trial), 1.0),
        )
        collide_coulomb!(
            minus,
            grid,
            dt;
            coulomb_log = T[0 logarithm; logarithm 0],
            epsilon0 = 1,
            weight_to_number = 1,
            u_floor = T(0.001),
            workspace = mw,
            rng = AntitheticNormal(MersenneTwister(0x714 + trial), -1.0),
        )
        samples[trial] = ((energy(plus[1]) + energy(minus[1])) / 2 - initial) / Float64(dt)
    end
    # Maxwellian Landau energy-exchange integral: E[1/|g|]=sqrt(2/pi)/sigma,
    # E[va|g]=(Ta/ma)/sigma²*g. This is independent of the pairing code.
    ma, mb, Ta, Tb = 1.0, 4.0, 0.25, 2.0
    reduced = ma * mb / (ma + mb)
    sigma = sqrt(Ta / ma + Tb / mb)
    K = Float64(logarithm) / (8pi * reduced^2)
    expected = 2 * reduced * K / (ma + mb) * sqrt(2 / pi) * (Tb - Ta) / sigma^3
    measured = mean(samples)
    se = std(samples) / sqrt(trials)
    @info "Maxwellian thermal exchange" T dt measured expected se
    return measured, expected, se
end
@testset "Independent Maxwellian energy exchange and timestep refinement" begin
    for T in (Float32, Float64)
        coarse = thermal_rate(T, 0.02)
        fine = thermal_rate(T, 0.00125)
        finest = thermal_rate(T, 0.0003125)
        @test abs(finest[1] - finest[2]) <= 7finest[3] + 0.03abs(finest[2])
        @test abs(fine[1] - fine[2]) <= 0.6abs(coarse[1] - coarse[2]) + 7 * (fine[3] + coarse[3])
        @test abs(finest[1] - finest[2]) <= abs(fine[1] - fine[2]) + 7 * (fine[3] + finest[3])
    end
end
function thermal_state(species)
    density = sum(ps -> sum(ps.weight), species)
    massdensity = sum(ps -> ps.m * sum(ps.weight), species)
    bulk = ntuple(c -> sum(ps -> ps.m * sum(ps.weight .* ps.v[c]), species) / massdensity, 3)
    temperatures = map(
        ps ->
            ps.m *
            sum(p -> ps.weight[p] * sum(c -> (ps.v[c][p] - bulk[c])^2, 1:3), 1:nparticles(ps)) /
            (3sum(ps.weight)),
        species,
    )
    return temperatures, bulk, density
end
@testset "Unequal physical masses reach common temperature" begin
    for seed in (0x14, 0x91)
        T = Float64
        n = 1024
        grid = FourierGrid((1,), (1.0,))
        species = make_species(T, n)
        rng = MersenneTwister(seed)
        for (s, temperature) in enumerate((0.25, 2.0)), c = 1:3
            randn!(rng, species[s].v[c])
            species[s].v[c] .*= sqrt(temperature / species[s].m)
        end
        before = thermal_state(species)
        target = sum(before[1]) / 2
        initialenergy = sum(energy, species)
        workspace = CoulombCollisionWorkspace(species, grid)
        for step = 1:2000
            collide_coulomb!(
                species,
                grid,
                0.01;
                coulomb_log = 8pi,
                epsilon0 = 1,
                weight_to_number = 1,
                workspace,
                rng,
            )
        end
        after = thermal_state(species)
        @info "common physical temperature" seed initial = before[1] target measured = after[1]
        @test abs(after[1][1] - after[1][2]) <= 0.15target
        @test all(value -> abs(value - target) <= 0.1target, after[1])
        @test sum(energy, species) ≈ initialenergy rtol = 1e-10
        @test all(isapprox.(after[2], before[2]; atol = 1e-11, rtol = 1e-10))
    end
end

end
