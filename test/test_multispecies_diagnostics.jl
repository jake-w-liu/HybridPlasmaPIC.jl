using HybridPlasmaPIC, Test

function diagnostic_populations(::Type{T}, ::Val{D}) where {T,D}
    species = [ParticleSet{D,T}(s == 3 ? 0 : 3s; q = T(s), m = T(s^2)) for s = 1:3]
    for (s, ps) in enumerate(species)
        for p in eachindex(ps.weight)
            ps.weight[p] = T(p + s) / T(7)
            for d = 1:D
                ps.x[d][p] = T(p) / T(10)
            end
            for c = 1:3
                ps.v[c][p] = T((-1)^(p + c) * (p + c + s)) / T(13)
            end
        end
    end
    return species
end

# Independent 512-bit evaluation of the discrete physical sums. The field
# is B=(0.8, 0.125sin(kx), 0), so J=(0,0,0.125k cos(kx)).
function diagnostic_oracle(species, B, density, closure, g, de2)
    return setprecision(512) do
        energy, momentum = BigFloat(0), zeros(BigFloat, 3)
        for ps in species, p in eachindex(ps.weight)
            mass_weight = BigFloat(ps.m) * BigFloat(ps.weight[p])
            energy += mass_weight * sum(BigFloat(ps.v[c][p])^2 for c = 1:3) / 2
            for c = 1:3
                momentum[c] += mass_weight * BigFloat(ps.v[c][p])
            end
        end
        magnetic, electron, inertia = BigFloat(0), BigFloat(0), BigFloat(0)
        k = 2BigFloat(π) / BigFloat(g.L[1])
        for I in CartesianIndices(density)
            b2 = sum(BigFloat(B[c][I])^2 for c = 1:3)
            n = BigFloat(density[I])
            magnetic += b2 / 2
            if closure isa IsothermalElectrons
                electron += BigFloat(closure.Te) * n * log(n)
            elseif closure isa PolytropicElectrons
                gamma = BigFloat(closure.γ)
                electron += BigFloat(closure.pe0) * (n / BigFloat(closure.n0))^gamma / (gamma - 1)
            else
                electron +=
                    BigFloat(closure.Cperp) * n * sqrt(b2) + BigFloat(closure.Cpar) * n^3 / (2b2)
            end
            phase = 2BigFloat(π) * (I[1] - 1) / g.n[1]
            inertia += BigFloat(de2) * (BigFloat(0.125) * k * cos(phase))^2 / 2
        end
        volume = prod(BigFloat.(g.dx))
        magnetic *= volume
        electron *= volume
        inertia *= volume
        return (;
            energy,
            momentum = Tuple(momentum),
            magnetic,
            electron,
            inertia,
            total = energy + magnetic + electron + inertia,
        )
    end
end

@testset "Multispecies energy and momentum budgets" begin
    for T in (Float32, Float64), D = 1:3
        species = diagnostic_populations(T, Val(D))
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> T(2π), D))
        B = ntuple(_ -> zeros(T, g.n), 3)
        fill!(B[1], T(0.8))
        for I in CartesianIndices(g.n)
            B[2][I] = T(0.125) * sin(T(2π) * T(I[1] - 1) / T(g.n[1]))
        end
        density = fill(T(1.4), g.n)
        closures = (
            IsothermalElectrons(T(0.2)),
            PolytropicElectrons(T(0.2), one(T), T(5 / 3)),
            CGLElectrons(T(0.2), T(0.3), one(T), one(T)),
        )
        for closure in closures, de2 in (zero(T), T(0.2))
            expected = diagnostic_oracle(species, B, density, closure, g, de2)
            actual = energy_budget(species, B, density, closure, g; de2)
            momentum = momentum_budget(species, B, g)
            # Field reductions span up to 512 cells; allow accumulated roundoff
            # in each production precision, far below an omitted population.
            for (value, reference) in zip(
                values(actual),
                (
                    expected.energy,
                    expected.magnetic,
                    expected.electron,
                    expected.inertia,
                    expected.total,
                ),
            )
                @test isapprox(value, reference; rtol = 128eps(T), atol = 128eps(T))
                @test value isa T
            end
            @test isapprox(kinetic_energy(species), expected.energy; rtol = 16eps(T))
            @test all(
                isapprox(a, b; rtol = 16eps(T), atol = 16eps(T)) for
                (a, b) in zip(total_momentum(species), expected.momentum)
            )
            @test momentum.particle == momentum.total == total_momentum(species)
            @test all(x -> x isa T, momentum.total)
            # Shared field/electron reservoirs are counted exactly once.
            one_population = energy_budget([species[1]], B, density, closure, g; de2)
            scalar = energy_budget(species[1], B, density, closure, g; de2)
            @test one_population == scalar
            @test momentum_budget([species[1]], B, g) == momentum_budget(species[1], B, g)
            for field in (:magnetic, :electron_internal, :electron_inertia)
                @test getproperty(actual, field) == getproperty(scalar, field)
            end

            none = ParticleSet{D,T}[]
            empty_budget = energy_budget(none, B, density, closure, g; de2)
            @test kinetic_energy(none) === zero(T)
            @test total_momentum(none) === (zero(T), zero(T), zero(T))
            @test empty_budget.kinetic === zero(T)
            @test empty_budget.total ==
                  actual.magnetic + actual.electron_internal + actual.electron_inertia
            @test momentum_budget(none, B, g).total === (zero(T), zero(T), zero(T))
            @test energy_budget([species[3]], B, density, closure, g; de2) == empty_budget
        end

        viewed = map(species) do ps
            ParticleSet{D,T}(
                map(a -> view(a, :), ps.x),
                map(a -> view(a, :), ps.v),
                view(ps.weight, :),
                view(ps.id, :),
                view(ps.tag, :),
                ps.q,
                ps.m,
            )
        end
        @test kinetic_energy(viewed) == kinetic_energy(species)
        @test total_momentum(viewed) == total_momentum(species)
        @test kinetic_energy(view(species, 1:2)) == kinetic_energy(species)
        @test energy_budget(viewed, B, density, closures[1], g) ==
              energy_budget(species, B, density, closures[1], g)
        for invalid in (-one(T), T(Inf), T(NaN))
            @test_throws ArgumentError energy_budget(
                species,
                B,
                density,
                closures[1],
                g;
                de2 = invalid,
            )
        end
        bad_B = ntuple(_ -> zeros(T, ntuple(_ -> 3, D)), 3)
        @test_throws DimensionMismatch energy_budget(species, bad_B, density, closures[1], g)
        @test_throws DimensionMismatch momentum_budget(species, bad_B, g)
        @test_throws DimensionMismatch energy_budget(
            species,
            B,
            zeros(T, ntuple(_ -> 3, D)),
            closures[1],
            g,
        )
    end
end

@testset "Heterogeneous nonempty momentum collections remain supported" begin
    a = diagnostic_populations(Float32, Val(1))[1]
    b = diagnostic_populations(Float64, Val(2))[2]
    population = ParticleSet[a, b]
    expected = total_momentum(a) .+ total_momentum(b)
    @test total_momentum(population) == expected
    @test total_momentum(view(population, 2:-1:1)) == expected
    @test_throws ArgumentError total_momentum(ParticleSet[])
end
