module LocalCoulombTests

using HybridPlasmaPIC
using Test
using Random
using Serialization

const H = HybridPlasmaPIC

function population(::Type{T}, D, count; q = 1, m = 1, seed = 0x87) where {T}
    ps = ParticleSet{D,T}(count; q, m)
    for x in ps.x
        fill!(x, T(0.25))
    end
    random = MersenneTwister(seed)
    for v in ps.v
        randn!(random, v)
    end
    ps.tag .= UInt32.(1:count)
    return ps
end

state(ps) =
    (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m)

function conserved(species)
    momentum = zeros(BigFloat, 3)
    energy = BigFloat(0)
    for ps in species, p = 1:nparticles(ps)
        weightmass = BigFloat(ps.weight[p]) * BigFloat(ps.m)
        for c = 1:3
            velocity = BigFloat(ps.v[c][p])
            momentum[c] += weightmass * velocity
            energy += weightmass * velocity^2 / 2
        end
    end
    return momentum, energy
end

mutable struct FixedRNG <: AbstractRNG
    normal::Float64
    uniform::Vector{Float64}
    next::Int
end
FixedRNG(normal = 1.0, uniform = [0.125]) = FixedRNG(normal, uniform, 1)
Random.rand(::FixedRNG, r::UnitRange{Int}) = last(r)
Random.randn(r::FixedRNG, T::Union{Type{Float32},Type{Float64}}) = T(r.normal)
function Random.rand(r::FixedRNG, T::Union{Type{Float32},Type{Float64}})
    value = r.uniform[mod1(r.next, length(r.uniform))]
    r.next += 1
    return T(value)
end

mutable struct LateInvalidRNG <: AbstractRNG
    inner::MersenneTwister
    count::Int
end
Random.rand(r::LateInvalidRNG, range::UnitRange{Int}) = rand(r.inner, range)
Random.rand(r::LateInvalidRNG, T::Union{Type{Float32},Type{Float64}}) = rand(r.inner, T)
function Random.randn(r::LateInvalidRNG, T::Union{Type{Float32},Type{Float64}})
    r.count += 1
    return r.count == 2 ? T(Inf) : randn(r.inner, T)
end

function collide!(species, grid, dt; kwargs...)
    return collide_coulomb!(
        species,
        grid,
        dt;
        coulomb_log = one(eltype(grid.dx)),
        epsilon0 = 1,
        weight_to_number = 1,
        kwargs...,
    )
end

@testset "Local Coulomb physical pair oracle" begin
    for T in (Float32, Float64),
        D = 1:3,
        masses in ((1, 1), (1, 4), (1836, 1)),
        charge in (1, 2, -2)

        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        a = population(T, D, 1; m = masses[1], q = -1)
        b = population(T, D, 1; m = masses[2], q = charge)
        foreach(v -> fill!(v, zero(T)), a.v)
        foreach(v -> fill!(v, zero(T)), b.v)
        a.v[3][1] = one(T)
        dt = T(0.001)
        before = conserved((a, b))
        @test collide_coulomb!(
            (a, b),
            grid,
            dt;
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = FixedRNG(),
        ) == (a, b)
        # Independent physical formula, evaluated at higher precision. The
        # fixed normal is one and azimuth is pi/4; initial relative g is +z.
        ma, mb = BigFloat(a.m), BigFloat(b.m)
        reduced = ma * mb / (ma + mb)
        variance = (BigFloat(a.q) * BigFloat(b.q))^2 * BigFloat(dt) / (8 * big(π) * reduced^2)
        delta = sqrt(variance)
        cosine = (1 - delta^2) / (1 + delta^2)
        sine = 2 * delta / (1 + delta^2)
        relative = (sine / sqrt(big(2)), sine / sqrt(big(2)), cosine)
        center = (big(0), big(0), ma / (ma + mb))
        for c = 1:3
            @test a.v[c][1] ≈ T(center[c] + mb / (ma + mb) * relative[c]) rtol = 64eps(T) atol =
                64eps(T)
            @test b.v[c][1] ≈ T(center[c] - ma / (ma + mb) * relative[c]) rtol = 64eps(T) atol =
                64eps(T)
        end
        after = conserved((a, b))
        @test maximum(abs.(after[1] - before[1])) <= 128eps(T) * max(before[2], 1)
        @test abs(after[2] - before[2]) <= 128eps(T) * max(before[2], 1)
    end
end

@testset "Local Coulomb cells, conservation and untouched metadata" begin
    for T in (Float32, Float64), D = 1:3, count in (0, 1, 2, 3, 4, 7, 16)
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        a = population(T, D, count; seed = 41)
        b = population(T, D, count + 2; m = 4, q = 2, seed = 97)
        a.weight .= T(1.5)
        b.weight .= T(1.5)
        # Separate two occupied cells, retaining odd and sparse populations.
        for ps in (a, b), p = 1:nparticles(ps)
            ps.x[1][p] = isodd(p) ? T(0.25) : T(2.25)
        end
        original = state.((a, b))
        before = conserved((a, b))
        work = CoulombCollisionWorkspace((a, b), grid)
        rng = MersenneTwister(423)
        for _ = 1:3
            @test collide!((a, b), grid, T(0.02); rng, workspace = work) == (a, b)
        end
        after = conserved((a, b))
        @test maximum(abs.(after[1] - before[1])) <= 256eps(T) * max(before[2], 1)
        @test abs(after[2] - before[2]) <= 256eps(T) * max(before[2], 1)
        for (k, ps) in enumerate((a, b))
            @test ps.x == original[k][1]
            @test (ps.weight, ps.id, ps.tag, ps.q, ps.m) == original[k][3:7]
            @test work.populations[k].cell == cell_index(ps, grid)
        end
    end
    for T in (Float32, Float64), D = 1:3
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        a, b = population(T, D, 1), population(T, D, 1; seed = 14)
        b.x[1][1] = T(1.25)
        before = state.((a, b))
        collide!((a, b), grid, T(100); rng = MersenneTwister(4))
        @test isequal(state.((a, b)), before)
        # Every canonical cell remains unchanged by a periodic box translation.
        for p = 1:1
            b.x[1][p] = a.x[1][p]
        end
        shifted = deepcopy((a, b))
        for ps in shifted, d = 1:D
            ps.x[d] .+= T(2) * grid.L[d]
        end
        collide!((a, b), grid, T(0.02); rng = MersenneTwister(9))
        collide!(shifted, grid, T(0.02); rng = MersenneTwister(9))
        @test a.v == shifted[1].v
        @test b.v == shifted[2].v
        # Zero weights do not dilute cell counts or consume a pairing slot.
        z = population(T, D, 3)
        z.weight .= 0
        before = state(z)
        collide!((z,), grid, T(1); rng = MersenneTwister(8))
        @test isequal(state(z), before)
    end
end

@testset "Local Coulomb unequal-weight expectation for a fixed physical pair" begin
    for T in (Float32, Float64)
        grid = FourierGrid((4,), (T(4),))
        initial_a, initial_b = population(T, 1, 1; m = 1), population(T, 1, 1; m = 4)
        foreach(v -> fill!(v, zero(T)), initial_a.v)
        foreach(v -> fill!(v, zero(T)), initial_b.v)
        initial_a.v[3][1] = one(T)
        initial_a.weight[1], initial_b.weight[1] = T(4), one(T)
        before = conserved((initial_a, initial_b))
        sum_momentum, sum_energy = zeros(BigFloat, 3), big(0)
        for uniform in (0.125, 0.375, 0.625, 0.875)
            a, b = deepcopy((initial_a, initial_b))
            # Azimuth is held fixed; exactly one of four quadrature nodes
            # accepts the heavier macroparticle's probability 1/4.
            rng = FixedRNG(1.0, [0.125, uniform, uniform])
            collide_coulomb!(
                (a, b),
                grid,
                T(0.01);
                coulomb_log = T[0 1; 1 0],
                epsilon0 = 1,
                weight_to_number = 1,
                rng,
            )
            momentum, energy = conserved((a, b))
            sum_momentum .+= momentum
            sum_energy += energy
        end
        @test maximum(abs.(sum_momentum / 4 - before[1])) <= 128eps(T) * before[2]
        @test abs(sum_energy / 4 - before[2]) <= 128eps(T) * before[2]
    end
end

@testset "Local Coulomb validation and global rollback" begin
    for T in (Float32, Float64)
        grid = FourierGrid((4,), (T(4),))
        original = (population(T, 1, 4), population(T, 1, 5; m = 4, seed = 2))
        for mutation in (
            ps -> (ps.m = zero(T)),
            ps -> (ps.m = T(Inf)),
            ps -> (ps.q = T(NaN)),
            ps -> (ps.x[1][end] = T(Inf)),
            ps -> (ps.v[2][end] = T(NaN)),
            ps -> (ps.weight[end] = -one(T)),
            ps -> (ps.weight[end] = T(Inf)),
        )
            species = deepcopy(original)
            mutation(species[2])
            before = state.(species)
            rng = MersenneTwister(45)
            untouched = copy(rng)
            @test_throws ArgumentError collide!(species, grid, T(0.01); rng)
            @test isequal(state.(species), before)
            @test rand(rng, UInt64) == rand(untouched, UInt64)
        end
        for keyword in (
            (; epsilon0 = 0),
            (; epsilon0 = Inf),
            (; weight_to_number = 0),
            (; weight_to_number = NaN),
            (; u_floor = 0),
            (; coulomb_log = -1),
            (; coulomb_log = T[1 2; 1 1]),
            (; coulomb_log = fill(T(NaN), 2, 2)),
        )
            species = deepcopy(original)
            before = state.(species)
            @test_throws ArgumentError collide_coulomb!(
                species,
                grid,
                T(0.01);
                coulomb_log = 1,
                epsilon0 = 1,
                weight_to_number = 1,
                keyword...,
            )
            @test isequal(state.(species), before)
        end
        species = deepcopy(original)
        before = state.(species)
        workspace = CoulombCollisionWorkspace(species, grid)
        @test_throws ArgumentError collide!(
            species,
            grid,
            T(0.001);
            workspace,
            rng = LateInvalidRNG(MersenneTwister(1), 0),
        )
        @test isequal(state.(species), before)
        @test collide!(species, grid, T(0.001); workspace, rng = MersenneTwister(1)) === species
        @test_throws ArgumentError collide!((species[1], species[1]), grid, T(0.01))
        @test_throws DimensionMismatch collide!(
            species,
            grid,
            T(0.01);
            workspace = CoulombCollisionWorkspace((species[1],), grid),
        )
        @test_throws DimensionMismatch collide_coulomb!(
            species,
            grid,
            T(0.01);
            coulomb_log = ones(T, 3, 3),
            epsilon0 = 1,
            weight_to_number = 1,
        )
        @test_throws ArgumentError collide!(species, grid, -one(T))
        @test_throws ArgumentError collide!(species, grid, T(Inf))
        before = state.(species)
        @test collide!(species, grid, zero(T); workspace) === species
        @test isequal(state.(species), before)
        # Public workspace storage can be tampered with: reject overlap before
        # either the staging copy or any live particle mutation.
        bad = CoulombCollisionWorkspace(species, grid)
        first = bad.populations[1]
        bad.populations[1] =
            H._CoulombCellWork(species[1].v, first.cell, first.index, first.offset, first.cursor)
        @test_throws ArgumentError collide!(species, grid, T(0.01); workspace = bad)
        @test isequal(state.(species), before)
    end
end

@testset "Local Coulomb RNG continuation and reusable scratch" begin
    for T in (Float32, Float64), D = 1:3
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        species = [population(T, D, 7), population(T, D, 3; m = 4, q = 2, seed = 5)]
        workspace = CoulombCollisionWorkspace(species, grid)
        rng = MersenneTwister(817)
        for _ = 1:3
            collide!(species, grid, T(0.001); workspace, rng)
        end
        io = IOBuffer()
        serialize(io, (species, rng))
        seekstart(io)
        restored, restored_rng = deserialize(io)
        restored_workspace = CoulombCollisionWorkspace(restored, grid)
        for _ = 1:3
            collide!(species, grid, T(0.001); workspace, rng)
            collide!(restored, grid, T(0.001); workspace = restored_workspace, rng = restored_rng)
        end
        @test isequal(state.(species), state.(restored))
        @test rand(rng, UInt64) == rand(restored_rng, UInt64)
        # Growth changes only derived workspace capacities; no stale bin entries.
        ps = species[1]
        for array in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            resize!(array, 9)
        end
        for x in ps.x
            x[8:9] .= T(0.25)
        end
        for v in ps.v
            v[8:9] .= zero(T)
        end
        ps.weight[8:9] .= one(T)
        ps.id[8:9] .= UInt64[8, 9]
        ps.tag[8:9] .= UInt32(0)
        @test collide!(species, grid, T(0.001); workspace, rng) === species
        @test length(workspace.populations[1].index) == 9
        @test workspace.populations[1].cell == cell_index(ps, grid)
    end
end

@testset "Coulomb SI and normalized full collision equivalence" begin
    for T in (Float32, Float64), D = 1:3, floor_multiplier in (1, 10)
        units = PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = 1.6726e-27)
        normalization = coulomb_normalization(units)
        L0, T0, V0 =
            to_SI(1.0, :length, units), to_SI(1.0, :time, units), to_SI(1.0, :velocity, units)
        normalized_grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        si_grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4L0), D))
        normalized = (population(T, D, 1; q = -1, m = 1 / 1836), population(T, D, 1; q = 2, m = 4))
        # Distinct, representable velocities ensure this fixture scatters.
        # Identical loader seeds alone would exercise a zero-relative-speed no-op.
        for (ps, velocity) in zip(normalized, ((1.0, -2.0, 0.5), (-0.5, 0.75, -1.0))), c = 1:3
            ps.v[c][1] = T(velocity[c])
        end
        relative_speed = sqrt(sum(c -> (normalized[1].v[c][1] - normalized[2].v[c][1])^2, 1:3))
        floor_speed = floor_multiplier == 1 ? T(1e-3) : T(floor_multiplier) * relative_speed
        @test relative_speed > zero(T)
        si = deepcopy(normalized)
        for (a, b) in zip(normalized, si)
            set_density_weight!(a, 1.0, normalized_grid)
            for d = 1:D
                b.x[d] .= to_SI(a.x[d], :length, units)
            end
            for c = 1:3
                b.v[c] .= to_SI(a.v[c], :velocity, units)
            end
            b.m, b.q = a.m * units.mi, a.q * units.e
            # D-dimensional weights represent real count per omitted SI
            # transverse measure. Physical density is w_SI/prod(dx_SI).
            b.weight .= a.weight .* (units.n0 * L0^D)
        end
        dt = T(0.03)
        # Scaling the log compensates the active-floor case so both regimes
        # retain measurable, weak-angle scattering. The enlarged log is a unit
        # conversion stress, not a proposed physical logarithm.
        logarithm_value = 20.0 * floor_multiplier^3
        logarithm = [0.0 logarithm_value; logarithm_value 0.0]
        before_normalized, before_si =
            map(ps -> map(copy, ps.v), normalized), map(ps -> map(copy, ps.v), si)
        mu_stored =
            BigFloat(normalized[1].m) * BigFloat(normalized[2].m) /
            (BigFloat(normalized[1].m) + BigFloat(normalized[2].m))
        variance =
            (BigFloat(normalized[1].q) * BigFloat(normalized[2].q))^2 *
            BigFloat(logarithm[1, 2]) *
            BigFloat(dt) *
            BigFloat(normalization.weight_to_number) *
            BigFloat(normalized[2].weight[1]) / (
                8 *
                BigFloat(pi) *
                BigFloat(normalization.epsilon0)^2 *
                mu_stored^2 *
                BigFloat(prod(normalized_grid.dx)) *
                BigFloat(max(relative_speed, floor_speed))^3
            )
        @test 0 < variance < 1
        collide_coulomb!(
            normalized,
            normalized_grid,
            dt;
            coulomb_log = logarithm,
            normalization...,
            u_floor = floor_speed,
            rng = FixedRNG(),
        )
        collide_coulomb!(
            si,
            si_grid,
            dt * T0;
            coulomb_log = logarithm,
            epsilon0 = 1 / (units.mu0 * units.c^2),
            weight_to_number = 1,
            u_floor = to_SI(floor_speed, :velocity, units),
            rng = FixedRNG(),
        )
        for s = 1:2
            @test any(c -> normalized[s].v[c] != before_normalized[s][c], 1:3)
            @test any(c -> si[s].v[c] != before_si[s][c], 1:3)
        end
        for (a, b) in zip(normalized, si), c = 1:3
            @test a.v[c] ≈ to_normalized(b.v[c], :velocity, units) rtol = 512eps(T) atol = 512eps(T)
        end
        # Independently combine every dimensional factor in the generator.
        ma, mb = BigFloat(1 / 1836), big(4)
        mu = ma * mb / (ma + mb)
        stored_weight = BigFloat(4^D)
        normalized_rate =
            4 * BigFloat(normalization.weight_to_number) * stored_weight /
            (BigFloat(normalization.epsilon0)^2 * mu^2)
        physical_rate =
            (BigFloat(units.e) * BigFloat(2units.e))^2 *
            (stored_weight * BigFloat(units.n0) * BigFloat(L0)^D) *
            BigFloat(T0) / (
                (inv(BigFloat(units.mu0) * BigFloat(units.c)^2))^2 *
                (mu * BigFloat(units.mi))^2 *
                BigFloat(V0)^3 *
                BigFloat(L0)^D
            )
        @test normalized_rate ≈ physical_rate rtol = 2e-14
    end
end

end
