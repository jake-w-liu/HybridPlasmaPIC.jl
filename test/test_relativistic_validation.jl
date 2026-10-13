using HybridPlasmaPIC, Test, Random, Serialization

# Each stored tuple is strictly superluminal by an exact rational comparison,
# although the original divided/squared/summed Float calculation is below one.
const _RELATIVISTIC_BOUNDARY_WITNESSES = (
    (Float16, (0.6143, 0.621, -0.4868)),
    (Float32, (-0.86362565, -0.49116498, 0.113612115)),
    (Float64, (-0.474852586251129, 0.8763033506698198, -0.0812862776577999)),
)

_relativistic_exact_norm(v) = sum(abs2, Rational{BigInt}.(v))
function _relativistic_test_bytes(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

@testset "Relativistic validation resolves the stored light-speed boundary" begin
    permutations = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
    for (T, raw) in _RELATIVISTIC_BOUNDARY_WITNESSES
        v = T.(raw)
        @test _relativistic_exact_norm(v) > 1
        @test sum(abs2, v) < one(T)
        half_range = exponent(floatmax(T)) ÷ 2
        for D = 1:3, e in (-half_range, 0, half_range), order in permutations, mask = 0:7
            c = ldexp(one(T), e)
            u = ntuple(d -> ldexp(v[order[d]], e) * (iszero(mask & (1 << (d - 1))) ? 1 : -1), 3)
            ps = ParticleSet{D,T}(2)
            foreach(x -> fill!(x, T(0.2)), ps.x)
            for d = 1:3
                ps.v[d][2] = u[d]
            end
            @test _relativistic_exact_norm(u) > Rational{BigInt}(c)^2
            before = _relativistic_test_bytes(ps)
            @test_throws ArgumentError HybridPlasmaPIC._require_subluminal_velocities(ps, c, "test")
            @test _relativistic_test_bytes(ps) == before

            # Moving the largest component one stored ULP toward zero yields a
            # neighboring physical state. Preserve every one the old guard accepts.
            largest = argmax(abs.(u))
            uvalid = Base.setindex(u, copysign(prevfloat(abs(u[largest])), u[largest]), largest)
            @test _relativistic_exact_norm(uvalid) < Rational{BigInt}(c)^2
            for d = 1:3
                ps.v[d][2] = uvalid[d]
            end
            ratios = uvalid ./ c
            rounded = ratios[1]^2 + ratios[2]^2 + ratios[3]^2
            if rounded < one(T)
                @test isnothing(HybridPlasmaPIC._require_subluminal_velocities(ps, c, "test"))
            else
                @test_throws ArgumentError HybridPlasmaPIC._require_subluminal_velocities(
                    ps,
                    c,
                    "test",
                )
            end
        end

        rng = Xoshiro(927)
        ps = ParticleSet{1,T}(1)
        for _ = 1:2000
            direction = randn(rng, 3)
            direction ./= sqrt(sum(abs2, direction))
            u = T.(direction) .* (one(T) + T(rand(rng, -10:10)) * eps(T))
            for d = 1:3
                ps.v[d][1] = u[d]
            end
            rounded = u[1]^2 + u[2]^2 + u[3]^2
            valid = rounded < one(T) && _relativistic_exact_norm(u) < 1
            if valid
                @test isnothing(HybridPlasmaPIC._require_subluminal_velocities(ps, one(T), "test"))
            else
                @test_throws ArgumentError HybridPlasmaPIC._require_subluminal_velocities(
                    ps,
                    one(T),
                    "test",
                )
            end
        end
    end
end

@testset "Relativistic validation retains bits lost by division" begin
    speeds = (Float16(1.246), 1.3265638f0, 0.7619480428727838)
    for ((T, raw), c) in zip(_RELATIVISTIC_BOUNDARY_WITNESSES, speeds)
        u = T.(raw) .* c
        @test _relativistic_exact_norm(u) > Rational{BigInt}(c)^2
        @test _relativistic_exact_norm(u ./ c) < 1
        for D = 1:3
            ps = ParticleSet{D,T}(1)
            for d = 1:3
                ps.v[d][1] = u[d]
            end
            @test_throws ArgumentError HybridPlasmaPIC._require_subluminal_velocities(ps, c, "test")
        end
    end
end

@testset "Relativistic validation preserves BigFloat settings and stored bits" begin
    stored = setprecision(BigFloat, 256) do
        BigFloat.(_RELATIVISTIC_BOUNDARY_WITNESSES[3][2])
    end
    for bits in (11, 24, 53, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                ps = ParticleSet{1,BigFloat}(1)
                for d = 1:3
                    ps.v[d][1] = stored[d]
                end
                @test_throws ArgumentError HybridPlasmaPIC._require_subluminal_velocities(
                    ps,
                    one(BigFloat),
                    "test",
                )
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
                @test all(precision(ps.v[d][1]) == 256 for d = 1:3)
                ps.v[1][1] = BigFloat(0.5)
                ps.v[2][1] = ps.v[3][1] = zero(BigFloat)
                @test isnothing(
                    HybridPlasmaPIC._require_subluminal_velocities(ps, one(BigFloat), "test"),
                )
            end
        end
    end
end

@testset "Radiation rejects superluminal stored velocities before cooling" begin
    for (T, raw) in _RELATIVISTIC_BOUNDARY_WITNESSES, D = 1:3
        ps = ParticleSet{D,T}(2)
        foreach(x -> fill!(x, T(0.2)), ps.x)
        ps.v[1][1] = T(0.25)
        for d = 1:3
            ps.v[d][2] = T(raw[d])
        end
        before = _relativistic_test_bytes(ps)
        z = (zero(T), zero(T), zero(T))
        @test_throws ArgumentError apply_radiation_reaction!(
            ps,
            (zero(T), one(T), zero(T)),
            z,
            T(0.1);
            K = one(T),
        )
        @test _relativistic_test_bytes(ps) == before
        @test_throws ArgumentError apply_radiation_reaction!(ps, z, z, T(0.1); K = one(T))
        @test _relativistic_test_bytes(ps) == before
    end
end

@testset "Relativistic PIC and checkpoint boundary rejection is transactional" begin
    for (T, raw) in _RELATIVISTIC_BOUNDARY_WITNESSES[2:3], D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        for State in (D == 1 ? (EMPIC1D, EMPIC) : (EMPIC,)), mobile in (false, true)
            s = State(g, 2; relativistic = true, c = one(T), mobile, mi = T(16))
            electrons = ParticleSet{D,T}(2; q = -one(T))
            ions = mobile ? ParticleSet{D,T}(2; m = T(16)) : nothing
            for ps in (electrons, ions)
                ps === nothing && continue
                foreach(x -> fill!(x, T(0.2)), ps.x)
            end
            init_empic!(s, electrons, ions)
            mktempdir() do dir
                path = joinpath(dir, "valid.bin")
                save_checkpoint(path, s, electrons, ions)
                saved = read(path)
                for species in (mobile ? (:electrons, :ions) : (:electrons,))
                    ps = species === :electrons ? electrons : ions
                    for d = 1:3
                        ps.v[d][2] = T(raw[d])
                    end
                    before = _relativistic_test_bytes((s, electrons, ions))
                    @test_throws ArgumentError step_empic!(s, electrons, ions, T(0.001))
                    @test _relativistic_test_bytes((s, electrons, ions)) == before
                    @test_throws ArgumentError save_checkpoint(path, s, electrons, ions)
                    @test read(path) == saved
                    @test _relativistic_test_bytes((s, electrons, ions)) == before
                    foreach(v -> fill!(v, zero(T)), ps.v)

                    payload = deserialize(IOBuffer(saved))
                    for d = 1:3
                        getproperty(payload, species).v[d][2] = T(raw[d])
                    end
                    bad_path = joinpath(dir, "superluminal.bin")
                    serialize(bad_path, payload)
                    before = _relativistic_test_bytes((s, electrons, ions))
                    @test_throws ArgumentError load_checkpoint!(s, electrons, ions, bad_path)
                    @test _relativistic_test_bytes((s, electrons, ions)) == before
                end
                @test load_checkpoint!(s, electrons, ions, path) === s
            end
        end
    end
end
