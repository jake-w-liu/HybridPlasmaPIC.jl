module RelativisticSpeciesFactorTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC

# Independent axial E-only solution: momentum translates by q*duration*E/m.
# This does not invoke the package's Boris, priming or conversion owners.
function axial_reference(v::T, q, m, duration, electric, c) where {T}
    Q = Rational{BigInt}
    return setprecision(BigFloat, 512) do
        setrounding(BigFloat, RoundNearest) do
            V, C = BigFloat(v), BigFloat(c)
            impulse = BigFloat(Q(q) * Q(duration) * Q(electric) / Q(m))
            momentum = V / sqrt(1 - (V / C)^2) + impulse
            T(momentum / sqrt(1 + (momentum / C)^2))
        end
    end
end

function physical_bytes(es, species...)
    names = filter(!=(:transaction), fieldnames(typeof(es)))
    state = NamedTuple{names}(map(name -> getfield(es, name), names))
    io = IOBuffer()
    serialize(io, (state, species))
    return take!(io)
end

function transverse!(es, value)
    fill!(es isa EMPIC1D ? es.Ey : es.E[2], value)
    return nothing
end

field_values(es) = es isa EMPIC1D ? (es.Ex, es.Ey, es.Bz) : (es.E..., es.B...)

@testset "Relativistic stages retain original species factors" begin
    for (T, mass_power, field_power) in ((Float16, 10, 10), (Float32, 80, 60), (Float64, 600, 400))
        z, c = zero(T), one(T)
        zero3 = (z, z, z)
        # Overflowed quotient, underflowed quotient and a subnormal quotient.
        cases = (
            (one(T), nextfloat(z), 2floatmin(T), ldexp(one(T), -precision(T) - 8)),
            (
                nextfloat(z),
                ldexp(one(T), mass_power),
                ldexp(one(T), field_power),
                ldexp(one(T), field_power),
            ),
            (one(T), floatmax(T), one(T), T(0.125)),
        )
        for (q, m, dt, electric) in cases
            ps = ParticleSet{1,T}(1; q, m)
            factor = H._validated_boris_factor(ps)
            E, B = (z, electric, z), zero3
            expected = axial_reference(z, q, m, dt, electric, c)
            @test expected != z
            @test H._boris_rel(zero3..., E..., B..., factor, dt, c) == (z, expected, z)
            ep, bp = map(x -> T[x], E), map(x -> T[x], B)
            for planar in (false, true)
                v = ntuple(_ -> T[0], 3)
                if planar
                    H._prime_empic1d!(v, ep[1], ep[2], bp[3], factor, dt / 2, 1, true, c)
                else
                    H._prime_leapfrog!(v, ep, bp, factor, dt / 2, 1, true, c)
                end
                primed = axial_reference(z, q, m, -dt / 2, electric, c)
                @test map(only, v) == (z, primed, z)
            end
        end
        # A non-exact normal ratio cannot be a DD scalar-coefficient proof.
        ps = ParticleSet{1,T}(1; q = one(T), m = T(3))
        factor = H._validated_boris_factor(ps)
        E = (z, T(0.125), z)
        @test !H._rel_prime_exact_coefficient(factor)
        expected = axial_reference(z, ps.q, ps.m, -T(0.5), E[2], c)
        @test H._rel_prime_repair(zero3, zero3, E, zero3, factor, T(0.5), c) == (z, expected, z)
        ordinary = H._BorisSpeciesFactors(one(T), one(T), one(T))
        @test H._boris_rel(zero3..., E..., zero3..., ordinary, T(0.125), c) ==
              H._boris_rel(zero3..., E..., zero3..., one(T), T(0.125), c)
        neutral = H._BorisSpeciesFactors(z, nextfloat(z), z)
        @test H._boris_rel(zero3..., E..., zero3..., neutral, one(T), c) == zero3
        # Explicit qm is still the supplied scalar, including IEEE propagation.
        @test all(isnan, H._boris_rel(zero3..., E..., zero3..., T(Inf), one(T), c))
    end
end

@testset "Relativistic EMPIC light-ion phase, restart and rollback" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        g = FourierGrid((4,), (T(4),))
        mass, dt, c = nextfloat(zero(T)), 2floatmin(T), one(T)
        electric = ldexp(one(T), -precision(T) - 8)
        es = State(g, 4; mobile = true, mi = mass, relativistic = true, c)
        electrons = ParticleSet{1,T}(4; q = -one(T))
        ions = ParticleSet{1,T}(4; m = mass)
        electrons.x[1] .= T.(0:3)
        ions.x[1] .= electrons.x[1]
        init_empic!(es, electrons, ions)
        transverse!(es, electric)
        primed = axial_reference(zero(T), one(T), mass, -dt / 2, electric, c)
        expected = axial_reference(primed, one(T), mass, dt, electric, c)
        transverse!(es, one(T))
        before = physical_bytes(es, electrons, ions)
        @test_throws ArgumentError step_empic!(es, electrons, ions, dt)
        @test physical_bytes(es, electrons, ions) == before
        transverse!(es, electric)
        @test step_empic!(es, electrons, ions, dt) === es
        @test all(==(expected), ions.v[2])
        @test es.step[] == 1 && es.time[] == dt

        mktempdir() do directory
            path = joinpath(directory, "light-ion.bin")
            save_checkpoint(path, es, electrons, ions)
            fresh = State(
                FourierGrid((4,), (T(4),)),
                0;
                mobile = true,
                mi = mass,
                relativistic = true,
                c,
            )
            fresh_e = ParticleSet{1,T}(0; q = -one(T))
            fresh_i = ParticleSet{1,T}(0; m = mass)
            load_checkpoint!(fresh, fresh_e, fresh_i, path)
            @test fresh.transaction === nothing
            @test fresh_i.q == one(T) && fresh_i.m == mass
            es.n_sub = fresh.n_sub = 2
            kick_dt = dt / 2 + (2dt) / 2
            @test all(==(electric), field_values(es)[2])
            next_velocity = axial_reference(expected, one(T), mass, kick_dt, electric, c)
            step_empic!(es, electrons, ions, 2dt)
            step_empic!(fresh, fresh_e, fresh_i, 2dt)
            @test all(==(next_velocity), ions.v[2])
            for (loaded, original) in ((fresh_e, electrons), (fresh_i, ions))
                @test all(
                    name -> getfield(loaded, name) == getfield(original, name),
                    fieldnames(typeof(original)),
                )
            end
            @test field_values(fresh) == field_values(es)
            @test (fresh.time[], fresh.step[], fresh.last_dt[], fresh.last_sub_dt[]) ==
                  (es.time[], es.step[], es.last_dt[], es.last_sub_dt[])
        end
        # The original-factor result may round to c; retain C1 rejection/undo.
        transverse!(es, one(T))
        before = physical_bytes(es, electrons, ions)
        @test_throws ArgumentError step_empic!(es, electrons, ions, 2dt)
        @test physical_bytes(es, electrons, ions) == before
    end
end

@testset "Relativistic kicks preserve repeated-storage sequential reads" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        count, dt, c = 3, T(0.125), one(T)
        electric = ldexp(one(T), cld(precision(T), 2) - 3)
        v = ntuple(_ -> view(T[0], ones(Int, count)), 3)
        ps = ParticleSet{1,T}(
            (zeros(T, count),),
            v,
            ones(T, count),
            UInt64.(1:count),
            zeros(UInt32, count),
            one(T),
            one(T),
        )
        E, B = (zero(T), electric, zero(T)), ntuple(_ -> zero(T), 3)
        Ep, Bp = map(x -> fill(x, count), E), map(x -> fill(x, count), B)
        h = dt / 2
        scales = H._boris_batch_scales(v, Ep, Bp, h, 1:count)
        @test H._boris_rel_batch_certificate(h, one(T), dt, c, scales...)
        @test !H._euler_batch_storage_safe(v, v, Ep, Bp)
        expected = ntuple(_ -> zero(T), 3)
        for _ = 1:count
            expected = H._boris_rel(expected..., E..., B..., one(T), dt, c)
        end
        es = State(FourierGrid((4,), (T(4),)), count; relativistic = true, c)
        if State === EMPIC1D
            H._push_velocities!(es, ps, Ep[1], Ep[2], Bp[3], dt)
        else
            H._push_velocities!(es, ps, Ep, Bp, dt)
        end
        @test map(first, v) == expected
    end
end
end
