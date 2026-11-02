module TestYeeLongitudinalPIC
using HybridPlasmaPIC, Test, Serialization
import JET
include("yee_charge_oracle.jl")
const O = YeeChargeOracle
const H = HybridPlasmaPIC

function population(::Type{T}, xs; q = one(T), m = one(T), v = zero(T), w = one(T)) where {T}
    ps = ParticleSet{1,T}(length(xs); q, m)
    ps.x[1] .= xs
    ps.v[1] .= v
    ps.weight .= w
    return ps
end
particle_snapshot(ps) =
    (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m)
function snapshot(st, species)
    return (
        map(particle_snapshot, species),
        map(copy, st.fields.E),
        map(copy, st.fields.B),
        st.time[],
        st.step[],
        st.initialized[],
        st.last_dt[],
        st.right_plate_charge[],
        map(copy, st.last_impulse),
        st.last_work[],
        st.last_energy_defect[],
    )
end

# Independent volume quadrature of the physical cloud against closed nodal hats.
function initial_reference(st, species)
    n, L, a = st.g.n[1], st.g.L[1], st.shape.half_width
    return [
        O.rat(st.left_plate_charge) +
        O.rat(st.background_charge) * (2e + 1) * O.rat(L) / (2n) +
        sum(
            (
                O.rat(ps.q) * O.rat(ps.weight[p]) * sum(O.cloud_nodes(ps.x[1][p], n, L, a)[1:e+1])
                for ps in species for p in eachindex(ps.weight)
            );
            init = zero(O.Q),
        ) for e = 0:n-1
    ]
end
function field_reference(field, x, st)
    return O.rat(st.g.dx[1]) * sum(
        O.rat(field[e+1]) * O.cloud_edge(x, e, st.g.n[1], st.g.L[1], st.shape.half_width) for
        e = 0:st.g.n[1]-1
    )
end
function step_reference(st, species, dt)
    T = eltype(st.g.L)
    old = [copy(ps.x[1]) for ps in species]
    u = [
        T[
            T(
                O.rat(ps.v[1][p]) +
                O.rat(dt) / 2 * O.rat(ps.q) / O.rat(ps.m) *
                field_reference(st.fields.E[1], ps.x[1][p], st),
            ) for p in eachindex(ps.weight)
        ] for ps in species
    ]
    y = [
        T[T(O.rat(x) + O.rat(dt) * O.rat(v)) for (x, v) in zip(xs, vs)] for (xs, vs) in zip(old, u)
    ]
    impulse = T.(O.cloud_impulse(species, old, y, st.g.n[1], st.g.L[1], st.shape.half_width))
    field = T[T(O.rat(e) - O.rat(i)) for (e, i) in zip(st.fields.E[1], impulse)]
    v = [
        T[
            T(
                O.rat(us[p]) +
                O.rat(dt) / 2 * O.rat(ps.q) / O.rat(ps.m) * field_reference(field, ys[p], st),
            ) for p in eachindex(ps.weight)
        ] for (ps, us, ys) in zip(species, u, y)
    ]
    return y, v, impulse, field
end

@testset "Longitudinal original-source initialization and Gauss plates" begin
    for T in (Float32, Float64), W in (Float32, Float64)
        g = YeeGrid((10,), (T(0.1),))
        a = W(0.023)
        species = [
            population(T, T[0.004, 0.041]; q = T(0.7), w = T[0.2, 0.5]),
            population(T, T[0.041, 0.096]; q = T(-0.4), m = T(2), w = T[0.8, 0.1]),
        ]
        st = YeeLongitudinalPIC(
            g;
            shape = YeeTriangularCloud(a),
            background_charge = W(0.37),
            left_plate_charge = W(-0.16),
        )
        exact = initial_reference(st, species)
        identities = st.fields.E[1], species[1].x[1], species[2].v[1]
        @test init!(st, species) === st
        for (candidate, ps) in zip(st._stage[].candidates, species)
            @test all(a -> a isa Vector, H._yl_particle_arrays(candidate))
            @test all(
                pair -> !Base.mightalias(pair...),
                zip(H._yl_particle_arrays(candidate), H._yl_particle_arrays(ps)),
            )
        end
        @test typeof(st.shape.half_width) === W
        @test typeof(st.background_charge) === W
        @test all(O.rounded(value, reference) for (value, reference) in zip(st.fields.E[1], exact))
        right =
            -O.rat(st.left_plate_charge) - O.rat(st.background_charge) * O.rat(g.L[1]) -
            sum(O.rat(ps.q) * sum(O.rat, ps.weight) for ps in species)
        @test st.right_plate_charge[] == right
        d = yee_longitudinal_diagnostics(st, species)
        charge = O.cloud_charge(species, [ps.x[1] for ps in species], 10, g.L[1], a)
        bg = O.rat(st.background_charge) * O.rat(g.L[1]) / 10
        charge[1] += bg / 2 + O.rat(st.left_plate_charge)
        charge[end] += bg / 2 + right
        charge[2:end-1] .+= bg
        residual = O.incidence(st.fields.E[1]) - charge
        @test all(
            O.rounded(value, reference) for (value, reference) in zip(d.gauss_residual, residual)
        )
        errors = abs.(O.rat.(st.fields.E[1]) - exact)
        budgets = [errors[1]; errors[1:end-1] + errors[2:end]; errors[end]]
        @test all(abs.(residual) .<= budgets)
        @test st.fields.E[1] === identities[1]
        @test species[1].x[1] === identities[2] && species[2].v[1] === identities[3]
        before = snapshot(st, species)
        species[2].m = -one(T)
        invalid = snapshot(st, species)
        @test_throws ArgumentError init!(st, species)
        @test isequal(snapshot(st, species), invalid)
        species[2].m = T(2)
        @test isequal(snapshot(st, species), before)
    end
    # Nodal charge may exceed IEEE range while its complete Gauss field is finite.
    for T in (Float32, Float64)
        ps = population(T, T[0.5]; q = floatmax(T), w = T(2))
        st = YeeLongitudinalPIC(
            YeeGrid((2,), (one(T),));
            shape = YeeTriangularCloud(T(0.125)),
            left_plate_charge = -floatmax(T),
        )
        exact = initial_reference(st, [ps])
        nodal = O.cloud_charge([ps], [ps.x[1]], 2, one(T), T(0.125))
        @test nodal[2] > O.rat(floatmax(T))
        init!(st, ps)
        @test all(isfinite, st.fields.E[1])
        @test all(O.rounded(value, reference) for (value, reference) in zip(st.fields.E[1], exact))
    end
    @test_throws ArgumentError YeeLongitudinalPIC(
        YeeGrid((4,), (1.0,); boundaries = (:periodic,));
        shape = YeeTriangularCloud(0.2),
    )
end

@testset "Independent synchronized two-kick sequence and continuity" begin
    for T in (Float32, Float64)
        g = YeeGrid((8,), (one(T),))
        species = [
            population(T, T[0.31, 0.57]; q = T(0.3), v = T[0.04, -0.03], w = T[0.2, 0.4]),
            population(T, T[0.43]; q = T(-0.5), m = T(2), v = T(0.02), w = T(0.1)),
        ]
        st = YeeLongitudinalPIC(g; shape = YeeTriangularCloud(T(0.21)), left_plate_charge = T(0.17))
        init!(st, species)
        dt = T(0.01)
        y, v, _, field = step_reference(st, species, dt)
        oldpositions = [copy(ps.x[1]) for ps in species]
        oldcharge = O.cloud_charge(species, oldpositions, 8, one(T), T(0.21))
        oldfield, oldK, oldU = copy(st.fields.E[1]), kinetic_energy(species), em_field_energy(st)
        @test step!(st, species, dt) === st
        for (ps, ys, vs) in zip(species, y, v)
            @test ps.x[1] ≈ ys rtol = 32eps(T)
            @test ps.v[1] ≈ vs rtol = 32eps(T)
        end
        actualpath =
            O.cloud_impulse(species, oldpositions, [ps.x[1] for ps in species], 8, one(T), T(0.21))
        @test all(O.rounded(value, exact) for (value, exact) in zip(st.last_impulse[1], actualpath))
        @test all(
            O.rounded(value, O.rat(old) - O.rat(current)) for
            (value, old, current) in zip(st.fields.E[1], oldfield, st.last_impulse[1])
        )
        @test st.fields.E[1] ≈ field rtol = 32eps(T)
        newcharge = O.cloud_charge(species, [ps.x[1] for ps in species], 8, one(T), T(0.21))
        # Exact physical endpoint oracle checks the measured independently rounded current.
        continuity = newcharge - oldcharge + O.incidence(st.last_impulse[1])
        @test maximum(abs, continuity) <=
              128 * O.rat(eps(T)) * max(maximum(abs, newcharge), one(O.Q))
        expectedwork =
            O.rat(g.dx[1]) / 2 * sum(
                O.rat(i) * (O.rat(a) + O.rat(b)) for
                (i, a, b) in zip(st.last_impulse[1], oldfield, st.fields.E[1])
            )
        @test abs(O.rat(st.last_work[]) - expectedwork) <= 16 * O.rounding_error(T, expectedwork)
        dK = O.rat(kinetic_energy(species)) - O.rat(oldK)
        @test O.rounded(st.last_energy_defect[], dK - O.rat(st.last_work[]))
        physical_delta =
            O.rat(kinetic_energy(species)) + O.rat(em_field_energy(st)) - O.rat(oldK) - O.rat(oldU)
        @test abs(physical_delta - O.rat(st.last_energy_defect[])) <=
              64 * O.rat(eps(T)) * max(abs(O.rat(oldK)) + abs(O.rat(oldU)), one(O.Q))
        @test !iszero(st.last_energy_defect[])
        accepted = snapshot(st, species)
        @test_throws ArgumentError step!(st, species, dt / 2)
        @test isequal(snapshot(st, species), accepted)
        init!(st, species)
        @test st.step[] == 0 && iszero(st.time[]) && iszero(st.last_dt[])
        @test iszero(st.last_work[]) && iszero(st.last_energy_defect[])
    end
end

@testset "Free flight, tracers and original charge/mass factors" begin
    for T in (Float32, Float64)
        ps = population(T, T[0.4]; q = zero(T), v = T(0.125), w = T(0.7))
        st = YeeLongitudinalPIC(
            YeeGrid((4,), (one(T),));
            shape = YeeTriangularCloud(T(0.2)),
            left_plate_charge = T(0.3),
        )
        init!(st, ps)
        field = copy(st.fields.E[1])
        step!(st, ps, T(0.125))
        @test ps.x[1] == T[0.415625] && ps.v[1] == T[0.125]
        @test st.fields.E[1] == field && all(iszero, st.last_impulse[1])
        tracer = population(T, T[0.5]; q = T(0.2), w = zero(T))
        init!(st, tracer)
        step!(st, tracer, T(0.01))
        @test tracer.v[1][1] > 0 && iszero(kinetic_energy(tracer))
        @test st.fields.E[1] == field && all(iszero, st.last_impulse[1])
    end
    # Original dt/2 can lie between subnormals or below the smallest one.
    # A large q/m makes prematurely rounding it a finite response error.
    for T in (Float32, Float64), quanta in (1, 3)
        scale = exponent(floatmin(T))
        q, m = ldexp(one(T), -scale), floatmin(T)
        ps = population(T, T[0.25]; q, m, w = zero(T))
        st = YeeLongitudinalPIC(
            YeeGrid((2,), (one(T),));
            shape = YeeTriangularCloud(T(0.125)),
            left_plate_charge = m,
        )
        init!(st, ps)
        dt = T(quanta) * nextfloat(zero(T))
        expected = O.rat(dt) * O.rat(q) * O.rat(m) / O.rat(ps.m)
        step!(st, ps, dt)
        @test O.rounded(ps.v[1][1], expected)
        @test ps.x[1] == T[0.25] && st.time[] === dt
    end
    # q/m overflows, but original-factor dt*q*E/m remains finite.
    ps = population(Float64, [0.5]; q = floatmax(Float64), m = floatmin(Float64), w = 0.0)
    st = YeeLongitudinalPIC(
        YeeGrid((4,), (1.0,));
        shape = YeeTriangularCloud(0.2),
        left_plate_charge = 1e-308,
    )
    init!(st, ps)
    @test isinf(ps.q / ps.m)
    step!(st, ps, 1e-310)
    @test isfinite(ps.v[1][1]) && ps.v[1][1] > 0
    # q/m underflows, but its original factors with E and dt have finite response.
    ps = population(Float64, [5e307]; q = nextfloat(0.0), m = floatmax(Float64), w = 0.0)
    st = YeeLongitudinalPIC(
        YeeGrid((4,), (1e308,));
        shape = YeeTriangularCloud(1e307),
        left_plate_charge = 1e308,
    )
    init!(st, ps)
    @test iszero(ps.q / ps.m)
    step!(st, ps, 1e150)
    @test isfinite(ps.v[1][1]) && ps.v[1][1] > 0
end

@testset "Nonzero rest control, temporal convergence and strict stability guard" begin
    # Within the left cell, this cloud has g=(2,0). With m=2 and b=-2,
    # its local force is -(x-1/4): an independently known nonzero omega=1
    # oscillator while its entire support remains inside that cell.
    function oscillator(dt, duration)
        ps = population(Float64, [0.3]; m = 2.0, v = 0.03)
        st = YeeLongitudinalPIC(
            YeeGrid((2,), (1.0,));
            shape = YeeTriangularCloud(0.1),
            background_charge = -2.0,
            left_plate_charge = 0.0,
        )
        init!(st, ps)
        for _ = 1:round(Int, duration / dt)
            step!(st, ps, dt)
        end
        expectedx = 0.25 + 0.05cos(duration) + 0.03sin(duration)
        expectedv = -0.05sin(duration) + 0.03cos(duration)
        return hypot(ps.x[1][1] - expectedx, ps.v[1][1] - expectedv)
    end
    errors = [oscillator(dt, 0.4) for dt in (0.1, 0.05, 0.025)]
    @test errors[1] > errors[2] > errors[3] > 0
    @test 3.5 < errors[1] / errors[2] < 4.5
    @test 3.5 < errors[2] / errors[3] < 4.5
    # h=dx=1/8, a=1/32 and sum(q^2*w/m)=1/8 give rest stiffness 1.
    # At dt=2 the rest bound is exactly 4, while the impact tube is 1/4<1/2.
    # Thus this rejection isolates the strict rest guard from wall admission.
    species = [population(Float64, [0.5]; m = 16.0), population(Float64, [0.5]; q = -1.0, m = 16.0)]
    st = YeeLongitudinalPIC(YeeGrid((8,), (1.0,)); shape = YeeTriangularCloud(0.03125))
    init!(st, species)
    @test all(iszero, st.fields.E[1])
    @test_throws ArgumentError step!(st, species, 2.0)
    @test H._yl_admit_step(st, species, 0.1) === nothing
end

@testset "Strict impact tube, source/alias rejection and late atomic failure" begin
    ps = population(Float64, [0.05]; q = 0.0, v = -0.5)
    st = YeeLongitudinalPIC(YeeGrid((4,), (1.0,)); shape = YeeTriangularCloud(0.2))
    init!(st, ps)
    before = snapshot(st, [ps])
    @test_throws ArgumentError step!(st, ps, 0.1)
    @test isequal(snapshot(st, [ps]), before)
    # A charged rest center with finite endpoint may still fail conservative admission.
    ps.v[1][1], ps.q = 0.0, 1.0
    init!(st, ps)
    before = snapshot(st, [ps])
    @test_throws ArgumentError step!(st, ps, 0.3)
    @test isequal(snapshot(st, [ps]), before)
    for mutate in (p -> (p.weight[1] = 2.0), p -> (p.id[1] += 1), p -> (p.q = -1.0))
        init!(st, ps)
        mutate(ps)
        changed = snapshot(st, [ps])
        @test_throws ArgumentError step!(st, ps, 0.001)
        @test isequal(snapshot(st, [ps]), changed)
    end
    aliased = ParticleSet{1,Float64}(
        (st.fields.E[1],),
        (zeros(4), zeros(4), zeros(4)),
        ones(4),
        UInt64[1, 2, 3, 4],
        zeros(UInt32, 4),
        1.0,
        1.0,
    )
    before = snapshot(st, [ps])
    @test_throws ArgumentError init!(st, aliased)
    @test isequal(snapshot(st, [ps]), before)
    # Finite phase space can have an unrepresentable required kinetic diagnostic.
    species = [
        population(Float64, [0.3]; q = 0.0, v = 0.1),
        population(Float64, [0.7]; q = 0.0, m = floatmax(Float64), v = 1.0, w = floatmax(Float64)),
    ]
    init!(st, species)
    before = snapshot(st, species)
    @test_throws ArgumentError step!(st, species, 0.01)
    @test isequal(snapshot(st, species), before)
    # Scratch can be reused after the late failure by explicit reinitialization.
    species[2].weight[1] = 0.0
    init!(st, species)
    step!(st, species, 0.01)
    @test st.step[] == 1
end

@testset "Longitudinal checkpoint/archive continuation and schema rejection" begin
    for scalar in (false, true), T in (Float32, Float64)
        g = YeeGrid((4,), (one(T),))
        cloud = YeeTriangularCloud(T === Float32 ? 0.2 : Float32(0.2))
        ps = population(T, T[0.4, 0.6]; q = T(0.2), v = T[0.01, -0.02], w = T[0.3, 0.4])
        particles = scalar ? ps : [ps]
        st = YeeLongitudinalPIC(
            g;
            shape = cloud,
            background_charge = T(-0.1),
            left_plate_charge = T(0.07),
        )
        init!(st, particles)
        step!(st, particles, T(0.01))
        # Preserve the actual evolved field, including an accumulated one-ulp
        # history perturbation that an initializer/projection would erase.
        st.fields.E[1][2] = nextfloat(st.fields.E[1][2])
        mktempdir() do dir
            path = joinpath(dir, "state.bin")
            save_checkpoint(path, st, particles)
            # Matched independent strided destination, preserving public identities.
            viewof(::Type{A}, n) where {A} = view(zeros(A, 2n), 1:2:2n)
            targetps = ParticleSet{1,T}(
                (viewof(T, 2),),
                ntuple(_ -> viewof(T, 2), 3),
                viewof(T, 2),
                viewof(UInt64, 2),
                viewof(UInt32, 2),
                zero(T),
                one(T),
            )
            target = YeeLongitudinalPIC(
                g;
                shape = cloud,
                background_charge = T(-0.1),
                left_plate_charge = T(0.07),
            )
            dest = scalar ? targetps : [targetps]
            ids = targetps.x[1], targetps.v[1], target.fields.E[1]
            @test load_checkpoint!(target, dest, path) === target
            @test isequal(snapshot(target, [targetps]), snapshot(st, [ps]))
            @test targetps.x[1] === ids[1] &&
                  targetps.v[1] === ids[2] &&
                  target.fields.E[1] === ids[3]
            candidate = target._stage[].candidates[1]
            @test all(a -> a isa Vector, H._yl_particle_arrays(candidate))
            @test all(
                pair -> !Base.mightalias(pair...),
                zip(H._yl_particle_arrays(candidate), H._yl_particle_arrays(targetps)),
            )
            for _ = 1:3
                step!(st, particles, T(0.01))
                step!(target, dest, T(0.01))
            end
            @test isequal(snapshot(target, [targetps]), snapshot(st, [ps]))
            archive = joinpath(dir, "archive.bin")
            archive_run(archive, st, particles; rng_seed = 7, boundary_desc = "closed event-free")
            loaded = load_archive(archive)
            @test loaded.state.integrator === :yee_longitudinal_pic
            original = deserialize(path)
            @test original.schema_version === 2
            @test original.force_policy === :original_geometry_boris_nearest_even
            @test original.force_policy_version === 1
            @test loaded.state.schema_version === 2
            for corrupt in (
                merge(original, (; schema_version = 1)),
                merge(original, (; force_policy = :rounded_gather)),
                merge(original, (; force_policy_version = 2)),
                merge(original, (; force_policy_version = true)),
                merge(original, (; force_policy_version = UInt(1))),
                (; (key => value for (key, value) in pairs(original) if key !== :force_policy)...),
                merge(original, (; shape_width_type = T)),
                merge(
                    original,
                    (;
                        shape_width_bits = original.shape_width_bits +
                                           one(original.shape_width_bits)
                    ),
                ),
                merge(original, (; arithmetic_policy = :unrestricted)),
                merge(original, (; gauss_policy = :project)),
                merge(original, (; velocity_phase = :half_time)),
                merge(original, (; right_plate_charge = original.right_plate_charge + 1)),
                merge(original, (; particle_layout = scalar ? :collection : :scalar)),
            )
                badpath = joinpath(dir, "corrupt.bin")
                open(io -> serialize(io, corrupt), badpath, "w")
                before = snapshot(target, [targetps])
                @test_throws ArgumentError load_checkpoint!(target, dest, badpath)
                @test isequal(snapshot(target, [targetps]), before)
            end
        end
    end
end

@testset "Empty collections and ordinary-vector checkpoint rebinding" begin
    g = YeeGrid((3,), (1.0,))
    st = YeeLongitudinalPIC(
        g;
        shape = YeeTriangularCloud(Float32(0.2)),
        background_charge = Float32(0.2),
        left_plate_charge = -0.1,
    )
    empty = ParticleSet{1,Float64}[]
    init!(st, empty)
    step!(st, empty, 0.01)
    @test st.step[] == 1 && all(iszero, st.last_impulse[1])
    mktempdir() do dir
        path = joinpath(dir, "empty.bin")
        save_checkpoint(path, st, empty)
        restored = YeeLongitudinalPIC(
            g;
            shape = st.shape,
            background_charge = st.background_charge,
            left_plate_charge = st.left_plate_charge,
        )
        load_checkpoint!(restored, ParticleSet{1,Float64}[], path)
        @test isequal(snapshot(restored, empty), snapshot(st, empty))
        @test isempty(restored._stage[].candidates)
        @test eltype(restored._stage[].candidates) === H._ElectrostaticHostParticle{1,Float64}
        ps = population(Float64, [0.4, 0.6]; q = -0.3, m = 2.0, w = [0.2, 0.4])
        ps.id .= UInt64[7, 11]
        ps.tag .= UInt32[3, 5]
        init!(st, [ps])
        step!(st, [ps], 0.01)
        save_checkpoint(path, st, [ps])
        targetps = ParticleSet{1,Float64}(0; q = 0, m = 1)
        identities = H._yl_particle_arrays(targetps)
        load_checkpoint!(restored, [targetps], path)
        @test isequal(snapshot(restored, [targetps]), snapshot(st, [ps]))
        @test all(pair -> pair[1] === pair[2], zip(identities, H._yl_particle_arrays(targetps)))
        @test all(a -> a isa Vector, H._yl_particle_arrays(restored._stage[].candidates[1]))
        step!(st, [ps], 0.01)
        step!(restored, [targetps], 0.01)
        @test isequal(snapshot(restored, [targetps]), snapshot(st, [ps]))
        payload = H._checkpoint_state(st, [ps])
        corrupt_record = merge(payload.species[1], (; q = 0.7))
        corrupt = merge(payload, (; species = [corrupt_record]))
        before = snapshot(restored, [targetps])
        @test_throws ArgumentError H._restore_checkpoint_state!(restored, [targetps], corrupt)
        @test isequal(snapshot(restored, [targetps]), before)
        # Incoming phase-space aliases are rejected before detach/publication.
        aliased_record = merge(
            payload.species[1],
            (; v = (payload.species[1].x[1], payload.species[1].v[2], payload.species[1].v[3])),
        )
        @test_throws ArgumentError H._restore_checkpoint_state!(
            restored,
            [targetps],
            merge(payload, (; species = [aliased_record])),
        )
        @test isequal(snapshot(restored, [targetps]), before)
    end
end

@testset "Stored transverse subnormals and rejected arithmetic contexts" begin
    for T in (Float32, Float64)
        ps = population(T, T[0.5]; q = zero(T))
        st = YeeLongitudinalPIC(YeeGrid((2,), (one(T),)); shape = YeeTriangularCloud(T(0.2)))
        init!(st, ps)
        ps.v[2][1] = nextfloat(zero(T))
        before = snapshot(st, [ps])
        @test_throws ArgumentError step!(st, ps, T(0.01))
        @test isequal(snapshot(st, [ps]), before)
        ps.v[2][1] = zero(T)
        st.fields.B[1][1] = nextfloat(zero(T))
        before = snapshot(st, [ps])
        @test_throws ArgumentError init!(st, ps)
        @test isequal(snapshot(st, [ps]), before)
        st.fields.B[1][1] = zero(T)
        U = T === Float32 ? UInt32 : UInt64
        # Bit pattern 1 is the minimal positive IEEE subnormal. Store it before
        # changing modes so this checks storage, independent of nextfloat.
        subnormal = T[reinterpret(T, one(U))]
        for (roundmode, flush) in ((RoundNearest, true), (RoundDown, false))
            saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
            before = snapshot(st, [ps])
            try
                Base.Rounding.setrounding_raw(Float64, Base.Rounding.to_fenv(roundmode))
                Base.set_zero_subnormals(flush)
                @test reinterpret(U, subnormal[1]) == one(U)
                @test !H._yl_stored_zero(subnormal[1])
                @test_throws ArgumentError init!(st, ps)
                @test_throws ArgumentError step!(st, ps, T(0.01))
            finally
                Base.Rounding.setrounding_raw(Float64, saved[1])
                Base.set_zero_subnormals(saved[2])
            end
            @test isequal(snapshot(st, [ps]), before)
        end
    end
end

@testset "Concrete longitudinal compiler diagnostics" begin
    for T in (Float32, Float64)
        st = YeeLongitudinalPIC(YeeGrid((2,), (one(T),)); shape = YeeTriangularCloud(T(0.2)))
        ps = population(T, T[0.5])
        for P in (typeof(ps), Vector{typeof(ps)})
            for (operation, types) in (
                (init!, Tuple{typeof(st),P}),
                (step!, Tuple{typeof(st),P,T}),
                (yee_longitudinal_diagnostics, Tuple{typeof(st),P}),
                (H._checkpoint_state, Tuple{typeof(st),P}),
            )
                report = JET.report_call(
                    operation,
                    types;
                    sourceinfo = :full,
                    stacktrace_types_limit = 0,
                )
                if !isempty(JET.get_reports(report))
                    show(IOContext(stdout, :limit => false), MIME"text/plain"(), report)
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end

end
