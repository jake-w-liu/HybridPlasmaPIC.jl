module TestYeeImpactLongitudinalPIC
using HybridPlasmaPIC, Test, Serialization
import JET
include("yee_charge_oracle.jl")
const O = YeeChargeOracle
const H = HybridPlasmaPIC

function population(::Type{T}, x, v; q = zero(T), m = one(T), w = one(T)) where {T}
    ps = ParticleSet{1,T}(length(x); q, m)
    ps.x[1] .= x
    ps.v[1] .= v
    ps.weight .= w
    return ps
end
function state(::Type{T}, particles; sigma = zero(T), width = T(0.125)) where {T}
    st = YeeImpactLongitudinalPIC(
        YeeGrid((4,), (one(T),));
        shape = YeeTriangularCloud(width),
        left_plate_charge = sigma,
        max_root_refinements = 128,
    )
    init!(st, particles)
    return st
end
function snapshot(st, particles)
    species = particles isa ParticleSet ? [particles] : particles
    return deepcopy((
        [
            (
                map(copy, ps.x),
                map(copy, ps.v),
                copy(ps.weight),
                copy(ps.id),
                copy(ps.tag),
                ps.q,
                ps.m,
            ) for ps in species
        ],
        map(copy, st.fields.E),
        map(copy, st.fields.B),
        H._yi_history_values(st),
        st.right_plate_charge[],
    ))
end

@testset "Exact ballistic cluster, source work and wall ledgers" begin
    for T in (Float32, Float64)
        ps = population(T, T[0.25, 0.75, 0.5], T[-1, 1, 0])
        st = state(T, ps)
        identities = (ps.x[1], ps.v[1], st.fields.E[1], st.last_impulse[1], st.last_work)
        oldK = O.rat(kinetic_energy(ps))
        @test oldK > 0
        step!(st, ps, T(0.5))
        @test ps.x[1] == T[0.25, 0.75, 0.5] && ps.v[1] == T[1, -1, 0]
        @test map(x -> x[], st.last_wall_count) == (1, 1)
        @test map(x -> x[], st.last_wall_impulse) == (-2 * O.rat(one(T)), 2 * O.rat(one(T)))
        @test iszero(st.last_electric_impulse[]) && iszero(st.last_energy_defect[])
        @test all(iszero, st.last_impulse[1]) && st.time[] == T(0.5)
        d = yee_impact_longitudinal_diagnostics(st, ps)
        @test d.wall_charge_transfer == d.wall_work_transfer == (zero(O.Q), zero(O.Q))
        @test all(iszero, d.gauss_residual)
        @test d.kinetic_energy == T(oldK) && iszero(d.last_physical_energy_change)
        @test all(
            pair -> pair[1] === pair[2],
            zip(identities, (ps.x[1], ps.v[1], st.fields.E[1], st.last_impulse[1], st.last_work)),
        )
        @test all(a -> a isa Vector, H._yl_particle_arrays(st._stage[].event_candidates[1]))
        @test !Base.mightalias(st._stage[].event_candidates[1].x[1], st._stage[].candidates[1].x[1])
        before = snapshot(st, ps)
        @test_throws ArgumentError step!(st, ps, T(0.25))
        @test isequal(snapshot(st, ps), before)
        # Reuse the accepted event state: a fresh zero-error tube clears,
        # while carried uncertainty beyond the nonselected wall distance rejects.
        cluster = (duration = O.rat(1 // 4), lower = O.rat(1 // 4), width = zero(O.Q))
        geom = (tau = O.rat(1 // 2),)
        @test H._yi_clear_remainder!(st, cluster, geom) == O.rat(1 // 4)
        st._stage[].carried[1][3] = O.rat(3 // 4)
        @test_throws ArgumentError H._yi_clear_remainder!(st, cluster, geom)
        @test isequal(snapshot(st, ps), before)
        st._stage[].carried[1][3] = zero(O.Q)
    end
end

# Independent physical-volume field/potential for the one-particle q=m=w=1
# fixture, never production A/g/root/converter.
# Composite Simpson has bounded work; paired 64/128 subdivisions control the
# reference relative to the measured trajectory signal, not BigFloat eps.
function physical_field(st, initial, x)
    n, L, a = st.g.n[1], st.g.L[1], st.shape.half_width
    nodes = O.cloud_nodes(x, n, L, a)
    initial_nodes = O.cloud_nodes(initial, n, L, a)
    return [
        O.rat(st.fields.E[1][e+1]) + sum(nodes[1:e+1]) - sum(initial_nodes[1:e+1]) for e = 0:n-1
    ]
end
potential(st, initial, x) = O.rat(st.g.dx[1]) / 2 * sum(abs2, physical_field(st, initial, x))
function travel_reference(st, x0, v0, terminal, subdivisions)
    U0 = potential(st, x0, x0)
    speed(x) = sqrt(BigFloat(O.rat(v0)^2 + 2 * (U0 - potential(st, x0, O.rat(x)))))
    function travel(x)
        dx = x / subdivisions
        total = inv(speed(zero(x))) + inv(speed(x))
        for j = 1:subdivisions-1
            total += (isodd(j) ? 4 : 2) / speed(dx * j)
        end
        return dx * total / 3
    end
    contact = travel(BigFloat(x0))
    @test 0 < contact < BigFloat(terminal)
    remaining = BigFloat(terminal) - contact
    lo, hi = zero(BigFloat), BigFloat(terminal) * 2
    @test travel(hi) > remaining
    for _ = 1:64
        middle = (lo + hi) / 2
        if travel(middle) < remaining
            lo = middle
        else
            hi = middle
        end
    end
    x = (lo + hi) / 2
    return (time = contact, x = x, v = speed(x))
end

@testset "Independent nonzero-force impact and temporal refinement" begin
    setprecision(BigFloat, 192) do
        T = Float64
        x0, v0, terminal = 1 / 128, -1.0, 1 / 64
        seed = population(T, [x0], [v0]; q = 1)
        initial = state(T, seed; sigma = -2.0)
        @test kinetic_energy(seed) > 0
        wall_field = physical_field(initial, x0, 0)
        force_at_wall =
            O.rat(initial.g.dx[1]) *
            sum(wall_field[e+1] * O.cloud_edge(0, e, 4, 1, 1 / 8) for e = 0:3)
        @test !iszero(force_at_wall)
        coarse = travel_reference(initial, x0, v0, terminal, 64)
        fine = travel_reference(initial, x0, v0, terminal, 128)
        oracle_error = abs(coarse.x - fine.x) + abs(coarse.v - fine.v)
        errors = BigFloat[]
        for steps in (4, 8, 16)
            ps = population(T, [x0], [v0]; q = 1)
            st = state(T, ps; sigma = -2.0)
            dt = terminal / steps
            for _ = 1:steps
                step!(st, ps, dt)
            end
            @test sum(x -> x[], st.cumulative_wall_count) == 1
            @test ps.v[1][1] > 0 && 0 < ps.x[1][1] < 1
            push!(errors, abs(BigFloat(ps.x[1][1]) - fine.x) + abs(BigFloat(ps.v[1][1]) - fine.v))
        end
        @test oracle_error < minimum(errors) / 32
        @test minimum(errors) > 1024 * eps(T)
        orders = (log2(errors[1] / errors[2]), log2(errors[2] / errors[3]))
        @test all(>(1.8), orders)
        @info "Impact temporal refinement" errors orders oracle_error
    end
    # Nonzero charged local currents/work are checked from actual accepted
    # constrained endpoints and each independently rounded stage source.
    for T in (Float32, Float64)
        ps = population(T, T[1/128, 0.5], T[-1, 0]; q = one(T), w = T[1, 1/128])
        st = state(T, ps; sigma = T(-2))
        x0, E0, K0 = copy(ps.x[1]), copy(st.fields.E[1]), H._kinetic_energy_rational([ps])
        step!(st, ps, T(1 / 64))
        event = st._stage[].event_candidates[1]
        @test st._stage[].selected_side[2] == 0 && 0 < event.x[1][2] < 1
        @test st._stage[].carried[1][2] > 0
        I1 = T.(O.cloud_impulse([ps], [x0], [event.x[1]], 4, one(T), T(1 / 8)))
        I2 = T.(O.cloud_impulse([ps], [event.x[1]], [ps.x[1]], 4, one(T), T(1 / 8)))
        @test any(!iszero, I1) && any(!iszero, I2)
        E1, E2 = st._event_fields.E[1], st.fields.E[1]
        @test all(O.rounded(e, O.rat(a) - O.rat(i)) for (e, a, i) in zip(E1, E0, I1))
        @test all(O.rounded(e, O.rat(a) - O.rat(i)) for (e, a, i) in zip(E2, E1, I2))
        exactwork =
            O.rat(st.g.dx[1]) / 2 * sum(
                O.rat(i) * (O.rat(a) + O.rat(b)) for
                (I, before, after) in ((I1, E0, E1), (I2, E1, E2)) for
                (i, a, b) in zip(I, before, after)
            )
        @test O.rounded(st.last_work[], exactwork)
        @test all(
            O.rounded(i, O.rat(a) + O.rat(b)) for (i, a, b) in zip(st.last_impulse[1], I1, I2)
        )
        @test O.rounded(st.last_energy_defect[], H._kinetic_energy_rational([ps]) - K0 - exactwork)
        @test st.last_wall_count[1][] == 1 && st.last_wall_impulse[1][] < 0
    end
end

@testset "Branches, late rejection and reusable scratch" begin
    q(x) = O.rat(x)
    @test H._yi_first_root((q(1), q(-4), q(4)), q(2), 1, 1) !== nothing
    @test_throws ArgumentError H._yi_first_root((q(1), q(-2), q(2)), q(2), 1, 1)
    @test_throws ArgumentError H._yi_first_root((q(1), q(-1), q(0)), q(1), 1, 1)
    for T in (Float32, Float64)
        ps = population(T, T[0.25], T[3])
        st = state(T, ps)
        before = snapshot(st, ps)
        @test_throws ArgumentError step!(st, ps, one(T))
        @test isequal(snapshot(st, ps), before)
        ps.v[1][1] = T(0.25)
        step!(st, ps, one(T))
        @test st.step[] == 1 && ps.x[1][1] == T(0.5)
        budgeted = YeeImpactLongitudinalPIC(
            YeeGrid((4,), (one(T),));
            shape = YeeTriangularCloud(T(0.125)),
            max_root_refinements = 1,
        )
        moving = population(T, T[0.25], T[-3])
        init!(budgeted, moving)
        budget_before = snapshot(budgeted, moving)
        @test_throws ArgumentError step!(budgeted, moving, T(0.25))
        @test isequal(snapshot(budgeted, moving), budget_before)

        # Close charged contacts require a representable common event endpoint;
        # overlapping physical enclosures alone cannot authorize a snap.
        close = population(T, T[1/128, 1/128+1/2^20], T[-1, -1]; q = one(T))
        unresolved = state(T, close; sigma = T(-2))
        close_before = snapshot(unresolved, close)
        @test_throws ArgumentError step!(unresolved, close, T(1 / 64))
        @test isequal(snapshot(unresolved, close), close_before)

        # Original subnormal duration: a zero-weight tracer keeps finite energy
        # while original q/m and the postreflection final kick cross overflow.
        # At the wall half the cloud is outside, so its gathered field is 8,
        # half the stored 16. The doubled charge reaches the late overflow.
        eta = reinterpret(T, T === Float32 ? UInt32(1) : UInt64(1))
        x = ldexp(one(T), exponent(eta) + exponent(floatmax(T)) - 1)
        tracer = population(
            T,
            T[x],
            T[-prevfloat(floatmax(T))];
            q = 2 * inv(floatmin(T)),
            m = floatmin(T),
            w = zero(T),
        )
        failed = state(T, tracer; sigma = T(16))
        original = snapshot(failed, tracer)
        @test_throws ArgumentError step!(failed, tracer, eta)
        @test isequal(snapshot(failed, tracer), original)
        @test failed._stage[].event_candidates[1].x[1][1] == zero(T)
        @test failed._stage[].candidates[1].x[1][1] > x
        tracer.v[1][1] = zero(T)
        step!(failed, tracer, eta)
        @test failed.step[] == 1 && all(isfinite, tracer.v[1])
    end
end

@testset "Detached complete restart and storage" begin
    for T in (Float32, Float64)
        ps = population(T, T[0.25], T[-1])
        st = state(T, ps; width = Float64(0.125))
        step!(st, ps, T(0.5))
        saved = deepcopy(H._checkpoint_state(st, ps))
        target = ParticleSet{1,T}(
            (view(zeros(T, 2), 1:2:2),),
            ntuple(_ -> view(zeros(T, 2), 1:2:2), 3),
            view(zeros(T, 2), 1:2:2),
            view(UInt64[1, 0], 1:2:2),
            view(UInt32[0, 0], 1:2:2),
            zero(T),
            one(T),
        )
        restored = state(T, population(T, T[0.5], T[0]); width = Float64(0.125))
        # A valid incoming snapshot replaces invalid destination history values.
        restored.last_max_snap[] = -one(O.Q)
        H._restore_checkpoint_state!(restored, target, saved)
        @test snapshot(restored, target) == snapshot(st, ps)
        @test all(a -> a isa Vector, H._yl_particle_arrays(restored._stage[].event_candidates[1]))
        # Both tangential leaves are physically zero but carry distinct bits:
        # a sequential live-source swap loses the saved negative zero.
        restored.fields.E[2] .= -zero(T)
        restored.fields.E[3] .= zero(T)
        swapped = merge(
            saved,
            (
                fields = (
                    E = (saved.fields.E[1], restored.fields.E[3], restored.fields.E[2]),
                    B = saved.fields.B,
                ),
            ),
        )
        H._restore_checkpoint_state!(restored, target, swapped)
        @test all(x -> isequal(x, zero(T)), restored.fields.E[2])
        @test all(x -> isequal(x, -zero(T)), restored.fields.E[3])
        H._restore_checkpoint_state!(restored, target, saved)
        rejected = snapshot(restored, target)
        for corrupt in (
            merge(saved, (schema_version = 2,)),
            merge(saved, (coverage_policy_version = true,)),
            merge(saved, (last_max_snap = -one(O.Q),)),
            merge(saved, (last_position_error_bound = BigInt(1) // BigInt(0),)),
            merge(saved, (max_root_refinements = 1,)),
            merge(saved, (shape_width_type = Float32,)),
            merge(saved, (cumulative_wall_count = (typemax(Int), 0),)),
            merge(saved, (right_plate_charge = saved.right_plate_charge + 1,)),
        )
            @test_throws ArgumentError H._restore_checkpoint_state!(restored, target, corrupt)
            @test isequal(snapshot(restored, target), rejected)
        end
        mktempdir() do directory
            checkpoint = joinpath(directory, "impact.bin")
            archive = joinpath(directory, "impact.archive")
            save_checkpoint(checkpoint, st, ps)
            load_checkpoint!(restored, target, checkpoint)
            archive_run(archive, st, ps; rng_seed = 17)
            H._restore_checkpoint_state!(restored, target, load_archive(archive).state)
            @test snapshot(restored, target) == snapshot(st, ps)
        end
        empty_species = ParticleSet{1,T}[]
        empty_state = state(T, empty_species)
        step!(empty_state, empty_species, T(0.01))
        @test isempty(empty_state._stage[].event_candidates)
        collection = [population(T, T[0.25], T[-1])]
        collection_state = state(T, collection)
        step!(collection_state, collection, T(0.5))
        collection_saved = deepcopy(H._checkpoint_state(collection_state, collection))
        resized = [ParticleSet{1,T}(0; q = 0)]
        resumed = state(T, ParticleSet{1,T}[])
        H._restore_checkpoint_state!(resumed, resized, collection_saved)
        @test snapshot(resumed, resized) == snapshot(collection_state, collection)
        @test_throws ArgumentError H._restore_checkpoint_state!(resumed, resized, saved)
        step!(resumed, resized, T(0.5))
        step!(collection_state, collection, T(0.5))
        @test snapshot(resumed, resized) == snapshot(collection_state, collection)
        # Known writable leaves/Refs still require complete independence.
        alias_ps = population(T, T[0.5], T[0])
        alias_state = state(T, alias_ps)
        alias_ps.weight = alias_ps.x[1]
        aliased = snapshot(alias_state, alias_ps)
        @test_throws ArgumentError step!(alias_state, alias_ps, T(0.01))
        @test isequal(snapshot(alias_state, alias_ps), aliased)
        fields = ntuple(i -> getfield(st, i), fieldcount(typeof(st)))
        index = findfirst(==(:cumulative_wall_count), fieldnames(typeof(st)))
        ref_alias = typeof(st)(Base.setindex(fields, st.last_wall_count, index)...)
        ref_before = snapshot(ref_alias, ps)
        @test_throws ArgumentError step!(ref_alias, ps, T(0.5))
        @test isequal(snapshot(ref_alias, ps), ref_before)
        calls = (
            (init!, Tuple{typeof(st),typeof(ps)}),
            (step!, Tuple{typeof(st),typeof(ps),T}),
            (yee_impact_longitudinal_diagnostics, Tuple{typeof(st),typeof(ps)}),
            (H._restore_checkpoint_state!, Tuple{typeof(st),typeof(ps),typeof(saved)}),
        )
        for (f, types) in calls
            report = JET.report_call(f, types)
            if !isempty(JET.get_reports(report))
                show(stdout, MIME"text/plain"(), report)
                println()
            end
            @test isempty(JET.get_reports(report))
        end
    end
end
end
