module TestYeeAbsorbingLongitudinalPIC
using HybridPlasmaPIC, Test, Serialization
import JET
include("yee_charge_oracle.jl")
const O = YeeChargeOracle
const H = HybridPlasmaPIC

function particles(::Type{T}, x, v; q = zero(T), m = one(T), weights = ones(T, length(x))) where {T}
    ps = ParticleSet{1,T}(length(x); q, m)
    ps.x[1] .= x
    ps.v[1] .= v
    ps.weight .= weights
    return ps
end
function state(::Type{T}, species; sigma = zero(T), width = Float64(1 / 16)) where {T}
    st = YeeAbsorbingLongitudinalPIC(
        YeeGrid((4,), (one(T),));
        shape = YeeTriangularCloud(width),
        left_plate_charge = sigma,
        max_root_refinements = 128,
    )
    init!(st, species)
    return st
end
function snapshot(st, input)
    species = input isa ParticleSet ? (input,) : input
    history = H._ya_history_values(st)
    if history.last_capture !== nothing
        record = history.last_capture
        names = fieldnames(typeof(record))
        capture = NamedTuple{names}(ntuple(i -> getfield(record, i), fieldcount(typeof(record))))
        history = merge(history, (last_capture = capture,))
    end
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
        history,
    ))
end
fieldenergy(E, dx) = O.rat(dx) / 2 * sum(abs2, O.rat.(E))
function death_reference(charge, side, n, L, a)
    nodes = O.cloud_nodes(side == 1 ? zero(L) : L, n, L, a)
    return [O.rat(charge) * ((side == 1 ? 0 : 1) - (1 - sum(nodes[1:e+1]))) for e = 0:n-1]
end

@testset "Independent finite-tail death continuity" begin
    for T in (Float32, Float64), a in (Float64(1 / 16), 2.0), side in (1, 2), charge in (-1, 1)
        g = YeeGrid((4,), (one(T),))
        raw = H._yee_cloud_death_sums(O.rat(charge), side, g, YeeTriangularCloud(a))
        expected = death_reference(charge, side, 4, one(T), a)
        @test raw[1:4] == expected
        nodes = O.rat(charge) .* O.cloud_nodes(side == 1 ? zero(T) : one(T), 4, one(T), a)
        transfer = -nodes
        transfer[side == 1 ? 1 : 5] += charge
        @test transfer + O.incidence(expected) == zeros(O.Q, 5)
        @test any(!iszero, expected)
    end
end

@testset "Stable capture and exact passive ledgers" begin
    for T in (Float32, Float64), side in (1, 2)
        ps = particles(T, side == 1 ? T[0.25, 0.6] : T[0.75, 0.4], side == 1 ? T[-1, 0] : T[1, 0])
        ps.tag .= UInt32[17, 29]
        st = state(T, ps)
        identities = H._yl_particle_arrays(ps)
        step!(st, ps, T(0.5))
        @test length(ps.weight) == 1 && ps.id == UInt64[2] && ps.tag == UInt32[29]
        @test all(pair -> pair[1] === pair[2], zip(identities, H._yl_particle_arrays(ps)))
        @test st.last_wall_count[side][] == 1 && sum(x -> x[], st.cumulative_wall_count) == 1
        @test st.last_wall_impulse[side][] == (side == 1 ? -1 : 1)
        @test st.last_wall_kinetic[side][] == st.last_wall_energy[side][] == 1 // 2
        @test iszero(st.last_wall_work[side][]) && iszero(st.last_total_energy_defect[])
        @test st.last_capture[].elapsed == 1 // 4 && st.last_capture[].external_step == 1
        @test st._initial_counts[] == [2]
        binding = st._binding[]
        step!(st, ps, T(0.5))
        @test st._binding[] === binding && st.last_capture[] === nothing
        @test sum(x -> x[], st.cumulative_wall_count) == 1
        @test all(iszero, yee_absorbing_longitudinal_diagnostics(st, ps).gauss_residual)
    end
    for T in (Float32, Float64)
        ps = particles(T, T[1/128, 1/8], T[-1, 0]; q = one(T), weights = T[1, 0])
        initial = deepcopy(ps)
        st = state(T, ps; sigma = T(-2))
        E0, K0 = copy(st.fields.E[1]), O.rat(kinetic_energy(ps))
        @test K0 > 0
        step!(st, ps, T(1 / 64))
        @test ps.id == UInt64[2] && ps.weight == T[0] && ps.v[1][1] < 0
        @test st.free_left_plate_charge[] == -1 && st.right_plate_charge[] == 1
        @test st.left_plate_charge == -2 && st.last_wall_charge[1][] == 1
        prex = T[0, st._stage[].event_candidates[1].x[1][1]]
        Ipre =
            T.(O.cloud_impulse([initial], [initial.x[1]], [prex], 4, one(T), st.shape.half_width))
        Ideath = T.(death_reference(1, 1, 4, one(T), st.shape.half_width))
        Epre, Epost = st._event_fields.E[1], st._death_fields.E[1]
        @test Epre[1] != Epost[1]
        @test all(O.rounded(e, O.rat(old) - O.rat(i)) for (e, old, i) in zip(Epre, E0, Ipre))
        @test all(O.rounded(e, O.rat(old) - O.rat(i)) for (e, old, i) in zip(Epost, Epre, Ideath))
        Wpre =
            O.rat(st.g.dx[1]) / 2 *
            sum(O.rat(i) * (O.rat(a) + O.rat(b)) for (i, a, b) in zip(Ipre, E0, Epre))
        Wdeath =
            O.rat(st.g.dx[1]) / 2 *
            sum(O.rat(i) * (O.rat(a) + O.rat(b)) for (i, a, b) in zip(Ideath, Epre, Epost))
        @test !iszero(Wdeath) && st.last_wall_work[1][] == Wdeath
        @test st.last_transport_work[] == Wpre
        @test st.last_wall_energy[1][] > 0 && st.last_wall_energy[1][] >= st.last_passivity_error[]
        @test st.last_wall_energy[1][] == st.last_wall_kinetic[1][] + Wdeath
        @test all(
            O.rounded(i, O.rat(a) + O.rat(b)) for (i, a, b) in zip(st.last_impulse[1], Ipre, Ideath)
        )
        @test O.rounded(st.last_energy_defect[], -K0 + st.last_wall_kinetic[1][] - Wpre)
        @test O.rounded(
            st.last_total_energy_defect[],
            -K0 + fieldenergy(st.fields.E[1], st.g.dx[1]) - fieldenergy(E0, st.g.dx[1]) +
            st.last_wall_energy[1][],
        )
    end
end

# Independent analytic-volume oracle for one q=m=w=1 capture in edge0 and
# one charged weight-zero tracer whose cloud stays strictly inside that edge.
# Simpson64/128 is a bounded reference check, not a BigFloat-epsilon campaign.
function capture_reference(initial_field, subdivisions)
    a, h, x0, terminal = BigFloat(1 / 16), BigFloat(1 / 4), BigFloat(1 / 128), BigFloat(1 / 64)
    cumulative(x) = (x + (a - x)^3 / (6a^2)) / h
    # The semidiscrete reference starts from the actual stored initial field,
    # retaining its Gauss-rounding residual rather than projecting it away.
    invariant = BigFloat(O.rat(initial_field)) + cumulative(x0)
    electric(x) = invariant - cumulative(x)
    potential(x) = h / 2 * (electric(x)^2 + 3)
    speed(x) = sqrt(1 + 2 * (potential(x0) - potential(x)))
    function integral(f)
        dx = x0 / subdivisions
        total = f(zero(x0)) + f(x0)
        for j = 1:subdivisions-1
            total += (isodd(j) ? 4 : 2) * f(dx * j)
        end
        return dx * total / 3
    end
    hit = integral(x -> inv(speed(x)))
    @test 0 < hit < terminal
    return (
        hit = hit,
        incoming = -speed(zero(x0)),
        survivor = integral(x -> electric(x) / speed(x)) + invariant * (terminal - hit),
    )
end
@testset "Independent nonzero force-jump refinement" begin
    setprecision(BigFloat, 192) do
        seed = particles(Float64, [1 / 128, 1 / 8], [-1.0, 0]; q = 1.0, weights = [1.0, 0])
        initial = state(Float64, seed; sigma = -2.0)
        coarse = capture_reference(initial.fields.E[1][1], 64)
        fine = capture_reference(initial.fields.E[1][1], 128)
        oracle_error =
            abs(coarse.hit - fine.hit) +
            abs(coarse.incoming - fine.incoming) +
            abs(coarse.survivor - fine.survivor)
        errors = BigFloat[]
        for steps in (4, 8, 16)
            ps = particles(Float64, [1 / 128, 1 / 8], [-1.0, 0]; q = 1.0, weights = [1.0, 0])
            st = state(Float64, ps; sigma = -2.0)
            dt, captured = (1 / 64) / steps, nothing
            for _ = 1:steps
                step!(st, ps, dt)
                st.last_capture[] === nothing || (captured = st.last_capture[])
            end
            @test captured !== nothing && sum(x -> x[], st.cumulative_wall_count) == 1
            hit = O.rat(dt) * (captured.external_step - 1) + captured.elapsed
            push!(
                errors,
                abs(BigFloat(hit) - fine.hit) +
                abs(BigFloat(captured.velocity) - fine.incoming) +
                abs(BigFloat(ps.v[1][1]) - fine.survivor),
            )
        end
        orders = (log2(errors[1] / errors[2]), log2(errors[2] / errors[3]))
        @test oracle_error < minimum(errors) / 32
        @test minimum(errors) > 1024 * eps(Float64)
        @test all(>(1.8), orders)
        @info "Capture force-jump refinement" errors orders oracle_error
    end
end

@testset "Passive/late-compaction failure and reusable binding" begin
    for T in (Float32, Float64)
        terminal = particles(T, T[0.25], T[-1])
        terminal_state = state(T, terminal)
        terminal_before = snapshot(terminal_state, terminal)
        @test_throws ArgumentError step!(terminal_state, terminal, T(0.25))
        @test isequal(snapshot(terminal_state, terminal), terminal_before)
        simultaneous = particles(T, T[0.25, 0.75], T[-1, 1])
        simultaneous_state = YeeAbsorbingLongitudinalPIC(
            YeeGrid((4,), (one(T),));
            shape = YeeTriangularCloud(T(1 / 16)),
            max_root_refinements = 4,
        )
        init!(simultaneous_state, simultaneous)
        simultaneous_before = snapshot(simultaneous_state, simultaneous)
        @test_throws ArgumentError step!(simultaneous_state, simultaneous, T(0.5))
        @test isequal(snapshot(simultaneous_state, simultaneous), simultaneous_before)
        ps = particles(T, T[1/4096], T[-0.25]; q = one(T))
        st = state(T, ps; sigma = T(16))
        before, binding = snapshot(st, ps), st._binding[]
        @test_throws ArgumentError step!(st, ps, T(1 / 512))
        @test isequal(snapshot(st, ps), before) && st._binding[] === binding
        @test any(!iszero, st._death_fields.E[1])
        eta = reinterpret(T, T === Float32 ? UInt32(1) : UInt64(1))
        x = ldexp(one(T), exponent(eta) + exponent(floatmax(T)) - 1)
        captured = particles(T, T[x], T[-floatmax(T)]; m = floatmin(T), weights = T[floatmin(T)])
        fast = prevfloat(prevfloat(prevfloat(prevfloat(floatmax(T)))))
        survivor =
            particles(T, T[0.5], T[fast]; q = inv(floatmin(T)), m = floatmin(T), weights = T[0])
        species = [captured, survivor]
        failed = state(T, species; sigma = T(16), width = T(1 / 8))
        original, committed = snapshot(failed, species), failed._binding[]
        @test kinetic_energy(captured) > 0
        @test_throws ArgumentError step!(failed, species, eta)
        @test isempty(failed._stage[].event_candidates[1].weight)
        @test length(failed._binding[].candidates[1].weight) == 1
        @test isequal(snapshot(failed, species), original) && failed._binding[] === committed
        survivor.v[1][1] = zero(T)
        step!(failed, species, eta)
        @test isempty(captured.weight) && failed.step[] == 1 && all(isfinite, survivor.v[1])
        collection_saved = deepcopy(H._checkpoint_state(failed, species))
        targets = [ParticleSet{1,T}(0), ParticleSet{1,T}(0)]
        restored = state(T, targets; sigma = T(16), width = T(1 / 8))
        H._restore_checkpoint_state!(restored, targets, collection_saved)
        @test snapshot(restored, targets) == snapshot(failed, species)
        @test isempty(targets[1].weight) && length(targets[2].weight) == 1
    end
end

@testset "Strict detached capture restart" begin
    for T in (Float32, Float64)
        ps = particles(T, T[1/128, 1/8], T[-1, 0]; q = one(T), weights = T[1, 0])
        st = state(T, ps; sigma = T(-2))
        step!(st, ps, T(1 / 64))
        saved = deepcopy(H._checkpoint_state(st, ps))
        target = ParticleSet{1,T}(0)
        resumed = state(T, ParticleSet{1,T}(0); sigma = T(-2))
        resumed.free_left_plate_charge[] = O.rat(-999)
        H._restore_checkpoint_state!(resumed, target, saved)
        @test snapshot(resumed, target) == snapshot(st, ps)
        for corrupt in (
            merge(saved, (schema_version = 2,)),
            merge(saved, (passivity_policy_version = true,)),
            merge(saved, (free_left_plate_charge = saved.free_left_plate_charge + 1,)),
            merge(saved, (initial_counts = saved.initial_counts .+ 1,)),
            merge(saved, (last_wall_energy = (-one(O.Q), zero(O.Q)),)),
            merge(saved, (shape_width_type = Float32,)),
        )
            before = snapshot(resumed, target)
            @test_throws ArgumentError H._restore_checkpoint_state!(resumed, target, corrupt)
            @test isequal(snapshot(resumed, target), before)
        end
        view_target = ParticleSet{1,T}(
            (view(zeros(T, 2), 1:2:2),),
            ntuple(_ -> view(zeros(T, 2), 1:2:2), 3),
            view(zeros(T, 2), 1:2:2),
            view(UInt64[1, 0], 1:2:2),
            view(UInt32[0, 0], 1:2:2),
            zero(T),
            one(T),
        )
        before = snapshot(resumed, view_target)
        @test_throws ArgumentError H._restore_checkpoint_state!(resumed, view_target, saved)
        @test isequal(snapshot(resumed, view_target), before)
        resumed.fields.E[2] .= -zero(T)
        swapped = merge(
            saved,
            (
                fields = (
                    E = (saved.fields.E[1], resumed.fields.E[3], resumed.fields.E[2]),
                    B = saved.fields.B,
                ),
            ),
        )
        H._restore_checkpoint_state!(resumed, target, swapped)
        @test all(x -> isequal(x, -zero(T)), resumed.fields.E[3])
        H._restore_checkpoint_state!(resumed, target, saved)
        mktempdir() do directory
            checkpoint, archive =
                joinpath(directory, "capture.bin"), joinpath(directory, "capture.archive")
            save_checkpoint(checkpoint, st, ps)
            load_checkpoint!(resumed, target, checkpoint)
            archive_run(archive, st, ps; rng_seed = 19)
            H._restore_checkpoint_state!(resumed, target, load_archive(archive).state)
        end
        step!(st, ps, T(1 / 64))
        step!(resumed, target, T(1 / 64))
        @test snapshot(resumed, target) == snapshot(st, ps)
        empty = ParticleSet{1,T}[]
        empty_state = state(T, empty)
        step!(empty_state, empty, T(1 / 64))
        @test isempty(empty_state._initial_counts[]) && empty_state.step[] == 1
        @test_throws ArgumentError H._restore_checkpoint_state!(empty_state, empty, saved)
        calls = (
            (init!, Tuple{typeof(st),typeof(ps)}),
            (step!, Tuple{typeof(st),typeof(ps),T}),
            (yee_absorbing_longitudinal_diagnostics, Tuple{typeof(st),typeof(ps)}),
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
