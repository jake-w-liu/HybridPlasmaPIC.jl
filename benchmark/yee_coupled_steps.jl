# Remote observation through the existing yee-electromagnetic-pic manual suite.
# Reports complete public steps, allocation traffic and GC metadata, without
# performance thresholds. Initialization/checkpoints/logging are not timed.
module ClosedYeeStepCostObservation
using HybridPlasmaPIC, SHA
const H = HybridPlasmaPIC

bits(x::Float32) = reinterpret(UInt32, x)
bits(x::Float64) = reinterpret(UInt64, x)

function ordinary(Model, ::Type{T}, transverse) where {T}
    g = YeeGrid((8,), (one(T),))
    ps = ParticleSet{1,T}(2; q = T(0.3), m = T(1.2))
    ps.x[1] .= T[0.34, 0.63]
    ps.weight .= T[0.02, 0.03]
    ps.v[1] .= T[0.012, -0.017]
    config = (;
        shape = YeeTriangularCloud(T(0.19)),
        background_charge = T(-0.015),
        left_plate_charge = T(0.01),
    )
    st = Model === YeeElectromagneticPIC ? Model(g; config..., c = T(0.7)) : Model(g; config...)
    if transverse
        ps.v[2] .= T[0.04, -0.02]
        ps.v[3] .= T[-0.03, 0.025]
        for i = 1:9
            x = T(i - 1) / 8
            st.fields.E[2][i] = T(0.04) * sinpi(x)
            st.fields.E[3][i] = T(0.03) * sinpi(2x)
        end
        st.fields.E[2][end] = st.fields.E[3][end] = zero(T)
        st.fields.B[2] .= T(0.025)
        st.fields.B[3] .= T(-0.02)
    end
    init!(st, ps)
    return st, ps, T(0.002)
end

function range_case(Model, name)
    T = Float32
    # These are the existing directed original-force fixtures, not a size sweep.
    n, a, field, x, v, q, m, dt = if name === :underflow
        (
            2,
            ldexp(1.0, 100),
            nextfloat(zero(T)),
            T(0.25),
            zero(T),
            ldexp(one(T), 126),
            ldexp(one(T), -126),
            ldexp(one(T), -56),
        )
    elseif name === :overflow
        (
            7,
            ldexp(1.0, -10),
            floatmax(T),
            T(0.5),
            zero(T),
            ldexp(one(T), -126),
            ldexp(one(T), 126),
            ldexp(one(T), -20),
        )
    elseif name === :cancellation
        (7, ldexp(1.0, -10), one(T), T(0.5), -one(T), T(256), one(T), T(1 / 128))
    elseif name === :magnetic_overflow
        @assert Model === YeeElectromagneticPIC
        (
            7,
            ldexp(1.0, -10),
            zero(T),
            T(0.5),
            one(T),
            ldexp(one(T), -100),
            one(T),
            ldexp(one(T), -40),
        )
    else
        error("unknown directed range case: $name")
    end
    ps = ParticleSet{1,T}(1; q, m)
    ps.x[1][1], ps.v[1][1], ps.weight[1] = x, v, zero(T)
    st = Model(YeeGrid((n,), (one(T),)); shape = YeeTriangularCloud(a), left_plate_charge = field)
    if name === :magnetic_overflow
        ps.v[3][1] = one(T)
        st.fields.B[2] .= floatmax(T)
    end
    init!(st, ps)
    return st, ps, dt
end

function two_steps!(st, ps, dt)
    step!(st, ps, dt)
    step!(st, ps, dt)
    return st
end
sample(st::S, ps::P, dt::T) where {S,P,T} = @timed two_steps!(st, ps, dt)
source(ps) = (ps.q, ps.m, copy(ps.weight), copy(ps.id), copy(ps.tag))
storage(st, ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, st.fields.E..., st.fields.B...)

function observe(factory, name, bclass)
    warm, particles, dt = factory()
    step!(warm, particles, dt)
    sample(warm, particles, dt) # Compile this actual two-step wrapper; discard the observation.
    times, traffic = Float64[], Int[]
    expected = nothing
    for repetition = 1:3
        st, ps, dt = factory()
        step!(st, ps, dt) # All measured samples begin at the identical fixed-dt phase.
        @assert st.step[] == 1 && isequal(st.last_dt[], dt)
        before = H._checkpoint_state(st, ps)
        @assert before.schema_version === 2
        @assert before.force_policy === :original_geometry_boris_nearest_even
        @assert before.force_policy_version === 1
        original, owners = source(ps), storage(st, ps)
        start_time = st.time[]
        config = (st.g, st.shape, st.background_charge, st.left_plate_charge)
        GC.gc()
        measured = sample(st, ps, dt)
        # These checks and checkpoint-state comparisons are outside measurement.
        @assert measured.value === st && st.step[] == 3
        @assert isfinite(st.time[]) && st.time[] > start_time && isequal(st.last_dt[], dt)
        @assert isequal(source(ps), original)
        @assert isequal(config, (st.g, st.shape, st.background_charge, st.left_plate_charge))
        @assert all(pair -> pair[1] === pair[2], zip(owners, storage(st, ps)))
        @assert all(a -> all(isfinite, a), (ps.x..., ps.v..., st.fields.E..., st.fields.B...))
        saved = H._checkpoint_state(st, ps) # Also validates required finite histories.
        if repetition == 1
            expected = deepcopy(saved)
        else
            @assert isequal(saved, expected)
        end
        push!(times, measured.time)
        push!(traffic, measured.bytes)
        gc = measured.gcstats
        println((;
            purpose = :closed_yee_complete_step_cost,
            acceptance = false,
            profile = name,
            model = nameof(typeof(st)),
            T = eltype(ps.weight),
            layout = :scalar,
            ncell = st.g.n[1],
            particles = length(ps.weight),
            weight_class = all(iszero, ps.weight) ? :zero_weight_tracer : :positive_weight,
            initial_B = bclass,
            width_type = typeof(st.shape.half_width),
            width_bits = bits(st.shape.half_width),
            dt_type = typeof(dt),
            dt_bits = bits(dt),
            sample = repetition,
            measured_steps = 2,
            start_step = 1,
            end_step = st.step[],
            start_time,
            end_time = st.time[],
            seconds = measured.time,
            allocated_traffic_bytes = measured.bytes,
            gc_seconds = measured.gctime,
            compile_seconds = measured.compile_time,
            recompile_seconds = measured.recompile_time,
            lock_conflicts = measured.lock_conflicts,
            malloc = gc.malloc,
            realloc = gc.realloc,
            poolalloc = gc.poolalloc,
            bigalloc = gc.bigalloc,
            freecall = gc.freecall,
            gc_pauses = gc.pause,
            gc_full_sweeps = gc.full_sweep,
        ))
        flush(stdout)
    end
    # Exactly three raw samples are retained, including any compilation or GC.
    println((;
        purpose = :closed_yee_cost_summary,
        acceptance = false,
        profile = name,
        model = nameof(typeof(warm)),
        T = eltype(particles.weight),
        initial_B = bclass,
        samples = 3,
        measured_steps = 2,
        minimum_seconds = minimum(times),
        median_seconds = sort(times)[2],
        median_allocated_traffic_bytes = sort(traffic)[2],
    ))
    flush(stdout)
    return nothing
end

function run_observation()
    root = dirname(@__DIR__)
    revision = readchomp(Cmd(`git rev-parse HEAD`; dir = root))
    github_sha = get(ENV, "GITHUB_SHA", nothing)
    @assert github_sha === nothing || github_sha == revision
    println((;
        purpose = :closed_yee_complete_step_cost,
        acceptance = false,
        source_sha = revision,
        github_sha,
        script_sha256 = bytes2hex(sha256(read(@__FILE__))),
        julia = VERSION,
        julia_build_commit = Base.GIT_VERSION_INFO.commit,
        machine = Sys.MACHINE,
        architecture = Sys.ARCH,
        kernel = Sys.KERNEL,
        cpu_name = Sys.CPU_NAME,
        cpu_models = unique([cpu.model for cpu in Sys.cpu_info()]),
        logical_cpus = Sys.CPU_THREADS,
        default_threads = Threads.nthreads(:default),
        interactive_threads = Threads.nthreads(:interactive),
        gc_threads = Threads.ngcthreads(),
        julia_num_threads = get(ENV, "JULIA_NUM_THREADS", nothing),
        julia_num_gc_threads = get(ENV, "JULIA_NUM_GC_THREADS", nothing),
        rounding = rounding(Float64),
        flush_to_zero = get_zero_subnormals(),
        runner_os = get(ENV, "RUNNER_OS", nothing),
        runner_arch = get(ENV, "RUNNER_ARCH", nothing),
        suite = get(ENV, "IMPLEMENTATION_SUITE", nothing),
        phase = :warmed_matched_fixed_dt,
        profiles = 13,
        samples = 3,
        measured_steps = 2,
        unmeasured_steps = 1,
        bytes_meaning = :allocation_traffic,
        peak_memory = :unmeasured,
    ))
    flush(stdout)
    count = 0
    for T in (Float32, Float64),
        (Model, transverse) in
        ((YeeLongitudinalPIC, false), (YeeElectromagneticPIC, false), (YeeElectromagneticPIC, true))

        observe(() -> ordinary(Model, T, transverse), :ordinary, transverse ? :nonzero : :zero)
        count += 1
    end
    for Model in (YeeLongitudinalPIC, YeeElectromagneticPIC),
        name in (:underflow, :overflow, :cancellation)

        observe(() -> range_case(Model, name), name, :zero)
        count += 1
    end
    observe(
        () -> range_case(YeeElectromagneticPIC, :magnetic_overflow),
        :magnetic_overflow,
        :nonzero,
    )
    @assert count + 1 == 13
    println((;
        purpose = :closed_yee_complete_step_cost_complete,
        profiles = 13,
        raw_samples = 39,
        acceptance = false,
        peak_memory = :unmeasured,
    ))
    flush(stdout)
end
run_observation()
end
