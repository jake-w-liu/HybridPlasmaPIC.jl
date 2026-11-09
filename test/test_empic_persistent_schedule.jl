module EMPICPersistentScheduleTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(::Type{T}, Model; D = 1, cap = T[1/16, 5/32, 0], M = 3,
    budget = 64, adaptive = nothing, lower = nothing, limits = nothing,
    weighted = false, relativistic = false) where {T}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
    species = [ParticleSet{D,T}(1; q = q, m = m) for (q, m) in ((-2,3),(2,5),(0,1))]
    for (i, p) in enumerate(species)
        for d in 1:D
            p.x[d][1] = T(1/4 + (i-1)/16)
        end
        p.weight[1] = weighted && i < 3 ? one(T) : zero(T)
        i == 3 && continue
        p.v[1][1], p.v[2][1] = T(i/64), T(-i/128)
    end
    es = Model(g, species; mobile = Bool[true,true,false], species_dt = cap,
        n_sub = M, event_budget = budget, species_adaptive = adaptive,
        species_dt_min = lower, schedule_limits = limits, shape = CIC(), relativistic)
    if es isa EMPIC
        for (a, value) in zip(es.E, (T(1/8),T(-1/16),T(1/32)))
            fill!(a, value)
        end
        fill!(es.B[3], T(1/8))
    else
        fill!(es.Ex, T(1/8)); fill!(es.Ey, T(-1/16)); fill!(es.Bz, T(1/8))
    end
    init_empic!(es, species)
    return es, species
end

phase_value(p) = (time = p.time, primed = p.primed, interval = p.interval,
    last = p.last, next = p.next, report = p.report)
function completed(es, species)
    s = es.collection.schedule
    return deepcopy((fields = H._empic_event_fields(es),
        species = [H._checkpoint_particle_state(p) for p in species],
        clock = (es.time[],es.step[],es.last_dt[],es.last_sub_dt[]),
        proposal = s.proposal, phase = phase_value(s.phase),
        initialized = es.collection.initialized[], n_sub = es.n_sub,
        c = es.c, relativistic = es.relativistic))
end
setfuture(es; kw...) = es isa EMPIC ? set_empic_schedule!(es; kw...) : set_empic1d_schedule!(es; kw...)

cross3(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
function boris(v, alpha, h, E, B)
    u = v .+ (alpha*h/2) .* E
    t = (alpha*h/2) .* B
    s = (2/(1+sum(abs2,t))) .* t
    return u .+ cross3(u .+ cross3(u,t),s) .+ (alpha*h/2) .* E
end

# Independent rational calendar and uniform-field Boris reference. Zero weights
# remove all charge/current feedback while retaining nonzero force and drift.
function reference(species, D)
    return (time = Ref(zero(Q)), primed = Ref(false), H = zeros(eltype(first(species)),length(species)),
        last = fill(zero(Q),length(species)), next = fill(zero(Q),length(species)),
        x = [Tuple(Q(p.x[d][1]) for d in 1:D) for p in species],
        v = [Tuple(Q(p.v[c][1]) for c in 1:3) for p in species])
end
function reference_step!(r, species, duration::T, cap, M, E, B) where {T}
    origin, finish = r.time[], r.time[] + Q(duration)
    first = !r.primed[]
    if first
        for i in 1:2
            r.H[i] = cap[i]
            alpha, v = Q(species[i].q)/Q(species[i].m), r.v[i]
            r.v[i] = v .- (Q(cap[i])/2*alpha) .* (E .+ cross3(v,B))
            r.last[i], r.next[i] = origin, origin + Q(cap[i])
        end
    end
    edges, kicks, base = Q[origin], Tuple{Int,Q}[], 1
    time = origin
    while time < finish
        for i in 1:2
            if first || r.next[i] == time
                selected = cap[i]
                consumed = r.H[i]/T(2) + selected/T(2)
                r.v[i] = boris(r.v[i], Q(species[i].q)/Q(species[i].m), Q(consumed), E, B)
                r.H[i], r.last[i], r.next[i] = selected, time, time + Q(selected)
                push!(kicks,(i,time))
            end
        end
        while base <= M && origin + Q(duration)*base/M <= time
            base += 1
        end
        edge = min(finish, r.next[1], r.next[2], origin + Q(duration)*base/M)
        actual = Q(T(edge-time))
        for i in 1:2
            r.x[i] = map((x,v) -> mod(x+actual*v,one(Q)), r.x[i],r.v[i][1:length(r.x[i])])
        end
        push!(edges,edge)
        time, first = edge, false
    end
    r.time[], r.primed[] = finish, true
    return edges, kicks
end

@testset "Persistent nondividing coupled calendar and force reference" begin
    for T in (Float32,Float64), Model in (EMPIC,EMPIC1D), D in 1:3
        (D > 1 && (T === Float32 || Model === EMPIC1D)) && continue
        es, species = fixture(T,Model; D)
        r = reference(species,D)
        E = Model === EMPIC ? Q.(map(first,es.E)) : (Q(first(es.Ex)),Q(first(es.Ey)),zero(Q))
        B = Model === EMPIC ? Q.(map(first,es.B)) : (zero(Q),zero(Q),Q(first(es.Bz)))
        for (duration,M,cap) in ((T(3/32),3,T[1/16,5/32,0]),
            (T(1/16),2,T[3/32,3/16,0]), (T(1/32),1,T[3/32,3/16,0]))
            pending = deepcopy(es.collection.schedule.phase.next)
            setfuture(es; n_sub=M,species_dt=cap)
            @test es.collection.schedule.phase.next == pending
            edges,kicks = reference_step!(r,species,duration,cap,M,E,B)
            step_empic!(es,species,duration)
            p = es.collection.schedule.phase
            @test p.time == r.time[] && p.next == r.next && p.last == r.last
            @test p.report.edges == edges
            @test [(k.species,k.time) for k in p.report.kicks] == kicks
            for i in 1:2, c in 1:3
                @test isapprox(species[i].v[c][1],T(r.v[i][c]); rtol=256eps(T),atol=256eps(T))
            end
            for i in 1:2, d in 1:D
                @test isapprox(species[i].x[d][1],T(r.x[i][d]); rtol=256eps(T),atol=256eps(T))
            end
            d = electromagnetic_diagnostics(es,species)
            @test d.velocity_time_offsets[1:2] == [(p.last[i]+p.next[i])/2-p.time for i in 1:2]
            @test iszero(d.velocity_time_offsets[3])
            @test d.public_clock_defect == Q(es.time[])-p.time
            before = completed(es,species)
            step_empic!(es,species,zero(T))
            @test isequal(completed(es,species),before)
        end
    end
end

@testset "Persistent admission, fractional caps, empty rows and tiny event gaps" begin
    for T in (Float32,Float64)
        for bad in (T[-1,1,0],T[1,1,1])
            @test_throws ArgumentError fixture(T,EMPIC; cap=bad)
        end
        @test_throws ArgumentError fixture(T,EMPIC; budget=true)
        @test_throws ArgumentError fixture(T,EMPIC; cap=[Q(1,BigInt(2)^2000),Q(1),Q(0)])
        es,species = fixture(T,EMPIC; cap=T[0.3*(1/32),5/32,0],M=2)
        step_empic!(es,species,T(1/32))
        @test es.collection.schedule.phase.time == Q(T(1/32))
        @test es.collection.schedule.phase.report.kicks[1].selected == T(0.3*(1/32))
        tiny = T(2)*floatmin(T)
        blocked,populations = fixture(T,EMPIC; cap=T[tiny,nextfloat(tiny),0],M=1)
        before = completed(blocked,populations)
        @test_throws ArgumentError step_empic!(blocked,populations,T(2)*tiny)
        @test isequal(completed(blocked,populations),before)
        for Model in (EMPIC,EMPIC1D)
            empty = ParticleSet{1,T}[]
            vacuum = Model(FourierGrid((1,),(one(T),)),empty;
                species_dt=T[],event_budget=8,n_sub=3)
            init_empic!(vacuum,empty)
            step_empic!(vacuum,empty,T(1/32))
            @test isempty(vacuum.collection.schedule.phase.interval)
            @test length(vacuum.collection.schedule.phase.report.intervals) == 3
        end
    end
end

@testset "Dividing old owner remains the all-fast kernel path" begin
    for Model in (EMPIC,EMPIC1D)
        event, particles = fixture(Float64,Model; cap=[1/128,1/128,0],M=4)
        oldparticles = deepcopy(particles)
        old = Model(FourierGrid(event.g.n,event.g.L),oldparticles;
            mobile=Bool[true,true,false],n_sub=4,shape=CIC())
        foreach(copyto!,H._empic_collection_seed_fields(old),H._empic_collection_seed_fields(event))
        init_empic!(old,oldparticles)
        step_empic!(event,particles,1/32)
        step_empic!(old,oldparticles,1/32)
        @test isequal(H._empic_physical_fields(event),H._empic_physical_fields(old))
        @test [H._checkpoint_particle_state(p) for p in particles] ==
            [H._checkpoint_particle_state(p) for p in oldparticles]
    end
end

@testset "Sampled adaptive decisions and transactional failure/retry" begin
    limits = (rotation=1/64,electric=1/64,travel=1/64,growth=2.0)
    es,species = fixture(Float64,EMPIC; cap=[1/4,5/32,0],adaptive=Bool[true,false,false],
        lower=[1/1024,0,0],limits)
    full_v = deepcopy(species[1].v)
    step_empic!(es,species,1/32)
    k = first(es.collection.schedule.phase.report.kicks)
    @test k.bootstrap && k.V1 == sum(Q(abs(v[1])) for v in full_v)
    @test k.selected < k.upper && k.halvings > 0
    selected = k.upper
    while !(Q(selected)*abs(Q(k.q))*k.B1 <= Q(limits.rotation)*Q(k.m) &&
        Q(selected)*abs(Q(k.q))*k.E1 <= Q(limits.electric)*Q(k.m)*Q(es.c) &&
        Q(selected)*k.V1 <= Q(limits.travel)*minimum(Q,es.g.dx))
        selected /= 2
    end
    @test selected == k.selected
    before = completed(es,species)
    @test_throws ArgumentError setfuture(es; species_dt=[-1,1,0])
    @test isequal(completed(es,species),before)
    setfuture(es; event_budget=1,n_sub=1)
    limited = completed(es,species)
    @test_throws ArgumentError step_empic!(es,species,1.0)
    @test isequal(completed(es,species),limited)
    setfuture(es; event_budget=128)
    step_empic!(es,species,1/64)
    @test es.step[] == 2
    committed = completed(es,species)
    species[1].weight[1] = species[2].weight[1] = 1.0
    @test_throws ArgumentError step_empic!(es,species,0.0)
    species[1].weight[1] = committed.species[1].weight[1]
    species[2].weight[1] = committed.species[2].weight[1]
    @test isequal(completed(es,species),committed)
end

@testset "Every mobile row supplies the final actual midpoint current" begin
    es,species = fixture(Float64,EMPIC; weighted=true,cap=[1/8,5/32,0],M=3)
    step_empic!(es,species,1/32)
    @test all(k -> k.bootstrap,es.collection.schedule.phase.report.kicks)
    @test length(es.collection.schedule.phase.report.intervals) == 3
    for c in 1:3
        expected = zeros(es.g.n)
        for i in 1:2
            midpoint = deepcopy(species[i])
            foreach(copyto!,midpoint.x,es.collection.midpoints[i])
            contribution = similar(expected)
            deposit_scalar!(contribution,midpoint,
                midpoint.q .* midpoint.weight .* midpoint.v[c],es.g,CIC())
            expected .+= contribution
        end
        expected ./= prod(es.g.dx)
        @test isapprox(es.Jraw[c],expected; rtol=128eps(Float64),atol=128eps(Float64))
    end
    report = es.collection.schedule.phase.report
    @test report.charge_source_scale == sum(abs(Q(p.q))*sum(Q,p.weight) for p in species)
    for c in 1:3, I in eachindex(es.J[c])
        @test report.impulse[c][I] == sum(Q(h)*Q(s.current[c][I]) for (h,s) in zip(report.intervals,report.segments))
        @test report.raw_impulse[c][I] == sum(Q(h)*Q(s.raw_current[c][I]) for (h,s) in zip(report.intervals,report.segments))
    end
    @test report.field_defect == report.field_after-report.field_before+report.work_used
    @test charge_conservation_residual(es,es.last_sub_dt[]) >= 0
end

function invented_duration(saved)
    # Keep source edges/fields fixed, but invent a different consumed h and
    # consistently adjust every dependent observation. Only canonical gap
    # materialization admission distinguishes this from an actual receipt.
    phase, report = saved.event_phase, saved.event_phase.report
    intervals, defects = copy(report.intervals),copy(report.duration_defects)
    intervals[1] = nextfloat(intervals[1])
    defects[1] = Q(intervals[1])-(report.edges[2]-report.edges[1])
    clock = report.public_before
    for h in intervals
        clock += h
    end
    impulse = map(a -> fill(zero(Q),size(a)),report.impulse)
    raw = report.raw_impulse === nothing ? nothing : map(a -> fill(zero(Q),size(a)),report.raw_impulse)
    used_work,raw_work = zero(Q), raw === nothing ? nothing : zero(Q)
    volume = prod(Q, ntuple(d -> saved.L[d]/saved.ncell[d],saved.D))
    for (h,s) in zip(intervals,report.segments), c in eachindex(impulse), I in eachindex(impulse[c])
        impulse[c][I] += Q(h)*Q(s.current[c][I])
        Ebar = (Q(s.E_before[c][I])+Q(s.E_after[c][I]))/2
        used_work += volume*Q(h)*Q(s.current[c][I])*Ebar
        if raw !== nothing
            raw[c][I] += Q(h)*Q(s.raw_current[c][I])
            raw_work += volume*Q(h)*Q(s.raw_current[c][I])*Ebar
        end
    end
    altered = merge(report,(intervals=intervals,duration_defects=defects,
        public_after=clock,impulse=impulse,raw_impulse=raw,work_used=used_work,work_raw=raw_work,
        field_defect=report.field_after-report.field_before+used_work,
        clock_defect=Q(clock)-Q(report.public_before)-(report.finish-report.start)))
    return merge(saved,(time=clock,last_sub_dt=last(intervals),event_phase=merge(phase,(report=altered,))))
end

@testset "Schema3 raw/archive continuation and corruption rejection" begin
    for Model in (EMPIC,EMPIC1D)
        es,species = fixture(Float64,Model; relativistic=true)
        step_empic!(es,species,3/32)
        setfuture(es; species_dt=[3/32,3/16,0],n_sub=2)
        saved = H._checkpoint_state(es,species)
        @test saved.schema_version === 3 && saved.schedule_policy === 2
        @test saved.phase_policy === :persistent_event
        mktempdir() do dir
            raw,archive = joinpath(dir,"event.chk"),joinpath(dir,"event.run")
            save_checkpoint(raw,es,species)
            archive_run(archive,es,species; rng_seed=19)
            target,destination = fixture(Float64,Model)
            # All overwritten old physical values and phase/proposal values may
            # be stale. Bound ordinary roles and constructor policy remain real.
            for a in H._empic_transaction_fields(target)
                fill!(a,NaN)
            end
            target.time[],target.last_dt[],target.last_sub_dt[] = NaN,NaN,NaN
            target.collection.schedule.phase.interval .= NaN
            target.collection.schedule.proposal.cap .= NaN
            load_checkpoint!(target,destination,raw)
            @test isequal(completed(target,destination),completed(es,species))
            archived = load_archive(archive).state
            @test archived.schema_version === 3
            load_archive!(target,destination,archive)
            for duration in (0.0,1/16,1/32)
                step_empic!(es,species,duration)
                step_empic!(target,destination,duration)
                @test isequal(completed(target,destination),completed(es,species))
            end
            for damage in (:deadline,:impulse,:duration,:source_alias,:unknown)
                incoming = deepcopy(saved)
                if damage === :deadline
                    incoming.event_phase.next[1] += Q(1/128)
                elseif damage === :impulse
                    incoming.event_phase.report.impulse[1][1] += one(Q)
                elseif damage === :duration
                    incoming = invented_duration(incoming)
                elseif damage === :source_alias
                    incoming = merge(incoming,(event_proposal=merge(incoming.event_proposal,
                        (cap=target.collection.schedule.proposal.cap,)),))
                else
                    incoming = merge(incoming,(schedule_policy=7,))
                end
                before = completed(target,destination)
                @test_throws ArgumentError H._restore_checkpoint_state!(target,destination,incoming)
                @test isequal(completed(target,destination),before)
            end
        end
    end
end
end # module
