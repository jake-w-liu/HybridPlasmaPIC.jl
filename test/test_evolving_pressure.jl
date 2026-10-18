using HybridPlasmaPIC, Test, LinearAlgebra

function ep_fields(f)
    return (
        f.n,
        f.ui...,
        f.B...,
        f.E...,
        f.J...,
        f.pe,
        f.gradp...,
        f.ninv,
        f.lapJ...,
        f.pforce...,
        f.cgl_dp,
    )
end

function ep_load(::Type{T}; γ = T(5 / 3), compressive = false, collection = false) where {T}
    g = FourierGrid((32,), (T(2pi),))
    ps = ParticleSet{1,T}(128)
    load_lattice_1d!(ps, zero(T), g.L[1])
    set_density_weight!(ps, one(T), g)
    if compressive
        @. ps.v[1] = sin(ps.x[1])
    else
        fill!(ps.v[1], T(0.2))
    end
    st = HybridStepper(g, HybridModel(EvolvingElectrons(γ)), CIC(), nparticles(ps))
    x = T[(i - 1) * g.dx[1] for i = 1:g.n[1]]
    st.fields.pe .= compressive ? ones(T, g.n) : (@. T(0.6) + T(0.08) * cos(x))
    populations = collection ? [ps] : ps
    init!(st, populations)
    return st, populations
end

@testset "Evolving pressure configuration and state ownership" begin
    for T in (Float32, Float64)
        @test closure_gamma(EvolvingElectrons(T(1.4))) == T(1.4)
        for bad in (zero(T), one(T), -one(T), T(Inf), T(NaN))
            @test_throws ArgumentError EvolvingElectrons(bad)
            @test_throws ArgumentError EvolvingElectrons{T}(bad)
        end
        c = EvolvingElectrons(T(5 / 3))
        for keyword in (:η, :ηH, :de2)
            @test_throws ArgumentError HybridModel(c; NamedTuple{(keyword,)}((0.01,))...)
        end
        g = FourierGrid((8,), (T(2pi),))
        @test_throws ArgumentError CAMCLStepper(g, HybridModel(c), CIC(), 0)
        @test_throws ArgumentError HallMHDModel(c)
        @test_throws ArgumentError electron_pressure!(zeros(T, 8), ones(T, 8), c)
        @test HybridStepper(g, HybridModel(IsothermalElectrons(T(0.2))), CIC(), 0).pressure_work ===
              nothing
        for collection in (false, true)
            st, ps = ep_load(T; collection)
            before = copy(st.fields.pe)
            init!(st, ps)
            @test st.fields.pe == before
            ohms_law!(st.fields, st.model, st.g)
            @test st.fields.pe == before
            ei = electron_internal_energy(st.fields.n, c, st.g; pe = st.fields.pe)
            @test ei ≈ sum(before) / (T(5 / 3) - one(T)) * prod(st.g.dx)
            @test_throws ArgumentError electron_internal_energy(st.fields.n, c, st.g)
            @test_throws ArgumentError energy_budget(ps, st.fields.B, st.fields.n, c, st.g)
            b = energy_budget(ps, st.fields.B, st.fields.n, c, st.g; pe = st.fields.pe)
            @test b.electron_internal == ei
            st.fields.pe .*= 2
            @test energy_budget(
                ps,
                st.fields.B,
                st.fields.n,
                c,
                st.g;
                pe = st.fields.pe,
            ).electron_internal ≈ 2ei
            copyto!(st.fields.pe, before)
            @test_throws ArgumentError electron_internal_energy(
                st.fields.n,
                IsothermalElectrons(T(0.2)),
                st.g;
                pe = before,
            )
            @test_throws DimensionMismatch electron_internal_energy(
                st.fields.n,
                c,
                st.g;
                pe = before[1:2],
            )
            for (dt, nb) in ((T(0.01), 1), (T(0.02), 2), (T(0.015), 4))
                step!(st, ps, dt; NB = nb)
                @test st.last_dt[] == dt
                @test all(p -> isfinite(p) && p >= zero(T), st.fields.pe)
            end
            @test maximum(abs, st.fields.pe .- before) > T(1e-5)
        end
        empty = ParticleSet{1,T}(0)
        st = HybridStepper(g, HybridModel(c), CIC(), 0)
        init!(st, empty)
        step!(st, empty, T(0.01); NB = 2)
        @test all(iszero, st.fields.pe)
        @test st.fields.floor_count[] == 8
    end
end

@testset "Evolving pressure manufactured stage equations" begin
    for T in (Float32, Float64)
        g = FourierGrid((32, 16), (T(2pi), T(2pi)))
        st = HybridStepper(g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
        x = reshape(T[(i - 1) * g.dx[1] for i = 1:32], :, 1)
        y = reshape(T[(j - 1) * g.dx[2] for j = 1:16], 1, :)
        fill!(st.fn, one(T))
        st.fields.B[3] .= T(0.2) .* sin.(y)
        st.fields.pe .= @. one(T) + T(0.1) * sin(x) * cos(y)
        for _ = 1:5
            HybridPlasmaPIC._advance_hybrid_fields!(st, T(0.1), 1)
        end
        # With frozen ui=0,n=1, Bz=0.2sin(y), the electron shear is
        # ue=(-0.2cos(y),0,0). Both the Hall and scalar-pressure curls vanish.
        expected = @. one(T) + T(0.1) * sin(x + T(0.1) * cos(y)) * cos(y)
        tol = T === Float32 ? T(3e-6) : T(2e-9)
        @test maximum(abs, st.fields.pe .- expected) < tol
        @test maximum(abs, st.fields.B[3] .- T(0.2) .* sin.(y)) < tol
        @test maximum(abs, st.fields.pe .- (@. one(T) + T(0.1) * sin(x) * cos(y))) > T(0.005)

        # Varying density makes div(J/n) nonzero despite div(J)=0.
        st.fn .= @. one(T) + T(0.1) * cos(x)
        fill!(st.fields.pe, T(0.7))
        HybridPlasmaPIC._ohm_ninv!(st.fields.ninv, st.fn, T(st.model.nfloor), st.fields.floor_count)
        HybridPlasmaPIC._bpressure_rhs!(
            st.k1,
            st.pressure_work.kp[1],
            st.fields.B,
            st.fields.pe,
            st,
            1,
            1,
        )
        divue = @. -T(0.02) * sin(x) * cos(y) / (one(T) + T(0.1) * cos(x))^2
        exact = @. -T(5 / 3) * T(0.7) * divue
        @test maximum(abs, st.pressure_work.kp[1] .- exact) < (T === Float32 ? T(2e-6) : T(3e-13))

        # Nonbarotropic pressure creates the Biermann battery from zero B.
        foreach(a -> fill!(a, zero(T)), st.fields.B)
        st.fields.pe .= @. one(T) + T(0.1) * sin(y)
        HybridPlasmaPIC._bpressure_rhs!(
            st.k1,
            st.pressure_work.kp[1],
            st.fields.B,
            st.fields.pe,
            st,
            1,
            1,
        )
        battery = @. T(0.01) * sin(x) * cos(y) / (one(T) + T(0.1) * cos(x))^2
        @test maximum(abs, st.k1[3] .- battery) < (T === Float32 ? T(2e-6) : T(3e-13))
        @test all(iszero, st.pressure_work.kp[1])
    end
end

@testset "Evolving pressure RK4 order under frozen constant moments" begin
    errors = Float64[]
    for nsteps in (2, 4, 8)
        g = FourierGrid((17,), (2pi,))
        st = HybridStepper(g, HybridModel(EvolvingElectrons()), CIC(), 0)
        x = [(i - 1) * g.dx[1] for i = 1:17]
        fill!(st.fn, 1.0)
        fill!(st.fui[1], 0.7)
        st.fields.pe .= @. 1 + 0.1 * cos(x)
        for _ = 1:nsteps
            HybridPlasmaPIC._advance_hybrid_fields!(st, 0.8 / nsteps, 1)
        end
        exact = @. 1 + 0.1 * cos(x - 0.7 * 0.8)
        push!(errors, maximum(abs, st.fields.pe .- exact))
    end
    @test all(r -> 14 < r < 18, errors[1:2] ./ errors[2:3])
    @info "evolving pressure frozen-moment RK4 errors" errors
end

@testset "Evolving pressure rejection and serial transaction" begin
    for T in (Float32, Float64), collection in (false, true)
        st, ps = ep_load(T; γ = T(3), compressive = true, collection)
        populations = collection ? ps : [ps]
        fields = map(copy, ep_fields(st.fields))
        identity = ep_fields(st.fields)
        particles = [(map(copy, p.x), map(copy, p.v)) for p in populations]
        clock = (st.time[], st.step[], st.last_dt[], st.fields.floor_count[])
        @test_throws ArgumentError step!(st, ps, T(2); NB = 1)
        @test all(isequal(a, b) for (a, b) in zip(ep_fields(st.fields), fields))
        @test all(a === b for (a, b) in zip(ep_fields(st.fields), identity))
        @test (st.time[], st.step[], st.last_dt[], st.fields.floor_count[]) == clock
        @test all(p.x == saved[1] && p.v == saved[2] for (p, saved) in zip(populations, particles))
        control, controlps = ep_load(T; γ = T(3), compressive = true, collection)
        step!(st, ps, T(0.001); NB = 2)
        step!(control, controlps, T(0.001); NB = 2)
        @test ep_fields(st.fields) == ep_fields(control.fields)
        for bad in (T(-1), T(NaN), T(Inf))
            st.fields.pe[1] = bad
            state = map(copy, ep_fields(st.fields))
            @test step!(st, ps, zero(T)) === st
            @test all(isequal(a, b) for (a, b) in zip(ep_fields(st.fields), state))
            @test_throws ArgumentError step!(st, ps, T(0.001))
            @test all(isequal(a, b) for (a, b) in zip(ep_fields(st.fields), state))
        end
        st.fields.pe = st.fields.n
        @test_throws ArgumentError step!(st, ps, T(0.001))
    end
end

@testset "Stored pressure public Ohm independence" begin
    for T in (Float32, Float64), target in (:n, :ui, :B, :E, :J, :gradp, :ninv, :lapJ)
        g = FourierGrid((7,), (T(2pi),))
        f = HybridFields{1,T}(g.n)
        fill!(f.n, one(T))
        f.pe .= T[1 + 0.1sin((i - 1) * g.dx[1]) for i = 1:7]
        old = getproperty(f, target)
        if old isa Tuple
            setproperty!(f, target, ntuple(c -> c == 1 ? f.pe : old[c], length(old)))
        else
            setproperty!(f, target, f.pe)
        end
        snapshot = map(copy, ep_fields(f))
        identity = ep_fields(f)
        floor = f.floor_count[]
        @test_throws ArgumentError ohms_law!(f, HybridModel(EvolvingElectrons(T(5 / 3))), g)
        @test all(isequal(a, b) for (a, b) in zip(ep_fields(f), snapshot))
        @test all(a === b for (a, b) in zip(ep_fields(f), identity))
        @test f.floor_count[] == floor
    end
    for T in (Float32, Float64), target in (:pforce, :cgl_dp)
        st, ps = ep_load(T)
        if target === :pforce
            st.fields.pforce = ntuple(_ -> zeros(T, st.g.n), 3)
        else
            st.fields.cgl_dp = zeros(T, st.g.n)
        end
        snapshot = map(copy, ep_fields(st.fields))
        @test_throws DimensionMismatch ohms_law!(st.fields, st.model, st.g)
        @test_throws DimensionMismatch step!(st, ps, T(0.001))
        @test all(isequal(a, b) for (a, b) in zip(ep_fields(st.fields), snapshot))
    end
    g = FourierGrid((7,), (Float32(2pi),))
    @test_throws ArgumentError HybridStepper(
        g,
        HybridModel(EvolvingElectrons(nextfloat(1.0))),
        CIC(),
        0,
    )
end

@testset "Pressure rollback after an accepted earlier subcycle" begin
    for T in (Float32, Float64), collection in (false, true)
        st, ps = ep_load(T; γ = T(3), compressive = true, collection)
        snapshot = map(copy, ep_fields(st.fields))
        identities = ep_fields(st.fields)
        populations = collection ? ps : (ps,)
        particles = [(map(copy, p.x), map(copy, p.v)) for p in populations]
        clock = (st.time[], st.step[], st.last_dt[], st.fields.floor_count[])
        err = try
            step!(st, ps, T(0.8); NB = 2)
            nothing
        catch error
            error
        end
        @test err isa ArgumentError
        @test occursin("subcycle 2", sprint(showerror, err))
        @test all(isequal(a, b) for (a, b) in zip(ep_fields(st.fields), snapshot))
        @test all(a === b for (a, b) in zip(ep_fields(st.fields), identities))
        @test (st.time[], st.step[], st.last_dt[], st.fields.floor_count[]) == clock
        @test all(
            isequal(p.x, s[1]) && isequal(p.v, s[2]) for (p, s) in zip(populations, particles)
        )
    end
end

@testset "Dynamic single population equals same-model collection" begin
    for T in (Float32, Float64)
        a, pa = ep_load(T)
        b, pb = ep_load(T; collection = true)
        for (dt, nb) in ((T(0.01), 1), (T(0.013), 2), (T(0.007), 4))
            step!(a, pa, dt; NB = nb)
            step!(b, pb, dt; NB = nb)
            @test ep_fields(a.fields) == ep_fields(b.fields)
            @test pa.x == only(pb).x && pa.v == only(pb).v
            @test (a.time[], a.step[], a.last_dt[]) == (b.time[], b.step[], b.last_dt[])
        end
        empty = ParticleSet{1,T}[]
        st = HybridStepper(a.g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
        init!(st, empty)
        step!(st, empty, T(0.01); NB = 2)
        @test all(iszero, st.fields.pe)
        @test st.fields.floor_count[] == prod(st.g.n)
    end
end
