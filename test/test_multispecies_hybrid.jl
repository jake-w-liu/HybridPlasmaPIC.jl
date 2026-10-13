using HybridPlasmaPIC, Test, LinearAlgebra

function multi_uniform(::Type{T}, ::Val{D}; shape = CIC(), integrator = HybridStepper) where {T,D}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
    densities = (T(0.8), T(0.3))
    velocities = ((zero(T), T(0.12), T(-0.04)), (zero(T), T(-0.08), T(0.09)))
    species = [ParticleSet{D,T}(prod(g.n) * s; q = T(s), m = T(s == 1 ? 1 : 4)) for s = 1:2]
    for (s, ps) in enumerate(species)
        p = 0
        for I in CartesianIndices(g.n), j = 1:s
            p += 1
            for d = 1:D
                fraction = d == 1 ? (T(j) - T(0.5)) / T(s) : T(0.5)
                ps.x[d][p] = (T(I[d] - 1) + fraction) * g.dx[d]
            end
            for c = 1:3
                ps.v[c][p] = velocities[s][c]
            end
        end
        ps.weight .= densities[s] * prod(g.L) / nparticles(ps)
    end
    st = integrator(g, HybridModel(IsothermalElectrons(0.1)), shape, 0)
    st.fields.B[1] .= T(0.8)
    return st, species, densities, velocities
end

multi_init!(st::HybridStepper, ps) = init!(st, ps)
multi_init!(st::CAMCLStepper, ps) = init_camcl!(st, ps)
multi_step!(st::HybridStepper, ps, dt; NB = 2) = step!(st, ps, dt; NB)
multi_step!(st::CAMCLStepper, ps, dt; NB = 2) = step_camcl!(st, ps, dt; NB)

function multi_public_state(st, species)
    fields = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    particles = map(species) do ps
        (ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    end
    cl = st isa CAMCLStepper ? (st.Blead, st.Blag, st.Bout, st.cl_h[], st.cl_live[]) : ()
    return deepcopy((fields, particles, st.time[], st.step[], st.last_dt[], cl))
end

@testset "Hybrid collection agrees with the scalar proton interface" begin
    for T in (Float32, Float64),
        D = 1:3,
        shape in (NGP(), CIC(), TSC()),
        St in (HybridStepper, CAMCLStepper)

        st, species, _, _ = multi_uniform(T, Val(D); shape, integrator = St)
        ps = species[1]
        for p in eachindex(ps.weight)
            ps.v[2][p] += T(0.002) * sin(ps.x[1][p])
        end
        reference = deepcopy(st)
        original = deepcopy(ps)
        multi_init!(st, [ps])
        multi_init!(reference, original)
        @test multi_public_state(st, [ps]) == multi_public_state(reference, [original])
        for dt in T.((0.001, 0.001, 0.0007, 0.0013))
            multi_step!(st, [ps], dt)
            multi_step!(reference, original, dt)
            @test multi_public_state(st, [ps]) == multi_public_state(reference, [original])
        end
    end
end

# Derived independently from the uniform Lorentz/Ohm equations. The total mass
# momentum is constant and the relative transverse velocity rotates at Ω.
function multi_uniform_oracle(species, densities, velocities, B, t)
    ρ = [Float64(species[s].m) * Float64(densities[s]) for s = 1:2]
    charge = [Float64(species[s].q) * Float64(densities[s]) for s = 1:2]
    a = charge ./ sum(charge)
    Ω = Float64(B) * (a[2] * species[1].q / species[1].m + a[1] * species[2].q / species[2].m)
    v0 = [Float64.(collect(v)) for v in velocities]
    center = (ρ[1] .* v0[1] + ρ[2] .* v0[2]) ./ sum(ρ)
    w = v0[1] - v0[2]
    θ = Ω * t
    rotated = [w[1], w[2] * cos(θ) + w[3] * sin(θ), w[3] * cos(θ) - w[2] * sin(θ)]
    integral = [
        w[1] * t,
        (w[2] * sin(θ) + w[3] * (1 - cos(θ))) / Ω,
        (w[3] * sin(θ) - w[2] * (1 - cos(θ))) / Ω,
    ]
    coefficients = (ρ[2] / sum(ρ), -ρ[1] / sum(ρ))
    return (
        [center + c .* rotated for c in coefficients],
        [center .* t + c .* integral for c in coefficients],
    )
end

function multi_uniform_error(::Type{T}, ::Val{D}, St, steps; variable = false) where {T,D}
    st, species, densities, velocities = multi_uniform(T, Val(D); integrator = St)
    initial_x = [deepcopy(ps.x) for ps in species]
    multi_init!(st, species)
    base = T(0.4 / steps)
    for i = 1:steps
        dt = variable ? base * T(isodd(i) ? 0.7 : 1.3) : base
        multi_step!(st, species, dt)
    end
    expected_v, _ =
        multi_uniform_oracle(species, densities, velocities, T(0.8), st.time[] - st.last_dt[] / 2)
    _, expected_x = multi_uniform_oracle(species, densities, velocities, T(0.8), st.time[])
    error = 0.0
    for (s, ps) in enumerate(species), c = 1:3
        error = max(error, maximum(abs.(ps.v[c] .- expected_v[s][c])))
    end
    for (s, ps) in enumerate(species), d = 1:D
        expected = mod.(initial_x[s][d] .+ expected_x[s][d], st.g.L[d])
        error = max(error, maximum(abs.(ps.x[d] .- expected)))
    end
    momentum_error = 0.0
    for c = 1:3
        initial_momentum = sum(species[s].m * densities[s] * velocities[s][c] for s = 1:2)
        final_momentum = sum(ps.m * sum(ps.weight .* ps.v[c]) / prod(st.g.L) for ps in species)
        momentum_error = max(momentum_error, abs(final_momentum - initial_momentum))
    end
    return (; error, momentum_error)
end

@testset "Two-ion uniform Lorentz/Ohm evolution and temporal order" begin
    for St in (HybridStepper, CAMCLStepper), D = 1:3, variable in (false, true)
        results = [multi_uniform_error(Float64, Val(D), St, n; variable) for n in (8, 16, 32)]
        errors = [r.error for r in results]
        momentum_errors = [r.momentum_error for r in results]
        @test errors[1] / errors[2] > 3.6
        @test errors[2] / errors[3] > 3.6
        @test errors[3] < 2e-5
        # The explicit carried-E predictor does not solve the coupled Boris
        # midpoint current implicitly. Momentum therefore has truncation error,
        # checked against the discrete force balance below and required to
        # converge at second order to the conserved continuum momentum.
        @test momentum_errors[1] / momentum_errors[2] > 3.6
        @test momentum_errors[2] / momentum_errors[3] > 3.6
        @test momentum_errors[3] < 1e-6
        @info "multispecies uniform convergence" integrator = St dimension = D variable errors =
            Tuple(errors) momentum_errors = Tuple(momentum_errors)
    end
    for St in (HybridStepper, CAMCLStepper), D = 1:3
        result = multi_uniform_error(Float32, Val(D), St, 16; variable = true)
        @test result.error < 1e-4
        @test result.momentum_error < 2e-5
    end
end

@testset "Hybrid population validation precedes all state mutation" begin
    for St in (HybridStepper, CAMCLStepper)
        st, species, _, _ = multi_uniform(Float64, Val(1); integrator = St)
        multi_init!(st, species)
        for corrupt! in (
            ps -> (ps.q = -1.0),
            ps -> (ps.m = 0.0),
            ps -> (ps.q = Inf),
            ps -> (ps.v[3][end] = NaN),
            ps -> (ps.x[1][end] = Inf),
            ps -> (ps.weight[end] = -1.0),
            ps -> resize!(ps.id, 1),
            ps -> resize!(ps.v[1], 1),
            ps -> (ps.v = (ps.v[1], ps.v[1], ps.v[3])),
        )
            bad = deepcopy(species)
            corrupt!(bad[2])
            before = multi_public_state(st, bad)
            @test_throws Union{ArgumentError,DimensionMismatch} multi_step!(st, bad, 0.001)
            @test isequal(multi_public_state(st, bad), before)
            @test_throws Union{ArgumentError,DimensionMismatch} multi_init!(st, bad)
            @test isequal(multi_public_state(st, bad), before)
            @test multi_step!(st, bad, 0.0) === st
            @test isequal(multi_public_state(st, bad), before)
        end
        alias = [species[1], species[1]]
        before = multi_public_state(st, alias)
        @test_throws ArgumentError multi_step!(st, alias, 0.001)
        @test multi_public_state(st, alias) == before
        for kind in (:repeated, :shared_position, :shared_component)
            a = ParticleSet{1,Float64}(0)
            b = ParticleSet{1,Float64}(0)
            if kind === :repeated
                b = a
            elseif kind === :shared_position
                b.x = a.x
            else
                a.v = (a.x[1], a.v[2], a.v[3])
            end
            bad = [a, b]
            before = multi_public_state(st, bad)
            @test_throws ArgumentError multi_init!(st, bad)
            @test multi_public_state(st, bad) == before
            @test_throws ArgumentError multi_step!(st, bad, 0.001)
            @test multi_public_state(st, bad) == before
        end
        st.step[] = typemax(Int)
        before = multi_public_state(st, species)
        @test_throws ArgumentError multi_step!(st, species, 0.001)
        @test multi_public_state(st, species) == before
    end
end

@testset "Two-ion momentum obeys the discrete Boris force balance" begin
    for St in (HybridStepper, CAMCLStepper), dt in (0.05, 0.025, 0.0125)
        st, species, densities, velocities = multi_uniform(Float64, Val(1); integrator = St)
        multi_init!(st, species)
        # C*v = v×B for the prescribed uniform B=(0.8,0,0). The Cayley
        # matrix below is derived from the centered Lorentz equation and
        # does not call the production Boris implementation.
        C = [0.0 0.0 0.0; 0.0 0.0 0.8; 0.0 -0.8 0.0]
        charge = [species[s].q * densities[s] for s = 1:2]
        mass = [species[s].m * densities[s] for s = 1:2]
        qm = [ps.q / ps.m for ps in species]
        v0 = [collect(v) for v in velocities]
        electric(v) = -C * sum(charge[s] * v[s] for s = 1:2) / sum(charge)
        E = electric(v0)
        minus = [v0[s] - dt / 2 * qm[s] * (E + C * v0[s]) for s = 1:2]
        for step = 1:2
            plus = [
                (I - dt / 2 * qm[s] * C) \ ((I + dt / 2 * qm[s] * C) * minus[s] + dt * qm[s] * E) for s = 1:2
            ]
            expected_change = dt * sum(charge[s] * (E + C * (plus[s] + minus[s]) / 2) for s = 1:2)
            before = sum(mass[s] * minus[s] for s = 1:2)
            multi_step!(st, species, dt)
            measured = [[sum(ps.v[c]) / nparticles(ps) for c = 1:3] for ps in species]
            change = sum(mass[s] * measured[s] for s = 1:2) - before
            @test isapprox(change, expected_change; atol = 2e-16, rtol = 1e-9)
            @test all(isapprox(measured[s], plus[s]; atol = 1e-15, rtol = 1e-14) for s = 1:2)
            base_E = electric(plus)
            predicted = [plus[s] + dt / 2 * qm[s] * (base_E + C * plus[s]) for s = 1:2]
            E = electric(predicted)
            @test isapprox([st.fields.E[c][1] for c = 1:3], E; atol = 1e-15, rtol = 1e-14)
            minus = plus
        end
    end
end

@testset "Hybrid collections retain empty and view-backed populations" begin
    for T in (Float32, Float64), D = 1:3, St in (HybridStepper, CAMCLStepper)
        st, species, _, _ = multi_uniform(T, Val(D); integrator = St)
        empty = ParticleSet{D,T}(0; q = T(2), m = T(4))
        pop = [species[1], empty]
        multi_init!(st, pop)
        multi_step!(st, pop, T(0.001))
        @test nparticles(empty) == 0
        @test st.step[] == 1
        none = typeof(species)()
        multi_init!(st, none)
        multi_step!(st, none, T(0.001))
        @test all(iszero, st.fields.n)
        @test all(all(iszero, c) for c in st.fields.ui)
        @test st.step[] == 1
        # Independent sliced arrays preserve the caller's population objects.
        viewed = map(species) do ps
            r = 1:nparticles(ps)
            ParticleSet{D,T}(
                map(a -> view(a, r), ps.x),
                map(a -> view(a, r), ps.v),
                view(ps.weight, r),
                view(ps.id, r),
                view(ps.tag, r),
                ps.q,
                ps.m,
            )
        end
        reference = deepcopy(st)
        original = deepcopy(species)
        multi_init!(st, viewed)
        multi_init!(reference, original)
        multi_step!(st, viewed, T(0.001))
        multi_step!(reference, original, T(0.001))
        @test multi_public_state(st, viewed) == multi_public_state(reference, original)
        @test length(st.work) == sum(nparticles, species)
    end
end

@testset "Hybrid collections reject empty wrappers sharing component or species storage" begin
    for T in (Float32, Float64), D = 1:3, St in (HybridStepper, CAMCLStepper)
        st, species, _, _ = multi_uniform(T, Val(D); integrator = St)
        multi_init!(st, species)
        multi_step!(st, species, T(0.001))
        for same_population in (false, true)
            empty = ParticleSet{D,T}(0)
            other = ParticleSet{D,T}(0)
            owner = same_population ? empty : other
            buffer = T === Float32 ? owner.tag : owner.id
            wrapped = ParticleSet{D,T}(
                empty.x,
                empty.v,
                reinterpret(T, buffer),
                empty.id,
                empty.tag,
                empty.q,
                empty.m,
            )
            population = [species[1], wrapped, other]
            before = multi_public_state(st, population)
            @test_throws ArgumentError multi_init!(st, population)
            @test multi_public_state(st, population) == before
            @test_throws ArgumentError multi_step!(st, population, T(0.001))
            @test multi_public_state(st, population) == before
        end
    end
end

function multi_wave_error(St, steps; variable = false)
    N, L, B0, amplitude = 16, 2π, 0.8, 1e-7
    g = FourierGrid((N,), (L,))
    density = (0.8, 0.3)
    species = [ParticleSet{1,Float64}(N; q = s, m = s == 1 ? 1 : 4) for s = 1:2]
    x = collect(0:(N-1)) .* g.dx[1]
    initial = amplitude .* ComplexF64[0.4+0.1im, -0.3im, -0.2+0.2im, 0.3, 0.5im, 0.4]
    for (s, ps) in enumerate(species)
        ps.x[1] .= x
        ps.weight .= density[s] * L / N
        for c = 2:3
            ps.v[c] .= real.(initial[2(s-1)+c-1] .* exp.(im .* x))
        end
    end
    st = St(g, HybridModel(IsothermalElectrons(0.0)), NGP(), 0)
    st.fields.B[1] .= B0
    st.fields.B[2] .= real.(initial[5] .* exp.(im .* x))
    st.fields.B[3] .= real.(initial[6] .* exp.(im .* x))

    # Independent Fourier-mode system: J=ik×B, E=−U×B0+J×B0/ne,
    # dv_s/dt=(q_s/m_s)(E+v_s×B0), dB/dt=−ik×E. NGP keeps
    # these tiny-amplitude particles on their original nodes, avoiding a
    # spatial interpolation error floor in the temporal convergence oracle.
    charges = [species[s].q * density[s] for s = 1:2]
    ne = sum(charges)
    function linear_rhs(z)
        uy = sum(charges[s] * z[2s-1] for s = 1:2) / ne
        uz = sum(charges[s] * z[2s] for s = 1:2) / ne
        ey = -B0 * uz + im * B0 * z[5] / ne
        ez = B0 * uy + im * B0 * z[6] / ne
        out = zeros(ComplexF64, 6)
        for s = 1:2
            qm = species[s].q / species[s].m
            out[2s-1] = qm * (ey + B0 * z[2s])
            out[2s] = qm * (ez - B0 * z[2s-1])
        end
        out[5], out[6] = im * ez, -im * ey
        return out
    end
    basis = Matrix{ComplexF64}(I, 6, 6)
    generator = hcat((linear_rhs(basis[:, c]) for c = 1:6)...)
    multi_init!(st, species)
    for i = 1:steps
        dt = 0.4 / steps * (variable ? (isodd(i) ? 0.7 : 1.3) : 1.0)
        multi_step!(st, species, dt; NB = variable && isodd(i) ? 3 : 2)
    end
    expected_v = exp(generator * (st.time[] - st.last_dt[] / 2)) * initial
    expected_b = exp(generator * st.time[]) * initial
    mode(a) = 2sum(a .* exp.(-im .* x)) / N
    error = maximum(abs(mode(species[s].v[c]) - expected_v[2(s-1)+c-1]) for s = 1:2, c = 2:3)
    error = max(
        error,
        abs(mode(st.fields.B[2]) - expected_b[5]),
        abs(mode(st.fields.B[3]) - expected_b[6]),
    )
    return error / amplitude
end

@testset "Multispecies transverse waves match an independent linear system" begin
    for St in (HybridStepper, CAMCLStepper), variable in (false, true)
        errors = [multi_wave_error(St, n; variable) for n in (8, 16, 32)]
        @test errors[1] / errors[2] > 3.4
        @test errors[2] / errors[3] > 3.4
        @test errors[3] < 5e-4
        @info "multispecies wave convergence" integrator = St variable errors
    end
end

@testset "Hybrid charge moments floor the combined density" begin
    for T in (Float32, Float64), D = 1:3, St in (HybridStepper, CAMCLStepper)
        st, species, density, velocity = multi_uniform(T, Val(D); integrator = St)
        st.model = HybridModel(st.model.closure; nfloor = 1.0)
        charge = [species[s].q * density[s] for s = 1:2]
        expected_n = sum(charge)
        expected_u = ntuple(c -> sum(charge[s] * velocity[s][c] for s = 1:2) / expected_n, 3)
        multi_init!(st, species)
        @test all(abs.(st.fields.n .- expected_n) .< 20eps(T))
        for c = 1:3
            @test all(abs.(st.fields.ui[c] .- expected_u[c]) .< 20eps(T))
        end
        @test st.fields.floor_count[] == 0
        @test all(abs.(st.fields.E[2] .+ T(0.8) * expected_u[3]) .< 20eps(T))
        @test all(abs.(st.fields.E[3] .- T(0.8) * expected_u[2]) .< 20eps(T))
    end
end

@testset "Hybrid population order and CGL compatibility" begin
    for T in (Float32, Float64), D = 1:3, St in (HybridStepper, CAMCLStepper)
        st, species, _, _ = multi_uniform(T, Val(D); integrator = St)
        reference, reversed = deepcopy(st), reverse(deepcopy(species))
        multi_init!(st, species)
        multi_init!(reference, reversed)
        for dt in T.((0.001, 0.0007, 0.0013))
            multi_step!(st, species, dt)
            multi_step!(reference, reversed, dt)
        end
        for (a, b) in zip(
            (st.fields.n, st.fields.ui..., st.fields.B..., st.fields.E...),
            (
                reference.fields.n,
                reference.fields.ui...,
                reference.fields.B...,
                reference.fields.E...,
            ),
        )
            @test isapprox(a, b; atol = 64eps(T), rtol = 64eps(T))
        end
        for (a, b) in zip(species, reverse(reversed)),
            (x, y) in zip((a.x..., a.v...), (b.x..., b.v...))

            @test isapprox(x, y; atol = 64eps(T), rtol = 64eps(T))
        end
    end
    for T in (Float32, Float64), D = 1:3
        initial, species, _, _ = multi_uniform(T, Val(D))
        model = HybridModel(CGLElectrons(0.1, 0.12, 1.0, 1.0); de2 = 0.02)
        st = HybridStepper(initial.g, model, CIC(), 0)
        st.fields.B[1] .= T(0.8)
        reference, original = deepcopy(st), deepcopy(species[1])
        init!(st, [species[1]])
        init!(reference, original)
        for dt in T.((0.001, 0.0007, 0.0013))
            step!(st, [species[1]], dt; NB = 2)
            step!(reference, original, dt; NB = 2)
            @test multi_public_state(st, [species[1]]) == multi_public_state(reference, [original])
        end
        @test_throws ArgumentError CAMCLStepper(initial.g, model, CIC(), 0)
    end
end

@testset "Hybrid collection workspaces resize after particle splitting" begin
    for St in (HybridStepper, CAMCLStepper)
        st, species, _, _ = multi_uniform(Float64, Val(1); integrator = St)
        multi_init!(st, species)
        multi_step!(st, species, 0.001)
        reference, original = deepcopy(st), deepcopy(species)
        ps = species[2]
        split = ParticleSet{1,Float64}(2nparticles(ps); q = ps.q, m = ps.m)
        for p in eachindex(ps.weight), child = 1:2
            i = 2(p - 1) + child
            split.x[1][i] = ps.x[1][p]
            for c = 1:3
                split.v[c][i] = ps.v[c][p]
            end
            split.weight[i] = ps.weight[p] / 2
        end
        species[2] = split
        multi_step!(st, species, 0.0013)
        multi_step!(reference, original, 0.0013)
        @test length(st.work) == sum(nparticles, species)
        for (a, b) in zip(
            (st.fields.n, st.fields.ui..., st.fields.B..., st.fields.E...),
            (
                reference.fields.n,
                reference.fields.ui...,
                reference.fields.B...,
                reference.fields.E...,
            ),
        )
            @test isapprox(a, b; atol = 1e-13, rtol = 1e-13)
        end
        for c = 1:3, child = 1:2
            @test isapprox(split.v[c][child:2:end], original[2].v[c]; atol = 1e-13, rtol = 1e-13)
        end
    end
end
