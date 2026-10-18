using HybridPlasmaPIC, LinearAlgebra, Test
# Direct, odd-grid discrete Fourier differentiation, independent of package FFTs.
function ep_dft_derivative(n)
    x = 2pi .* (0:n-1) ./ n
    modes = -(n ÷ 2):(n÷2)
    return [real(sum(im * k * cis(k * (xi - xj)) for k in modes)) / n for xi in x, xj in x]
end
function ep_reference_rhs(z, n, ui, γ, Dx, Dy)
    dx(a) = Dx * a
    dy(a) = a * transpose(Dy)
    bx, by, bz, p = eachslice(z; dims = 3)
    j = (dy(bz), -dx(bz), dx(by) - dy(bx))
    u = ntuple(c -> ui[c] .- j[c] ./ n, 3)
    e = (
        -(u[2] .* bz - u[3] .* by) .- dx(p) ./ n,
        -(u[3] .* bx - u[1] .* bz) .- dy(p) ./ n,
        -(u[1] .* by - u[2] .* bx),
    )
    return cat(
        -dy(e[3]),
        dx(e[3]),
        -dx(e[2]) + dy(e[1]),
        -u[1] .* dx(p) .- u[2] .* dy(p) .- γ .* p .* (dx(u[1]) + dy(u[2]));
        dims = 3,
    )
end
function ep_midpoint_reference(z0, n, ui, γ, Dx, Dy, duration, steps)
    z = copy(z0)
    h = duration / steps
    for _ = 1:steps
        k = ep_reference_rhs(z, n, ui, γ, Dx, Dy)
        z .+= h .* ep_reference_rhs(z .+ (h / 2) .* k, n, ui, γ, Dx, Dy)
    end
    return z
end
function ep_joint_probe()
    g = FourierGrid((7, 5), (2pi, 2pi))
    st = HybridStepper(g, HybridModel(EvolvingElectrons()), CIC(), 0)
    x = reshape((0:6) .* g.dx[1], :, 1)
    y = reshape((0:4) .* g.dx[2], 1, :)
    n = @. 1 + 0.1 * cos(x) * cos(y)
    ui = (repeat(0.12 .* sin.(y), 7, 1), repeat(0.08 .* cos.(x), 1, 5), @. 0.04 * sin(x) * cos(y))
    z0 = cat(
        repeat(0.2 .+ 0.03 .* cos.(y), 7, 1),
        repeat(0.1 .+ 0.02 .* sin.(x), 1, 5),
        (@. 0.4 + 0.1 * sin(x) + 0.05 * cos(y)),
        (@. 0.7 + 0.04 * cos(x) + 0.03 * sin(y));
        dims = 3,
    )
    dx, dy = ep_dft_derivative(7), ep_dft_derivative(5)
    # Periodic skew-adjoint differentiation transfers electron pressure work:
    # sum(dp/dt)/(gamma-1) = -sum(p*div(ue)). This is a semidiscrete
    # source identity, not a claim of exact particle timestep conservation.
    reference_rate = ep_reference_rhs(z0, n, ui, 5 / 3, dx, dy)
    j = (
        z0[:, :, 3] * transpose(dy),
        -dx * z0[:, :, 3],
        dx * z0[:, :, 2] - z0[:, :, 1] * transpose(dy),
    )
    ue = ntuple(c -> ui[c] .- j[c] ./ n, 3)
    divergence = dx * ue[1] + ue[2] * transpose(dy)
    work = -sum(z0[:, :, 4] .* divergence)
    @test abs(sum(reference_rate[:, :, 4]) / (5 / 3 - 1) - work) < 128eps(Float64)
    st.fn .= n
    for c = 1:3
        st.fui[c] .= ui[c]
        st.fields.B[c] .= z0[:, :, c]
    end
    st.fields.pe .= z0[:, :, 4]
    HybridPlasmaPIC._ohm_ninv!(st.fields.ninv, st.fn, st.model.nfloor, st.fields.floor_count)
    HybridPlasmaPIC._bpressure_rhs!(
        st.k1,
        st.pressure_work.kp[1],
        st.fields.B,
        st.fields.pe,
        st,
        1,
        1,
    )
    production_rate = cat(st.k1..., st.pressure_work.kp[1]; dims = 3)
    @test maximum(abs, production_rate .- reference_rate) < 128eps(Float64)
    @test abs(sum(st.pressure_work.kp[1]) / (5 / 3 - 1) - work) < 128eps(Float64)
    @test abs(work) > 1e-3
    refs = [ep_midpoint_reference(z0, n, ui, 5 / 3, dx, dy, 0.4, k) for k in (512, 1024, 2048)]
    extrap1 = (4refs[2] - refs[1]) / 3
    extrap2 = (4refs[3] - refs[2]) / 3
    @test maximum(abs, extrap1 .- extrap2) < 2e-11
    @test maximum(abs, dx * sin.(repeat(x, 1, 5)) .- cos.(repeat(x, 1, 5))) < 2e-14
    errors = []
    for steps in (2, 4, 8)
        st.fn .= n
        for c = 1:3
            st.fui[c] .= ui[c]
            st.fields.B[c] .= z0[:, :, c]
        end
        st.fields.pe .= z0[:, :, 4]
        for _ = 1:steps
            HybridPlasmaPIC._advance_hybrid_fields!(st, 0.4 / steps, 1)
        end
        got = cat(st.fields.B..., st.fields.pe; dims = 3)
        push!(errors, maximum(abs, got .- extrap2))
    end
    @test all(r -> 12 < r < 20, errors[1:2] ./ errors[2:3])
    @test errors[end] < 2e-9
    @info "joint B/pressure RK4 vs direct-DFT extrapolated midpoint reference" errors reference_delta =
        maximum(abs, extrap1 .- extrap2)
end
function ep_acoustic_error(steps; multispecies = false, variable = false, nb = 1)
    N = 17
    g = FourierGrid((N,), (2pi,))
    x = (0:N-1) .* g.dx[1]
    qs = multispecies ? (1.0, 2.0) : (1.0,)
    ms = multispecies ? (1.0, 4.0) : (1.0,)
    ns = multispecies ? (0.8, 0.3) : (1.0,)
    ne = sum(qs .* ns)
    γ = 5 / 3
    p0 = 0.7
    amplitude = 1e-7
    u0 = multispecies ? ComplexF64[0.31+0.07im, -0.23+0.11im] : ComplexF64[0.31+0.07im]
    z0 = vcat(u0, ComplexF64[0.6-0.13im])
    count = length(qs)
    mat = zeros(ComplexF64, count + 1, count + 1)
    populations = ParticleSet{1,Float64}[]
    for s = 1:count
        ps = ParticleSet{1,Float64}(N; q = qs[s], m = ms[s])
        ps.x[1] .= x
        set_density_weight!(ps, ns[s], g)
        ps.v[1] .= amplitude .* real.(u0[s] .* cis.(x))
        push!(populations, ps)
        mat[s, end] = -im * qs[s] / (ms[s] * ne)
        mat[end, s] = -im * γ * p0 * qs[s] * ns[s] / ne
    end
    st = HybridStepper(g, HybridModel(EvolvingElectrons(γ)), NGP(), N * count)
    st.fields.pe .= p0 .+ amplitude .* real.(z0[end] .* cis.(x))
    ps = multispecies ? populations : only(populations)
    init!(st, ps)
    duration = 0.8
    for k = 1:steps
        dt = duration / steps * (variable ? (isodd(k) ? 0.7 : 1.3) : 1.0)
        step!(st, ps, dt; NB = nb)
    end
    ionref = exp(mat * (st.time[] - st.last_dt[] / 2)) * z0
    pref = exp(mat * st.time[]) * z0
    mode(a) = 2sum(a .* cis.(-x)) / N / amplitude
    error = maximum(
        vcat(
            [abs(mode(populations[s].v[1]) - ionref[s]) for s = 1:count],
            abs(mode(st.fields.pe .- p0) - pref[end]),
        ),
    )
    displacement =
        maximum(maximum(min.(abs.(p.x[1] .- x), g.L[1] .- abs.(p.x[1] .- x))) for p in populations)
    @test 0 < displacement < g.dx[1] / 2
    @test maximum(abs, st.fields.n .- ne) < 8eps(Float64) * ne
    return error
end
@testset "Evolving pressure joint stages vs independent discrete DFT" begin
    ep_joint_probe()
end

@testset "Evolving pressure moving-ion acoustic convergence" begin
    # Linearized cold-ion pressure equations give
    # omega^2 = gamma*p0*sum(q_s^2*n_s/m_s)/ne^2 for k=1.
    # One NGP particle per cell moves without crossing a cell boundary, keeping
    # the spatial map fixed. Amplitude 1e-7 separates temporal error from the
    # nonlinear O(amplitude^2) correction. Expected ions use their half time.
    for multi in (false, true), variable in (false, true), nb in (1, 2, 4)
        errors = [ep_acoustic_error(k; multispecies = multi, variable, nb) for k in (8, 16, 32)]
        @test all(r -> 3.4 < r < 4.6, errors[1:2] ./ errors[2:3])
        @test errors[end] < 5e-5
    end
end

function pressure_energy_error(steps; variable = false)
    N = 33
    g = FourierGrid((N,), (2pi,))
    x = (0:N-1) .* g.dx[1]
    ps = ParticleSet{1,Float64}(N)
    ps.x[1] .= x
    set_density_weight!(ps, 1.0, g)
    ps.v[1] .= 0.03 .* cos.(x)
    st = HybridStepper(g, HybridModel(EvolvingElectrons()), NGP(), N)
    st.fields.pe .= 0.7 .+ 0.04 .* sin.(x)
    init!(st, ps)
    initial =
        energy_budget(ps, st.fields.B, st.fields.n, st.model.closure, g; pe = st.fields.pe).total
    duration = 0.4
    for k = 1:steps
        dt = duration / steps * (variable ? (isodd(k) ? 0.7 : 1.3) : 1.0)
        step!(st, ps, dt; NB = 2)
    end
    p = copy(st.fields.pe)
    vprev = copy(ps.v[1])
    prevdt = st.last_dt[]
    nextdt = duration / steps * (variable ? 0.7 : 1.0)
    step!(st, ps, nextdt; NB = 2)
    vc = (nextdt .* vprev .+ prevdt .* ps.v[1]) ./ (prevdt + nextdt)
    energy = sum(ps.weight .* vc .^ 2) / 2 + sum(p) / (5 / 3 - 1) * prod(g.dx)
    # Particle displacement stays within the original nearest-grid-point cell;
    # deposition/gather therefore keep the exact same semidiscrete spatial map.
    shift = maximum(min.(abs.(ps.x[1] .- x), g.L[1] .- abs.(ps.x[1] .- x)))
    return abs(energy - initial), shift
end
@testset "Evolving pressure common-time ideal energy convergence" begin
    for variable in (false, true)
        results = [pressure_energy_error(k; variable) for k in (8, 16, 32)]
        errors = first.(results)
        @test all(r -> 3.4 < r < 4.6, errors[1:2] ./ errors[2:3])
        @test all(last(r) < pi / 33 for r in results)
        @test errors[end] < 1e-7
        @info "common-time ideal energy error" variable errors
    end
end


@testset "Evolving pressure genuine 3D mixed-grid translation" begin
    for T in (Float32, Float64)
        g = FourierGrid((7, 9, 11), (T(2pi), T(2pi), T(2pi)))
        st = HybridStepper(g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
        x = reshape(T.(0:6) .* g.dx[1], 7, 1, 1)
        y = reshape(T.(0:8) .* g.dx[2], 1, 9, 1)
        z = reshape(T.(0:10) .* g.dx[3], 1, 1, 11)
        phase = @. x + 2y + 3z
        velocity = (T(0.13), T(-0.07), T(0.11))
        omega = velocity[1] + 2velocity[2] + 3velocity[3]
        fill!(st.fn, one(T))
        for c = 1:3
            fill!(st.fui[c], velocity[c])
        end
        st.fields.pe .= @. one(T) + T(0.07) * cos(phase)
        initial = copy(st.fields.pe)
        h = T(0.025)
        for _ = 1:12
            HybridPlasmaPIC._advance_hybrid_fields!(st, h, 1)
        end
        exact = @. one(T) + T(0.07) * cos(phase - omega * (12h))
        # The analytic Fourier phase uses all three spatial derivatives. The
        # Float64 tolerance exceeds the RK4 phase truncation (~3e-13 here);
        # Float32 includes FFT/accumulation roundoff over twelve steps.
        tolerance = T === Float32 ? T(64) * eps(T) : T(4e-12)
        @test maximum(abs, st.fields.pe .- exact) < tolerance
        @test maximum(abs, st.fields.pe .- initial) > T(0.005)
        @test maximum(maximum(abs, b) for b in st.fields.B) < T(64) * eps(T)
        @test all(p -> isfinite(p) && p > zero(T), st.fields.pe)
    end
end
