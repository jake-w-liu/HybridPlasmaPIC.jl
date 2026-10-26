using Test, Random
using HybridPlasmaPIC
isdefined(@__MODULE__, :_composition_fixture) || include("collisional_fixtures.jl")

@testset "Collisional production invariants in a uniform cell" begin
    # One equal-weight, equal-mass pair has no spatial gradients. Ohm's law
    # makes E=-U×B, so the bulk force vanishes. The first Euler backup followed
    # by Boris increases the perpendicular relative radius by sqrt(1+t^2),
    # t=Omega*h/2; subsequent Boris rotations and collisions preserve it.
    for T in (Float32, Float64), R in (Xoshiro, MersenneTwister), omega in (0.0, 0.2)
        h = T(0.03)
        g = FourierGrid((1,), (one(T),))
        p = ParticleSet{1,T}(2; q = one(T), m = one(T))
        p.x[1] .= (T(0.25), T(0.75))
        bulk = (T(0.2), T(-0.1), T(0.05))
        for c = 1:3
            p.v[c][1] = bulk[c] + (c == 1 ? T(0.5) : zero(T))
            p.v[c][2] = bulk[c] - (c == 1 ? T(0.5) : zero(T))
        end
        inner = HybridStepper(g, HybridModel(IsothermalElectrons(zero(T))), NGP(), 2)
        inner.fields.B[3][1] = T(omega)
        st = CollisionalHybridStepper(
            inner,
            [p];
            dt = h,
            rng = R(619),
            coulomb_log = one(T),
            epsilon0 = one(T),
            weight_to_number = one(T),
        )
        init!(st, [p])
        expected_radius = sqrt(1 + (Float64(T(omega)) * Float64(h) / 2)^2)
        expected_E = (-Float64(bulk[2]) * T(omega), Float64(bulk[1]) * T(omega), 0.0)
        expected_energy = sum(abs2, Float64.(bulk)) + expected_radius^2 / 4
        for k = 1:16
            step!(st, [p], h)
            tolerance = 256eps(T) * k
            relative = ntuple(c -> Float64(p.v[c][1]) - Float64(p.v[c][2]), 3)
            measured_bulk = ntuple(c -> (Float64(p.v[c][1]) + Float64(p.v[c][2])) / 2, 3)
            energy = sum(sum(Float64(v)^2 for v in component) for component in p.v) / 2
            @test abs(hypot(relative...) - expected_radius) <= tolerance
            @test all(abs(measured_bulk[c] - bulk[c]) <= tolerance for c = 1:3)
            @test all(abs(st.inner.fields.E[c][1] - expected_E[c]) <= tolerance for c = 1:3)
            @test abs(energy - expected_energy) <= tolerance
            @test st.inner.fields.B[3][1] === T(omega)
            @test st.source_count == st.inner.step[] == k
            @test all(x -> 0 <= x < 1, p.x[1])
        end
    end
end

@testset "Collisional carried E drives the next kick" begin
    # After one step the composition and a stale-E variant (canonical step and
    # collision, no recomputation) share x, v, B and RNG and differ only in the
    # carried E. The next Boris kick is affine in E: dv = (h/2)(q/m)(I+R)dE with a
    # rotation R, |R-I| <= h*Omega. So the next-step velocity difference is
    # h(q/m)dE within (sqrt(3)/2)h*Omega*|h(q/m)dE|_inf plus rounding. Its size is
    # set by the collisional change of the charge-weighted bulk velocity, not by
    # the precision, and need not exceed a fixed multiple of eps(Float32).
    H = HybridPlasmaPIC
    for D = 1:3, T in (Float32, Float64), strong in (false, true)
        composed, cps = _composition_feedback_fixture(Val(D), T, Xoshiro; strong)
        actual, aps = _composition_feedback_fixture(Val(D), T, Xoshiro; strong)
        stale, sps = _composition_feedback_fixture(Val(D), T, Xoshiro; strong)
        step!(composed, cps, composed.dt)
        _composition_component!(actual, aps)
        step!(stale.inner, sps, stale.dt)
        collide_coulomb!(
            sps,
            stale.grid,
            stale.dt;
            coulomb_log = stale.coulomb_log,
            epsilon0 = stale.epsilon0,
            weight_to_number = stale.weight_to_number,
            u_floor = stale.u_floor,
            rng = stale.rng,
            workspace = stale.workspace,
        )
        @test isequal(_composition_physical(composed, cps), _composition_physical(actual, aps))
        @test all(isequal(p.x, q.x) && isequal(p.v, q.v) for (p, q) in zip(aps, sps))
        E = actual.inner.fields.E
        dE = ntuple(c -> Float64.(E[c]) .- Float64.(stale.inner.fields.E[c]), 3)
        Escale = max(1.0, maximum(maximum(abs, a) for a in E))
        @test maximum(maximum(abs, a) for a in dE) > 128eps(T) * Escale
        B = map(copy, actual.inner.fields.B)
        x = [map(copy, p.x) for p in aps]
        step!(actual.inner, aps, actual.dt)
        step!(stale.inner, sps, stale.dt)
        h = Float64(actual.dt)
        g = actual.grid
        misfit = signal = 0.0
        for (a, (p, q)) in enumerate(zip(aps, sps))
            qm = Float64(p.q) / Float64(p.m)
            for i in eachindex(p.weight)
                node = ntuple(d -> mod(floor(Int, x[a][d][i] / g.dx[d] + 0.5), g.n[d]) + 1, D)
                I = CartesianIndex(node)
                omega = qm * sqrt(sum(Float64(B[c][I])^2 for c = 1:3))
                predicted = ntuple(c -> h * qm * dE[c][I], 3)
                observed = ntuple(c -> Float64(p.v[c][i]) - Float64(q.v[c][i]), 3)
                speed = maximum(c -> abs(Float64(p.v[c][i])), 1:3)
                bound = h * omega * maximum(abs, predicted) + 8eps(T) * (1 + speed)
                misfit = max(misfit, maximum(abs, observed .- predicted) / bound)
                signal = max(signal, maximum(abs, predicted))
            end
        end
        @test misfit <= 1
        @test signal > 0
    end
end
