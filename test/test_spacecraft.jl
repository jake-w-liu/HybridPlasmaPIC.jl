using HybridPlasmaPIC, Test, LinearAlgebra

@testset "spacecraft / shock-frame diagnostics" begin

    @testset "gather_at: linear field a + b*x" begin
        n = 64
        L = 10.0
        g = FourierGrid((n,), (L,))
        dx = g.dx[1]
        a, b = 1.3, 0.7
        # nodes at (i-1)*dx
        xnodes = [(i - 1) * dx for i = 1:n]
        field = a .+ b .* xnodes
        # Interior points (away from the periodic wrap at the right edge) must
        # recover a + b*xpos to roundoff.
        for xpos in range(0.0, stop = (n - 2) * dx, length = 37)
            @test gather_at(field, g, xpos) ≈ a + b * xpos atol = 1e-12
        end
        # Exactly on a node returns that node's value.
        @test gather_at(field, g, xnodes[10]) ≈ field[10] atol = 1e-12
        # Periodic wrap: position past the last node blends node n with node 1.
        s = (n - 0.5) * dx                # halfway between node n and node 1 (wrapped)
        @test gather_at(field, g, s) ≈ 0.5 * (field[n] + field[1]) atol = 1e-12
        # Negative position wraps too.
        @test gather_at(field, g, -0.5 * dx) ≈ 0.5 * (field[n] + field[1]) atol = 1e-12
    end

    @testset "SyntheticProbe: sample! and advance!" begin
        n = 50
        L = 5.0
        g = FourierGrid((n,), (L,))
        dx = g.dx[1]
        xnodes = [(i - 1) * dx for i = 1:n]
        a, b = 2.0, -0.4
        field = a .+ b .* xnodes

        x0 = 1.234
        probe = SyntheticProbe(x0)
        @test probe.x == x0
        v = sample!(probe, field, g, 0.0)
        @test v ≈ a + b * x0 atol = 1e-12
        @test probe.t == [0.0]
        @test probe.val ≈ [a + b * x0] atol = 1e-12

        # Moving probe: advance then sample again.
        vx, dt = 0.5, 0.1
        newx = advance!(probe, vx, dt)
        @test newx ≈ x0 + vx * dt
        @test probe.x ≈ x0 + vx * dt
        sample!(probe, field, g, dt)
        @test length(probe.t) == 2
        @test probe.t[2] ≈ dt
        @test probe.val[2] ≈ a + b * (x0 + vx * dt) atol = 1e-12
    end

    @testset "shock_frame" begin
        @test shock_frame(3.0, 1.0) == 2.0
        @test shock_frame(-1.0, 2.0) == -3.0
        @test shock_frame(5.0f0, 2.0f0) === 3.0f0
    end

    @testset "dehoffmann_teller_velocity: residual flow ∥ B" begin
        # u with a genuine perpendicular component relative to B.
        for (u, B) in [
            ((1.0, 2.0, -0.5), (0.0, 0.0, 3.0)),
            ((2.0, -1.0, 0.7), (1.0, 1.0, 1.0)),
            ((-0.3, 0.4, 2.0), (2.0, -1.0, 0.5)),
        ]
            V = dehoffmann_teller_velocity(u, B)
            r = (u[1] - V[1], u[2] - V[2], u[3] - V[3])   # residual flow
            # (u - V_HT) must be parallel to B  ⇒  (u - V_HT) × B ≈ 0
            cx = r[2] * B[3] - r[3] * B[2]
            cy = r[3] * B[1] - r[1] * B[3]
            cz = r[1] * B[2] - r[2] * B[1]
            @test norm((cx, cy, cz)) < 1e-12
            # V_HT itself is perpendicular to B (it is the ⊥ part of u).
            @test abs(V[1] * B[1] + V[2] * B[2] + V[3] * B[3]) < 1e-12
        end

        # If u is already parallel to B, V_HT = 0 (no perpendicular flow).
        Vpar = dehoffmann_teller_velocity((2.0, 2.0, 2.0), (1.0, 1.0, 1.0))
        @test norm(Vpar) < 1e-12

        # Degenerate B = 0 ⇒ zero (undefined frame, handled gracefully).
        @test dehoffmann_teller_velocity((1.0, 2.0, 3.0), (0.0, 0.0, 0.0)) == (0.0, 0.0, 0.0)
    end

    @testset "classify_reflected: hand-built set" begin
        # Convention: upstream at +x, shock at x_shock moving +x at Vs.
        x_shock = 5.0
        Vs = 1.0
        ps = ParticleSet{1,Float64}(5)
        #            x      vx      expected reflected?
        # 1: upstream, vx-Vs > 0 (vx=3 > 1)            → reflected
        # 2: upstream, vx-Vs < 0 (vx=0.5 < 1)          → not (still moving downstream in frame)
        # 3: downstream (x<x_shock), vx-Vs > 0         → not (not upstream of shock)
        # 4: upstream, vx-Vs = 0 exactly (vx=1)        → not (strict >)
        # 5: upstream, large +vx                       → reflected
        ps.x[1] .= [6.0, 7.0, 2.0, 8.0, 9.0]
        ps.v[1] .= [3.0, 0.5, 3.0, 1.0, 4.0]
        flags = classify_reflected(ps, x_shock, Vs)
        @test flags == [true, false, false, false, true]
    end

end
