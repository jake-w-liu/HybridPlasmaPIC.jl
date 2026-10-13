using HybridPlasmaPIC, Test, LinearAlgebra

@testset "HallMHDModel parameter validation" begin
    model = HallMHDModel(IsothermalElectrons(0.5); Ti = 0.2, η = 0.01, ηH = 0.02, nfloor = 1e-4)
    @test model.closure.Te == 0.5
    @test model.Ti == 0.2
    @test model.η == 0.01
    @test model.ηH == 0.02
    @test model.nfloor == 1e-4

    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); Ti = NaN)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); Ti = -0.1)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); η = NaN)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); η = -0.1)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); ηH = NaN)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); ηH = -0.1)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); nfloor = 0.0)
    @test_throws ArgumentError HallMHDModel(IsothermalElectrons(0.0); nfloor = NaN)
end

@testset "Hall-MHD Ohm's law matches hybrid field algebra" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    closure = IsothermalElectrons(0.4)
    hm = HallMHDState(g, HallMHDModel(closure; η = 0.03, ηH = 0.01, nfloor = 1e-5))
    hf = HybridFields{1,T}((n,))

    hm.fields.n .= @. 1.0 + 0.1 * cos(2 * x)
    hm.fields.ui[1] .= 0.2
    hm.fields.ui[2] .= @. 0.1 * sin(x)
    hm.fields.B[3] .= @. 1.0 + 0.05 * cos(3 * x)

    copyto!(hf.n, hm.fields.n)
    for c = 1:3
        copyto!(hf.ui[c], hm.fields.ui[c])
        copyto!(hf.B[c], hm.fields.B[c])
    end

    hall_mhd_ohms_law!(hm)
    ohms_law!(hf, HybridModel(closure; η = 0.03, ηH = 0.01, nfloor = 1e-5), g)

    for c = 1:3
        @test maximum(abs, hm.fields.E[c] .- hf.E[c]) < 1e-12
        @test maximum(abs, hm.fields.J[c] .- hf.J[c]) < 1e-12
    end
    @test hm.fields.floor_count[] == hf.floor_count[]
end

@testset "Hall-MHD supports CGL electron closure" begin
    T = Float64
    n = 32
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    cgl = HallMHDState(g, HallMHDModel(CGLElectrons(0.4, 0.4, 1.0, 1.0); η = 0.02, nfloor = 1e-5))
    iso = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.4); η = 0.02, nfloor = 1e-5))

    for st in (cgl, iso)
        fill!(st.fields.n, 1.0)
        st.fields.ui[1] .= @. 0.2 * sin(x)
        st.fields.ui[2] .= @. 0.1 * cos(x)
        st.fields.B[1] .= 1.0
        hall_mhd_ohms_law!(st)
    end

    @test length(cgl.fields.pforce[1]) == n
    for c = 1:3
        @test maximum(abs, cgl.fields.E[c] .- iso.fields.E[c]) < 1e-12
        @test maximum(abs, cgl.fields.J[c] .- iso.fields.J[c]) < 1e-12
    end
end

@testset "Hall-MHD continuity RHS: advected density wave" begin
    T = Float64
    n = 64
    L = 2π
    k = 2π / L
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); Ti = 0.0))
    A = 0.2
    U = 0.3
    st.fields.n .= @. 1 + A * cos(k * x)
    fill!(st.fields.ui[1], U)
    rhs = hall_mhd_rhs!(st)
    @test maximum(abs, rhs.dn .- (@. U * A * k * sin(k * x))) < 1e-12
end

@testset "Hall-MHD momentum RHS: magnetic pressure force" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); Ti = 0.0))
    fill!(st.fields.n, 1.0)
    k = 3.0
    st.fields.B[3] .= cos.(k .* x)
    rhs = hall_mhd_rhs!(st)
    @test maximum(abs, rhs.du[1] .- (@. k * sin(k * x) * cos(k * x))) < 1e-10
    @test maximum(abs, rhs.du[2]) < 1e-12
    @test maximum(abs, rhs.du[3]) < 1e-12
end

@testset "Hall-MHD ion-pressure force" begin
    T = Float64
    n = 64
    L = 2π
    k = 2π / L
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    Ti = 0.7
    A = 0.1
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); Ti = Ti))
    st.fields.n .= @. 1 + A * cos(k * x)
    rhs = hall_mhd_rhs!(st)
    expected = @. Ti * A * k * sin(k * x) / (1 + A * cos(k * x))
    @test maximum(abs, rhs.du[1] .- expected) < 1e-12
end

@testset "Hall-MHD ion-acoustic dispersion carries electron pressure" begin
    # Warm-electron oracle: with uniform B ∥ k the seeded density mode oscillates at
    # ω = k√(Ti+Te) (isothermal ions + electrons), so the closure's −∇p_e/n must
    # accelerate the bulk fluid — a Te=0 momentum test cannot see that coupling.
    n = 64
    L = 2π
    k = 2π / L
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    for (Te, Ti) in ((0.5, 0.0), (0.5, 0.3))
        st = HallMHDState(g, HallMHDModel(IsothermalElectrons(Te); Ti = Ti))
        ε = 1e-3
        st.fields.n .= @. 1 + ε * cos(k * x)
        st.fields.B[1] .= 1.0
        hall_mhd_ohms_law!(st)
        proj(s) = sum((s.fields.n[i] - 1) * cos(k * x[i]) for i = 1:n)
        ts = [0.0]
        as = [proj(st)]
        for _ = 1:1000                              # t = 20, ≳2 acoustic periods
            step_hall_mhd!(st, 0.02)
            push!(ts, st.time[])
            push!(as, proj(st))
        end
        # fit ω from the linearly interpolated zero crossings of a(t) ∝ cos(ωt)
        crossings = Float64[]
        for i = 1:length(as)-1
            if (as[i] > 0) != (as[i+1] > 0)
                push!(crossings, ts[i] + (ts[i+1] - ts[i]) * as[i] / (as[i] - as[i+1]))
            end
        end
        @test length(crossings) >= 4
        ωfit = π * (length(crossings) - 1) / (crossings[end] - crossings[1])
        @test isapprox(ωfit, k * sqrt(Ti + Te); rtol = 1e-3)   # measured rel.err ≈ 8e-8
    end
end

@testset "Hall-MHD CGL pressure-tensor force enters momentum" begin
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st = HallMHDState(g, HallMHDModel(CGLElectrons(0.4, 0.3, 1.0, 1.0); Ti = 0.0))
    st.fields.n .= @. 1 + 0.1 * cos(x)
    st.fields.B[1] .= 1.0                           # uniform B ⇒ J = 0 ⇒ du = -(∇·P_e)/n
    hall_mhd_ohms_law!(st)
    rhs = hall_mhd_rhs!(st)
    @test maximum(abs, rhs.du[1]) > 1e-3
    for c = 1:3
        @test maximum(abs, rhs.du[c] .- (@. -st.fields.pforce[c] / st.fields.n)) < 1e-13
    end
end

@testset "Hall-MHD RK4 uniform equilibrium and mutation safety" begin
    T = Float64
    g = FourierGrid((16, 8), (2π, 2π))
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.5); Ti = 0.2))
    fill!(st.fields.n, 1.0)
    fill!(st.fields.B[3], 1.0)
    hall_mhd_ohms_law!(st)

    n0 = copy(st.fields.n)
    u0 = ntuple(c -> copy(st.fields.ui[c]), 3)
    B0 = ntuple(c -> copy(st.fields.B[c]), 3)

    @test_throws ArgumentError step_hall_mhd!(st, NaN)
    @test_throws ArgumentError step_hall_mhd!(st, -0.1)
    @test st.fields.n == n0
    @test all(st.fields.ui[c] == u0[c] for c = 1:3)
    @test all(st.fields.B[c] == B0[c] for c = 1:3)

    @test step_hall_mhd!(st, 0.0) === st
    @test st.step[] == 0
    @test st.time[] == 0.0
    @test st.fields.n == n0
    @test all(st.fields.ui[c] == u0[c] for c = 1:3)
    @test all(st.fields.B[c] == B0[c] for c = 1:3)

    for _ = 1:5
        step_hall_mhd!(st, 0.05)
    end
    @test maximum(abs, st.fields.n .- n0) < 1e-12
    @test all(maximum(abs, st.fields.ui[c] .- u0[c]) < 1e-12 for c = 1:3)
    @test all(maximum(abs, st.fields.B[c] .- B0[c]) < 1e-12 for c = 1:3)
    @test st.step[] == 5
    @test st.time[] ≈ 0.25
end

@testset "Hall-MHD step conserves periodic total mass" begin
    T = Float64
    n = 64
    L = 2π
    k = 2π / L
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); Ti = 0.0))
    st.fields.n .= @. 1 + 0.05 * cos(k * x)
    st.fields.ui[1] .= @. 0.1 * sin(k * x)
    mass0 = sum(st.fields.n) * g.dx[1]
    for _ = 1:10
        step_hall_mhd!(st, 0.01)
    end
    mass1 = sum(st.fields.n) * g.dx[1]
    @test abs(mass1 - mass0) < 1e-11
    @test all(isfinite, st.fields.n)
end

@testset "Hall-MHD polytropic closure survives transient RK4 predictor undershoot" begin
    # An RK4 predictor stage can hand the RHS a density transiently below 0 even when the
    # committed state is healthy; the polytropic pressure `(n/n0)^γ` must not throw a
    # DomainError on a negative base (non-integer γ). Feed the RHS exactly such a stage
    # state: it must evaluate finite instead of throwing. (The RHS is tested directly
    # because grid-sampled advective undershoot needs u·k·dt beyond the RK4 stability
    # limit; the original stepped scenario survived its 5 steps only while the inert
    # pre-fix p_e left those unstable modes exactly unseeded, and died at step 9.)
    g = FourierGrid((16, 4), (2π, 2π))
    st = HallMHDState(g, HallMHDModel(PolytropicElectrons(0.5, 1.0, 5 / 3); Ti = 0.0, η = 0.0))
    xs = [(i - 1) * g.dx[1] for i = 1:16]
    for j = 1:4, i = 1:16
        st.fields.n[i, j] = 1 + 1.01 * cos(xs[i])      # min cell −0.01: an undershot stage
        st.fields.B[3][i, j] = 1.0
        st.fields.ui[1][i, j] = 0.5
    end
    rhs = hall_mhd_rhs!(st)                            # negative base reaches the closure
    @test all(isfinite, rhs.dn)
    @test all(all(isfinite, rhs.du[c]) for c = 1:3)
    @test all(all(isfinite, rhs.dB[c]) for c = 1:3)

    # And a committed step from a healthy warm-polytropic state completes and stays
    # finite (whistler ω·dt ≈ 1.3 and acoustic/advective ω·dt ≪ 1 at dt = 0.01).
    st2 = HallMHDState(g, HallMHDModel(PolytropicElectrons(0.5, 1.0, 5 / 3); Ti = 0.1))
    for j = 1:4, i = 1:16
        st2.fields.n[i, j] = 1 + 0.3 * cos(xs[i])
        st2.fields.B[3][i, j] = 1.0
        st2.fields.ui[1][i, j] = 0.5
    end
    for _ = 1:5
        step_hall_mhd!(st2, 0.01)
    end
    @test all(isfinite, st2.fields.n)
    @test minimum(st2.fields.n) > 0
end

function _hall_fields_snapshot(f)
    return map(fieldnames(typeof(f))) do name
        value = getfield(f, name)
        value isa Base.RefValue ? value[] : deepcopy(value)
    end
end

@testset "Hall-MHD final projection failure leaves all public fields unchanged" begin
    for T in (Float32, Float64), D = 1:3
        # Projection removes a purely longitudinal magnetic mode. CGL p_parallel
        # is finite before projection but exceeds T's range after |B| collapses;
        # the failed final Ohm solve must not commit the already-valid RK state.
        pressure = T === Float32 ? T(1e30) : T(1e285)
        dt = T === Float32 ? T(1e-35) : T(1e-300)
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
        st = HallMHDState(g, HallMHDModel(CGLElectrons(zero(T), pressure, one(T), one(T))))
        st.fields.n .= one(T)
        for I in CartesianIndices(st.fields.n)
            st.fields.B[1][I] = cos((I[1] - 1) * g.dx[1])
        end
        hall_mhd_ohms_law!(st)
        before = _hall_fields_snapshot(st.fields)
        references = map(name -> getfield(st.fields, name), fieldnames(typeof(st.fields)))
        @test_throws ArgumentError step_hall_mhd!(st, dt)
        @test _hall_fields_snapshot(st.fields) == before
        @test st.time[] == zero(T)
        @test st.step[] == 0
        # A retry without projection is representable. Public array references
        # remain valid across the successful transaction.
        step_hall_mhd!(st, dt; project_b = false)
        @test st.time[] == dt
        @test st.step[] == 1
        @test all(
            getfield(st.fields, name) === reference for
            (name, reference) in zip(fieldnames(typeof(st.fields)), references)
        )
        @test all(all(isfinite, v) for v in (st.fields.ui..., st.fields.B..., st.fields.E...))
    end
end

@testset "Hall-MHD rejected density candidate preserves derived fields" begin
    g = FourierGrid((16,), (2π,))
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0)))
    for i = 1:g.n[1]
        st.fields.n[i] = 1 + 0.2cos((i - 1) * g.dx[1])
    end
    st.fields.ui[1] .= 1.0
    hall_mhd_ohms_law!(st)
    before = _hall_fields_snapshot(st.fields)
    @test_throws ArgumentError step_hall_mhd!(st, 10.0)
    @test _hall_fields_snapshot(st.fields) == before
    @test (st.time[], st.step[]) == (0.0, 0)
end

@testset "Hall-MHD density floor is an inclusive committed-state bound" begin
    for T in (Float32, Float64),
        D = 1:3,
        closure in (
            IsothermalElectrons(0.2),
            PolytropicElectrons(0.2, 1.0, 5 / 3),
            CGLElectrons(0.2, 0.15, 1.0, 1.0),
        )
        # Uniform n,u,B is an exact stationary solution for all three closures.
        # Shared Ohm's law floors only n < nfloor, so equality must also be a
        # valid RK4 result; use adjacent representable densities to test the bound.
        nfloor = T(1e-4)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        st = HallMHDState(g, HallMHDModel(closure; Ti = 0.3, nfloor))
        st.fields.B[1] .= one(T)
        for density in (nfloor, nextfloat(nfloor))
            st.fields.n .= density
            hall_mhd_ohms_law!(st)
            @test st.fields.floor_count[] == 0
            @test all(==(inv(density)), st.fields.ninv)
            rhs = hall_mhd_rhs!(st)
            @test all(all(iszero, a) for a in (rhs.dn, rhs.du..., rhs.dB...))
            before = _hall_fields_snapshot(st.fields)
            step_hall_mhd!(st, T(0.1))
            @test _hall_fields_snapshot(st.fields) == before
        end
        @test st.time[] == T(0.1) + T(0.1)
        @test st.step[] == 2

        st.fields.n .= prevfloat(nfloor)
        hall_mhd_ohms_law!(st)
        @test st.fields.floor_count[] == length(st.fields.n)
        @test all(==(inv(nfloor)), st.fields.ninv)
        before = _hall_fields_snapshot(st.fields)
        @test_throws ArgumentError step_hall_mhd!(st, T(0.1))
        @test _hall_fields_snapshot(st.fields) == before
        @test st.time[] == T(0.1) + T(0.1)
        @test st.step[] == 2
    end
end

@testset "Hall-MHD rejects nonfinite final electric fields before commit" begin
    for T in (Float32, Float64),
        D = 1:3,
        closure in (IsothermalElectrons(0.0), PolytropicElectrons(0.0, 1.0, 5 / 3))

        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
        st = HallMHDState(g, HallMHDModel(closure))
        st.fields.n .= one(T)
        large = T(2) * sqrt(floatmax(T))
        # Disjoint u_y and longitudinal B_x give exactly zero RHS. Projection
        # spreads B_x into the occupied velocity plane, overflowing u_y*B_x
        # only in the final Ohm solve; the RK4 physical candidate stays finite.
        for I in CartesianIndices(st.fields.n)
            if I[1] == 1
                st.fields.ui[2][I] = large
            else
                st.fields.B[1][I] = large
            end
        end
        hall_mhd_ohms_law!(st)
        @test all(all(isfinite, a) for a in st.fields.E)
        rhs = hall_mhd_rhs!(st)
        @test all(all(iszero, a) for a in (rhs.dn, rhs.du..., rhs.dB...))
        before = _hall_fields_snapshot(st.fields)
        references = map(name -> getfield(st.fields, name), fieldnames(typeof(st.fields)))
        @test_throws ErrorException step_hall_mhd!(st, T(0.01))
        @test _hall_fields_snapshot(st.fields) == before
        @test (st.time[], st.step[]) == (zero(T), 0)
        @test all(all(isfinite, a) for a in (st.stage.n, st.stage.ui..., st.stage.B...))
        @test any(any(!isfinite, a) for a in st.stage.E)

        # The same physical state is admissible when projection is disabled.
        step_hall_mhd!(st, T(0.01); project_b = false)
        @test _hall_fields_snapshot(st.fields) == before
        @test (st.time[], st.step[]) == (T(0.01), 1)
        @test all(
            getfield(st.fields, name) === reference for
            (name, reference) in zip(fieldnames(typeof(st.fields)), references)
        )

        # Scaling both fields down keeps the projected electric field finite;
        # the validation must accept that neighboring physical configuration.
        st.fields.ui[2] ./= T(4)
        st.fields.B[1] ./= T(4)
        step_hall_mhd!(st, T(0.01))
        @test all(all(isfinite, a) for a in (st.fields.ui..., st.fields.B..., st.fields.E...))
        @test maximum(abs, st.fields.E[3]) > zero(T)
        @test (st.time[], st.step[]) == (T(0.01) + T(0.01), 2)
    end
end

# A circularly polarized wave has u·k = 0 and constant |B|, so the nonlinear
# equations reduce exactly to the two-mode Hall/Alfvén matrix used below.
function _hall_helical_wave_state(D, T; eta = 0.0, hyper = 0.0)
    g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> T(2π), D))
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); η = eta, ηH = hyper))
    b0 = [d <= D ? inv(sqrt(Float64(D))) : 0.0 for d = 1:3]
    ref = D == 1 ? [0.0, 1.0, 0.0] : [1.0, -1.0, 0.0] / sqrt(2.0)
    e2 = cross(b0, ref)
    st.fields.n .= 1
    for I in CartesianIndices(st.fields.n)
        phase = sum((I[d] - 1) * g.dx[d] for d = 1:D)
        for c = 1:3
            st.fields.B[c][I] = b0[c] + 0.01 * (ref[c] * cos(phase) + e2[c] * sin(phase))
        end
    end
    return st, b0, ref, e2
end
function _hall_helical_wave_amplitude(st, e1, e2)
    g = st.g
    D = length(g.n)
    v = zeros(ComplexF64, 2)
    for I in CartesianIndices(st.fields.n)
        phase = sum((I[d] - 1) * g.dx[d] for d = 1:D)
        for (s, f) in enumerate((st.fields.ui, st.fields.B))
            v[s] += sum((e1[c] + im * e2[c]) * f[c][I] for c = 1:3) * exp(-im * phase)
        end
    end
    v ./= length(st.fields.n)
    return v
end
@testset "Hall oblique helical wave exact two-mode matrix" begin
    for T in (Float32, Float64), D = 1:3, (eta, hyper) in ((0.0, 0.0), (0.03, 0.002))
        st, b0, e1, e2 = _hall_helical_wave_state(D, T; eta, hyper)
        initial = _hall_helical_wave_amplitude(st, e1, e2)
        for _ = 1:20
            step_hall_mhd!(st, T(0.01))
        end
        k = sqrt(D)
        M = ComplexF64[0 im*k; im*k im*k^2-eta*k^2-hyper*k^4]
        expected = exp(0.2M) * initial
        measured = _hall_helical_wave_amplitude(st, e1, e2)
        error = norm(measured - expected) / norm(expected)
        @test error < (T == Float32 ? 1e-4 : 2e-7)
        @test maximum(abs, st.fields.n .- 1) < (T == Float32 ? 1e-5 : 1e-12)
    end
end
@testset "Hall wave RK4 time convergence" begin
    for D = 1:3
        errors = Float64[]
        for N in (8, 16, 32)
            st, b0, e1, e2 = _hall_helical_wave_state(D, Float64)
            initial = _hall_helical_wave_amplitude(st, e1, e2)
            for _ = 1:N
                step_hall_mhd!(st, 0.4 / N)
            end
            k = sqrt(D)
            M = ComplexF64[0 im*k; im*k im*k^2]
            push!(errors, norm(_hall_helical_wave_amplitude(st, e1, e2) - exp(0.4M) * initial))
        end
        rates = log2.(errors[1:2] ./ errors[2:3])
        @test all(rates .> 3.8)
    end
end

@testset "Hall nonlinear velocity advection oracle" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 16, D), ntuple(_ -> T(2π), D))
        st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0)))
        st.fields.n .= 1
        amplitude = T.((0.12, -0.17, 0.09))
        k = ntuple(d -> T(d), D)
        projection = sum(k[d] * amplitude[d] for d = 1:D)
        expected_n = zeros(T, g.n)
        expected_u = ntuple(_ -> zeros(T, g.n), 3)
        for I in CartesianIndices(st.fields.n)
            phase = sum(k[d] * (I[d] - 1) * g.dx[d] for d = 1:D)
            expected_n[I] = -projection * cos(phase)
            for c = 1:3
                st.fields.ui[c][I] = amplitude[c] * sin(phase)
                expected_u[c][I] = -amplitude[c] * projection * sin(phase) * cos(phase)
            end
        end
        rhs = hall_mhd_rhs!(st)
        tolerance = T == Float32 ? 2e-6 : 2e-14
        @test maximum(abs, rhs.dn - expected_n) < tolerance
        for c = 1:3
            @test maximum(abs, rhs.du[c] - expected_u[c]) < tolerance
            @test all(iszero, rhs.dB[c])
        end
    end
end
