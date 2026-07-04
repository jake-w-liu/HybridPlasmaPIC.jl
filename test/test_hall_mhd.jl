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
    # An RK4 predictor stage can drive a cell's density transiently below 0 even when the
    # committed state is healthy; the polytropic pressure `(n/n0)^γ` must not throw a
    # DomainError on a negative base (non-integer γ). The step must complete and stay finite.
    g = FourierGrid((16, 4), (2π, 2π))
    st = HallMHDState(g, HallMHDModel(PolytropicElectrons(0.5, 1.0, 5 / 3); Ti = 0.0, η = 0.0))
    xs = [(i - 1) * g.dx[1] for i = 1:16]
    for j = 1:4, i = 1:16
        st.fields.n[i, j] = 1 + 0.9 * cos(xs[i])       # min committed density ~0.1 ≫ nfloor
        st.fields.B[3][i, j] = 1.0
        st.fields.ui[1][i, j] = 50.0                   # strong flow ⇒ predictor undershoot
    end
    for _ = 1:5
        step_hall_mhd!(st, 0.02)
    end
    @test all(isfinite, st.fields.n)
    @test minimum(st.fields.n) > 0
end
