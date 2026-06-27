# PUSH-001 extras: oblique initial velocity (v∥ preserved, v⊥ gyrates) and a
# long-time (1000-gyroperiod) stability test.

using HybridPlasmaPIC, Test

@testset "Boris: oblique initial velocity" begin
    T = Float64
    B0 = 1.0
    dt = 0.05
    ps = ParticleSet{2,T}(1)
    vperp = 1.3
    vpar = 0.7
    ps.x[1][1] = 0
    ps.x[2][1] = 0
    ps.v[1][1] = vperp
    ps.v[2][1] = 0
    ps.v[3][1] = vpar   # v⊥ in x, v∥ along z=B
    E = (0.0, 0.0, 0.0)
    B = (0.0, 0.0, B0)
    for _ = 1:5000
        push_uniform!(ps, E, B, dt)
    end
    @test isapprox(ps.v[3][1], vpar; atol = 1e-12)                 # v∥ unchanged (no force ∥B)
    vperp_now = sqrt(ps.v[1][1]^2 + ps.v[2][1]^2)
    @test isapprox(vperp_now, vperp; atol = 1e-10)                 # |v⊥| conserved
end

@testset "Boris: 1000-gyroperiod stability" begin
    T = Float64
    B0 = 1.0
    Ωc = B0
    dt = 0.02
    ps = ParticleSet{2,T}(1)
    ps.v[1][1] = 1.0
    ps.v[2][1] = 0.0
    ps.v[3][1] = 0.0
    E = (0.0, 0.0, 0.0)
    B = (0.0, 0.0, B0)
    period = 2π / Ωc
    nsteps = round(Int, 1000 * period / dt)             # 1000 gyroperiods
    spdmax = 0.0
    spdmin = Inf
    for _ = 1:nsteps
        push_uniform!(ps, E, B, dt)
        s = sqrt(ps.v[1][1]^2 + ps.v[2][1]^2)
        spdmax = max(spdmax, s)
        spdmin = min(spdmin, s)
    end
    @test (spdmax - spdmin) < 1e-10                     # speed conserved to roundoff over 1000 orbits
    @test isfinite(ps.x[1][1]) && isfinite(ps.x[2][1])  # bounded, no blow-up
end
