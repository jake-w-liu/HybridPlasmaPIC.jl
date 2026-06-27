# test_shock_diag.jl — shock diagnostics on small CONSTRUCTED inputs (no sims).
# Verifies: surface spectrum peaks at the seeded k_y; transverse coherence of a
# sinusoid matches cos; CrossingLogger counts a hand-built crossing set exactly
# and energy_gain matches the analytic ΔKE; NIF boost removes the tangential
# flow; boundary_reflection_fraction matches a hand-counted set; boundary energy
# flux matches a hand-set field/moment state.

using HybridPlasmaPIC, Test

@testset "shock diagnostics" begin
    T = Float64

    # ---------------- shock_surface_spectrum: peaks at the seeded k_y ----------
    @testset "surface spectrum peaks at seeded k_y" begin
        nx, ny = 16, 32
        Lx, Ly = 10.0, 8.0
        sh = PerpShock2D(nx, ny, Lx, Ly; B0 = 1.0)
        # Build a Bz field whose per-column shock front x_s(y) is a pure m=3
        # sinusoid. shock_surface walks from the downstream peak outward and
        # records the first node where Bz drops below (B0+peak)/2. We set, per
        # column j, a compressed downstream block [1, ip(j)] at Bz=2 and upstream
        # at Bz=B0=1, so the front sits at node ip(j)+1 → x = (ip(j))*dx.
        mseed = 3
        dx = sh.sbp.dx
        # choose base block so the sinusoid stays interior
        amp_nodes = 2                      # ± nodes of ripple
        base = nx ÷ 2
        ipk = zeros(Int, ny)
        for j = 1:ny
            shift = round(Int, amp_nodes * cos(2π * mseed * (j - 1) / ny))
            ip = base + shift
            ipk[j] = ip
            for i = 1:nx
                sh.Bz[i, j] = i <= ip ? 2.0 : 1.0
            end
        end
        spec = shock_surface_spectrum(sh)
        # peak (excluding DC) should be at m = mseed
        Pnon = copy(spec.Ps)
        Pnon[1] = -Inf                      # ignore DC
        kpk = argmax(Pnon) - 1              # 0-based mode index of the peak
        @test kpk == mseed
        # ky value at the peak matches 2π·m/Ly
        @test isapprox(spec.ky[mseed+1], 2π * mseed / Ly; rtol = 1e-12)
        # DC bin ≈ 0 because the mean is removed
        @test spec.Ps[1] < 1e-18
        # front actually located at the expected nodes
        xs_expected = T[ipk[j] * dx for j = 1:ny]
        @test isapprox(spec.xs, xs_expected; atol = 1e-12)
    end

    # ---------------- transverse_coherence: sinusoid → cos --------------------
    @testset "transverse coherence of a sinusoid matches cos" begin
        nx, ny = 16, 64
        Lx, Ly = 10.0, 16.0
        sh = PerpShock2D(nx, ny, Lx, Ly; B0 = 1.0)
        mseed = 4
        base = nx ÷ 2
        amp = 3.0
        # Set x_s(y) directly by constructing Bz so the front node tracks a
        # sinusoid; but the coherence test only needs x_s, so build Bz with a
        # continuous threshold crossing. Easiest: front node = round(base + amp·cos).
        for j = 1:ny
            ip = base + round(Int, amp * cos(2π * mseed * (j - 1) / ny))
            for i = 1:nx
                sh.Bz[i, j] = i <= ip ? 2.0 : 1.0
            end
        end
        coh = transverse_coherence(sh)
        @test coh.C[1] ≈ 1.0
        # For x_s ∝ cos(2π m y/Ly), C_s(Δy) = cos(2π m Δy/Ly) (exact for the
        # rounded front too, since rounding is symmetric — compare to the
        # analytic autocorrelation of the ACTUAL discrete x_s).
        xs, m, _ = shock_surface(sh)
        f = xs .- m
        var0 = sum(abs2, f) / ny
        for lag = 0:ny-1
            acc = 0.0
            for jx = 1:ny
                jj = mod(jx - 1 + lag, ny) + 1
                acc += f[jx] * f[jj]
            end
            @test coh.C[lag+1] ≈ (acc / ny) / var0
        end
        # and the dominant shape is cosine at the seeded mode: correlation with
        # cos(2π m Δy/Ly) is near 1.
        cosref = T[cos(2π * mseed * coh.dy[k] / Ly) for k = 1:ny]
        corr = sum(coh.C .* cosref) / sqrt(sum(abs2, coh.C) * sum(abs2, cosref))
        @test corr > 0.99
    end

    # ---------------- CrossingLogger: exact count + analytic ΔKE --------------
    @testset "CrossingLogger counts crossings and ΔKE exactly" begin
        # 3 particles, flat surface at x_surface = 5.0.
        # p1: 4 → 6 (crosses up), then 6 → 4 (crosses down)  → 2 crossings
        # p2: 4 → 4.5 (no cross)                              → 0 crossings
        # p3: 6 → 4 (crosses down)                            → 1 crossing
        xsurf = 5.0
        ps = ParticleSet{1,T}(3)
        ps.id .= UInt64[10, 20, 30]
        ps.m = 1.0

        logger = CrossingLogger(T)

        # call 1 — registration (no crossings)
        ps.x[1] .= [4.0, 4.0, 6.0]
        ps.v[1] .= [1.0, 0.0, 0.0]
        ps.v[2] .= [0.0, 0.0, 0.0]
        ps.v[3] .= [0.0, 0.0, 0.0]
        @test log_crossings!(logger, ps, xsurf) == 0
        @test crossing_count(logger) == 0

        # call 2 — p1 4→6 (cross up), p3 6→4 (cross down); set velocities so the
        # KE changes are analytic.
        ps.x[1] .= [6.0, 4.5, 4.0]
        # p1 KE before = ½(1²)=0.5; now set v=(2,0,0) → KE=2.0, ΔKE_p1=+1.5
        # p3 KE before = ½(0)=0;    now set v=(0,3,0) → KE=4.5, ΔKE_p3=+4.5
        ps.v[1] .= [2.0, 0.0, 0.0]
        ps.v[2] .= [0.0, 0.0, 3.0]
        ps.v[3] .= [0.0, 0.0, 0.0]
        @test log_crossings!(logger, ps, xsurf) == 2
        @test crossing_count(logger) == 2
        @test energy_gain(logger) ≈ 1.5 + 4.5

        # call 3 — p1 6→4 (cross down). Set p1 v=(0,0,1)→KE=0.5, ΔKE since call2 =
        # 0.5 − 2.0 = −1.5. p2, p3 do not cross.
        ps.x[1] .= [4.0, 4.6, 3.5]
        ps.v[1] .= [0.0, 0.0, 0.0]
        ps.v[2] .= [0.0, 0.0, 0.0]
        ps.v[3] .= [1.0, 0.0, 0.0]
        @test log_crossings!(logger, ps, xsurf) == 1
        @test crossing_count(logger) == 3
        @test energy_gain(logger) ≈ (1.5 + 4.5) + (0.5 - 2.0)
    end

    @testset "CrossingLogger with per-particle (rippled) surface" begin
        # per-particle surface; only p1 crosses its own surface.
        ps = ParticleSet{1,T}(2)
        ps.id .= UInt64[1, 2]
        ps.m = 1.0
        logger = CrossingLogger(T)
        surf = [3.0, 7.0]
        ps.x[1] .= [2.0, 8.0]                  # p1 below 3, p2 above 7
        log_crossings!(logger, ps, surf)
        ps.x[1] .= [4.0, 8.0]                  # p1 now above 3 (cross), p2 stays above
        ps.v[2] .= [2.0, 0.0]                  # p1 KE 0→2 ⇒ ΔKE=+2
        @test log_crossings!(logger, ps, surf) == 1
        @test crossing_count(logger) == 1
        @test energy_gain(logger) ≈ 2.0
    end

    # ---------------- normal_incidence_frame ----------------------------------
    @testset "NIF boost removes the tangential flow" begin
        n_hat = (1.0, 0.0, 0.0)               # shock normal = x
        u = (-3.0, 1.7, -0.4)                 # inflow with tangential (y,z) parts
        B = (0.0, 0.0, 1.0)                   # perpendicular field (unused geometrically)
        Vnif = normal_incidence_frame(u, B, n_hat)
        # boost is purely tangential
        @test Vnif[1] ≈ 0.0 atol = 1e-14
        @test Vnif[2] ≈ 1.7
        @test Vnif[3] ≈ -0.4
        # flow seen in NIF is purely normal
        urel = (u[1] - Vnif[1], u[2] - Vnif[2], u[3] - Vnif[3])
        @test urel[1] ≈ -3.0
        @test hypot(urel[2], urel[3]) < 1e-13   # tangential component ≈ 0

        # oblique, non-unit normal n̂ = (1,1,0)/√2 given as (1,1,0)
        nh = (1.0, 1.0, 0.0)
        u2 = (2.0, 0.0, 5.0)
        V2 = normal_incidence_frame(u2, B, nh)
        ur2 = (u2[1] - V2[1], u2[2] - V2[2], u2[3] - V2[3])
        # residual flow must be parallel to nh (tangential part = 0): cross product ≈ 0
        cx = ur2[2] * nh[3] - ur2[3] * nh[2]
        cy = ur2[3] * nh[1] - ur2[1] * nh[3]
        cz = ur2[1] * nh[2] - ur2[2] * nh[1]
        @test hypot(cx, cy, cz) < 1e-13

        # degenerate normal
        @test normal_incidence_frame(u, B, (0.0, 0.0, 0.0)) == (0.0, 0.0, 0.0)
    end

    # ---------------- boundary_reflection_fraction ----------------------------
    @testset "boundary_reflection_fraction matches the count (1D)" begin
        N = 24
        Lx = 12.0
        sh = PerpShock(N, Lx; B0 = 1.0)
        dx = sh.s.dx                          # = Lx/(N-1)
        ncells = 3
        edge = Lx - ncells * dx
        ps = ParticleSet{1,T}(6)
        # band is [edge, Lx]. Place: 3 in band (2 moving +x = upstream/back), 3 out.
        ps.x[1] .= [Lx, Lx - dx, Lx - 2dx, edge - 0.5, 1.0, 5.0]
        ps.v[1] .= [1.0, 2.0, -1.0, 5.0, 1.0, 1.0]   # in-band: +,+,−
        frac = boundary_reflection_fraction(sh, ps; ncells = ncells)
        @test frac ≈ 2 / 3
        # no particles in band → 0
        ps2 = ParticleSet{1,T}(2)
        ps2.x[1] .= [0.0, 1.0]
        ps2.v[1] .= [1.0, 1.0]
        @test boundary_reflection_fraction(sh, ps2; ncells = ncells) == 0.0
        @test_throws ArgumentError boundary_reflection_fraction(sh, ps; ncells = 0)
        @test_throws ArgumentError boundary_reflection_fraction(sh, ps; ncells = -1)
    end

    @testset "boundary_reflection_fraction matches the count (2D)" begin
        nx, ny = 20, 8
        Lx, Ly = 10.0, 4.0
        sh = PerpShock2D(nx, ny, Lx, Ly; B0 = 1.0)
        dx = sh.sbp.dx
        ncells = 2
        edge = Lx - ncells * dx
        ps = ParticleSet{2,T}(4)
        ps.x[1] .= [Lx, Lx - dx, edge - 0.1, 2.0]   # first two in band
        ps.x[2] .= [1.0, 1.0, 1.0, 1.0]
        ps.v[1] .= [3.0, -1.0, 9.0, 9.0]            # in-band: +, −  → 1/2
        frac = boundary_reflection_fraction(sh, ps; ncells = ncells)
        @test frac ≈ 0.5
        @test_throws ArgumentError boundary_reflection_fraction(sh, ps; ncells = 0)
    end

    # ---------------- boundary_energy_flux (1D) -------------------------------
    @testset "boundary_energy_flux on a hand-set state" begin
        N = 8
        Lx = 7.0
        sh = PerpShock(N, Lx; B0 = 1.0)
        # hand-set fields/moments at the two boundary nodes (1 = wall, N = inflow)
        sh.Bz[1] = 2.0
        sh.Bz[N] = 1.0
        sh.Ey[1] = 0.5
        sh.Ey[N] = -3.0
        sh.n[1] = 4.0
        sh.n[N] = 1.0
        sh.ux[1] = 0.0
        sh.ux[N] = -3.0
        sh.uy[1] = 0.1
        sh.uy[N] = 0.2
        sh.pe[1] = 3.0
        sh.pe[N] = 2.0
        bf = boundary_energy_flux(sh)
        # entries are (inflow=node N, wall=node 1)
        # magnetic (Poynting) S_x = E_y B_z
        @test bf.magnetic[1] ≈ -3.0 * 1.0      # inflow
        @test bf.magnetic[2] ≈ 0.5 * 2.0       # wall
        # kinetic F_K = ½ n (ux²+uy²) ux
        kin_wall = 0.5 * 4.0 * (0.0^2 + 0.1^2) * 0.0
        kin_inflow = 0.5 * 1.0 * ((-3.0)^2 + 0.2^2) * (-3.0)
        @test bf.kinetic[1] ≈ kin_inflow
        @test bf.kinetic[2] ≈ kin_wall
        # inflow kinetic flux is negative (upstream plasma streams in −x)
        @test bf.kinetic[1] < 0
        # electron enthalpy flux γe/(γe−1)·pe·ux  (γe=5/3 ⇒ factor 2.5)
        c = sh.γe / (sh.γe - 1)
        @test bf.enthalpy[1] ≈ c * 2.0 * (-3.0)     # inflow
        @test bf.enthalpy[2] ≈ c * 3.0 * 0.0        # wall
        @test bf.total[1] ≈ bf.magnetic[1] + bf.kinetic[1] + bf.enthalpy[1]
        @test bf.total[2] ≈ bf.magnetic[2] + bf.kinetic[2] + bf.enthalpy[2]
    end

    @testset "boundary_energy_flux full kinetic-ion flux (ps) conserves the moment" begin
        # the deposited ion energy flux integrates (H-weighted) to the exact
        # particle sum Σ_p w·½|v|²·v_x — and includes the thermal part the bulk
        # estimate ½n|u|²u_x misses.
        N = 16
        Lx = 10.0
        sh = PerpShock(N, Lx; B0 = 1.0)
        Np = 200
        ps = ParticleSet{1,Float64}(Np)
        for p = 1:Np
            ps.x[1][p] = Lx * (p - 0.5) / Np           # spread across the box
            ps.v[1][p] = -2.0 + 0.3 * sin(0.7p)        # drift + deterministic spread
            ps.v[2][p] = 0.4 * cos(0.9p)
            ps.v[3][p] = 0.2 * sin(1.3p)
        end
        ps.weight .= shock_density_weight(1.0, Lx, Np)
        Fe = HybridPlasmaPIC._ion_energy_flux(sh, ps)
        # H-weighted grid integral == exact particle moment (CIC partition of unity)
        moment = sum(
            ps.weight[p] * 0.5 * (ps.v[1][p]^2 + ps.v[2][p]^2 + ps.v[3][p]^2) * ps.v[1][p] for
            p = 1:Np
        )
        @test sum(Fe[i] * sh.s.H[i] for i = 1:N) ≈ moment rtol = 1e-12
        @test all(isfinite, Fe)
    end
end
