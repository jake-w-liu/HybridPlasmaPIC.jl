module SlabGyroOracleTests
using HybridPlasmaPIC, Test, FFTW, Serialization
const H = HybridPlasmaPIC

# Independent implementations of the documented formulas. The TSC kernel is the
# standard three-point piecewise-quadratic shape; node v sits at x = v*dx and
# maps to array index mod(v,n)+1. For B = B0 z-hat the gyro ring is the full
# phase set X + rho*(cos th_j, sin th_j, 0), th_j = 2*pi*j/nring.
function tsc1d(s::T) where {T}
    c = round(Int, s)
    δ = s - c
    return (c - 1, (T(1 // 2) * (T(1 // 2) - δ)^2, T(3 // 4) - δ^2, T(1 // 2) * (T(1 // 2) + δ)^2))
end
wrapcell(v, n) = mod(v, n) + 1
cellvol(g) = prod(g.L) / length(g.cbuf)
ring_point(X, ρ, j, nring) = (X[1] + ρ * cospi(2j / nring), X[2] + ρ * sinpi(2j / nring), X[3])

function deposit_test(st, species)
    T = eltype(st.B0)
    g = st.g
    out = zeros(T, g.n)
    stencils(p) = ntuple(d -> tsc1d(mod(p[d], g.L[d]) / g.dx[d]), 3)
    for p in species, i in eachindex(p.weight)
        X = ntuple(d -> p.X[d][i], 3)
        ρ = sqrt(2 * p.m * p.mu[i] / (p.q^2 * st.B0))
        for j = 0:(st.nring-1)
            sten = stencils(ring_point(X, ρ, j, st.nring))
            for O in CartesianIndices((3, 3, 3))
                index = ntuple(d -> wrapcell(sten[d][1] + O[d] - 1, g.n[d]), 3)
                out[index...] +=
                    p.q * p.weight[i] * prod(d -> sten[d][2][O[d]], 1:3) / (st.nring * cellvol(g))
            end
        end
    end
    return out
end

# Field equation chi*(-Delta_perp)phi + alpha*P(phi) = rho, P removing the
# y-z plane average at each x (all modes with ky == kz == 0), zero DC.
symbol_test(st, I) =
    st.chi * (st.g.kfull[1][I[1]]^2 + st.g.kfull[2][I[2]]^2) +
    (st.g.midx[2][I[2]] != 0 || st.g.midx[3][I[3]] != 0 ? st.alpha : zero(st.alpha))

function field_test(st, rho)
    rhohat = fft(rho)
    phihat = similar(rhohat)
    for I in CartesianIndices(phihat)
        phihat[I] = all(==(1), Tuple(I)) ? zero(eltype(phihat)) : rhohat[I] / symbol_test(st, I)
    end
    return ifft(phihat)
end

# TSC-interpolated value and interpolant derivative at a point (mu = 0 rings
# collapse to the centre, so the ring mean equals the centre sample).
function sample_test(field, point, g)
    T = eltype(field)
    sten = ntuple(d -> tsc1d(mod(point[d], g.L[d]) / g.dx[d]), 3)
    value = zero(T)
    grad = ntuple(_ -> zero(T), 3)
    for O in CartesianIndices((3, 3, 3))
        index = ntuple(d -> wrapcell(sten[d][1] + O[d] - 1, g.n[d]), 3)
        value += field[index...] * prod(d -> sten[d][2][O[d]], 1:3)
    end
    for d = 1:3
        others = filter(!=(d), (1, 2, 3))
        c = sten[d][1] + 1
        δ = T(mod(point[d], g.L[d]) / g.dx[d] - c)
        total = zero(T)
        for a = 1:3, b = 1:3
            index = ntuple(
                j ->
                    j == d ? wrapcell(c, g.n[j]) :
                    wrapcell(sten[j][1] + (j == others[1] ? a : b) - 1, g.n[j]),
                3,
            )
            left = Base.setindex(index, wrapcell(c - 1, g.n[d]), d)
            right = Base.setindex(index, wrapcell(c + 1, g.n[d]), d)
            total +=
                (
                    (δ - T(1 // 2)) * (field[left...] - field[index...]) +
                    (δ + T(1 // 2)) * (field[right...] - field[index...])
                ) *
                sten[others[1]][2][a] *
                sten[others[2]][2][b]
        end
        grad = Base.setindex(grad, total / g.dx[d], d)
    end
    return value, grad
end

function snapshot(st, species)
    io = IOBuffer()
    serialize(
        io,
        (
            rho = copy(st.rho_gc),
            phi = copy(st.phi),
            clock = (st.time, st.step, st.last_dt),
            species = [deepcopy((p.X, p.vpar, p.mu, p.weight, p.id, p.tag)) for p in species],
        ),
    )
    return take!(io)
end

# Uniform in x and z at every y plane, so an x-directed ExB drift translates
# each plane's lattice onto itself and rho stays fixed; phi then depends only
# on y and the strictly-positive electron-density gate is comfortably met.
function drift_fixture(::Type{T}) where {T}
    g = FourierGrid((8, 6, 4), (T(8), T(6), T(4)))
    p = GyrocentreSet{T}(192; q = one(T), m = T(2))
    i = 0
    for j = 1:6, iz = 1:4, ix = 1:8
        i += 1
        p.X[1][i] = T(ix - 1)
        p.X[2][i] = T(j - 1)
        p.X[3][i] = T(iz - 1)
        p.weight[i] = T(1 + T(1 // 4) * cospi(2 * (j - 1) / 6))
    end
    fill!(p.mu, zero(T))
    fill!(p.vpar, T(1 // 16))
    species = [p]
    st = SlabGyrokineticPIC(g, species; B0 = one(T), Te = T(1 // 2), nring = 8)
    return st, species
end

@testset "Ring-TSC deposit and Larmor radius" begin
    for T in (Float32, Float64)
        g = FourierGrid((7, 5, 4), (T(7), T(5), T(4)))
        p = GyrocentreSet{T}(3; q = T(2), m = T(3))
        p.X[1] .= T.((0.37, 2.9, 5.1))
        p.X[2] .= T.((1.8, 0.4, 3.3))
        p.X[3] .= T.((2.1, 1.2, 0.6))
        p.mu .= T.((0.03, 0.0, 0.11))
        p.weight .= T.((0.8, 1.2, 0.5))
        fill!(p.vpar, T(1 // 8))
        # Te must exceed e*max|P(phi)| for the strictly-positive electron
        # density gate; this sparse-charge fixture gives large potentials.
        st = SlabGyrokineticPIC(g, [p]; B0 = T(2), Te = T(512), nring = 12)
        for i = 1:3
            @test st._radius[1][i] ≈ sqrt(2 * p.m * p.mu[i] / (p.q^2 * st.B0)) rtol = T(4) * eps(T)
        end
        init!(st, [p])
        expected = deposit_test(st, [p])
        scale = maximum(abs, expected)
        @test all(
            isapprox.(st.rho_gc, expected; atol = T(96) * eps(T) * scale, rtol = T(96) * eps(T)),
        )
        @test abs(sum(st.rho_gc) * cellvol(g) - sum(p.q .* p.weight)) <=
              T(256) * eps(T) * sum(p.q .* p.weight)
    end
end

@testset "Full-k field solve, E×B/parallel drift, free-energy diagnostics" begin
    for T in (Float32, Float64)
        st, species = drift_fixture(T)
        p = species[1]
        init!(st, species)
        # Independent spectral solve agrees with the stored potential, the
        # stored symbol table and the fresh residual receipt.
        phi_expected = field_test(st, st.rho_gc)
        scale = max(maximum(abs, phi_expected), one(T))
        @test all(isapprox.(st.phi, phi_expected; atol = T(512) * eps(T) * scale))
        for I in CartesianIndices(st._lambda)
            all(==(1), Tuple(I)) && continue
            @test st._lambda[I] ≈ symbol_test(st, I) rtol = T(4) * eps(T)
        end
        rhohat = fft(st.rho_gc)
        phihat = fft(st.phi)
        N = length(st.phi)
        residual = maximum(
            abs,
            (symbol_test(st, I) * phihat[I] - rhohat[I]) / N for
            I in CartesianIndices(phihat) if Tuple(I) != (1, 1, 1)
        )
        source_norm = maximum(abs, rhohat[I] / N for I in CartesianIndices(phihat))
        # The production receipt uses its own scaled FFT path; the independent
        # FFTW residual should land in the same tolerance range, not bitwise.
        @test residual <= 4 * (st.field_atol + st.field_rtol * source_norm)
        @test st.last_report.residual <= st.last_report.tolerance
        # Expected per-plane ExB drift from the independent TSC gradient.
        rates = [sample_test(st.phi, ntuple(d -> p.X[d][i], 3), st.g)[2][2] for i = 1:192]
        h = T(1 // 16)
        x0 = copy(p.X[1])
        z0 = copy(p.X[3])
        step!(st, species, h)
        for i = 1:192
            δx = mod(p.X[1][i] - x0[i], st.g.L[1])
            δx > st.g.L[1] / 2 && (δx -= st.g.L[1])
            @test abs(δx - (-rates[i] / st.B0) * h) <= T(1024) * eps(T) * h
            @test p.X[3][i] ≈ mod(z0[i] + h * p.vpar[i], st.g.L[3]) atol = T(64) * eps(T)
        end
        # Translating uniform x-columns must not change the deposited profile.
        @test all(isapprox.(st.rho_gc, deposit_test(st, species); atol = T(256) * eps(T)))
        # Independent quadratic-form diagnostics.
        d = slab_gyrokinetic_diagnostics(st, species)
        phihat = fft(st.phi)
        k2perp(I) = st.g.kfull[1][I[1]]^2 + st.g.kfull[2][I[2]]^2
        nonzonal(I) = st.g.midx[2][I[2]] != 0 || st.g.midx[3][I[3]] != 0
        factor = cellvol(st.g) / (2N)
        Upol = factor * sum(st.chi * k2perp(I) * abs2(phihat[I]) for I in CartesianIndices(phihat))
        Ue =
            factor * sum(
                nonzonal(I) ? st.alpha * abs2(phihat[I]) : zero(st.alpha) for
                I in CartesianIndices(phihat)
            )
        pairing = cellvol(st.g) * sum(st.rho_gc .* st.phi) / 2
        Kpar = sum(p.weight .* p.m .* p.vpar .^ 2) / 2
        Kmu = sum(p.weight .* p.mu) * st.B0
        amp = max(abs(Upol), abs(Ue), abs(pairing), Kpar + Kmu, one(T))
        @test d.energy.polarization ≈ Upol atol = T(4096) * eps(T) * amp * N
        @test d.energy.electron_response_free_energy ≈ Ue atol = T(4096) * eps(T) * amp * N
        @test d.energy.field_pairing_defect ≈ Upol + Ue - pairing atol = T(8192) * eps(T) * amp * N
        @test d.energy.parallel ≈ Kpar atol = T(64) * eps(T) * amp
        @test d.energy.magnetic_moment ≈ Kmu atol = T(64) * eps(T) * amp
        @test d.energy.total ≈ Kpar + Kmu + Upol + Ue atol = T(8192) * eps(T) * amp * N
        @test d.charge[1] ≈ sum(p.q .* p.weight)
        @test d.parallel_momentum ≈ sum(p.weight .* p.m .* p.vpar)
    end
end

@testset "Invariant and domain guards" begin
    for T in (Float32, Float64)
        st, species = drift_fixture(T)
        p = species[1]
        init!(st, species)
        before = snapshot(st, species)
        p.weight[1] *= 2
        @test_throws ArgumentError step!(st, species, T(1 // 16))
        p.weight[1] /= 2
        p.mu[1] = T(1 // 10)
        @test_throws ArgumentError step!(st, species, T(1 // 16))
        p.mu[1] = zero(T)
        @test snapshot(st, species) == before
        # Constructor domain validation.
        bad = GyrocentreSet{T}(4; q = one(T), m = one(T))
        bad.mu[1] = -one(T)
        @test_throws ArgumentError SlabGyrokineticPIC(
            st.g,
            [bad];
            B0 = one(T),
            Te = one(T),
            nring = 4,
        )
        good = GyrocentreSet{T}(4; q = one(T), m = one(T))
        @test_throws ArgumentError SlabGyrokineticPIC(st.g, [good]; B0 = zero(T), Te = one(T))
        @test_throws ArgumentError SlabGyrokineticPIC(st.g, [good]; B0 = one(T), Te = zero(T))
        @test_throws ArgumentError SlabGyrokineticPIC(
            st.g,
            [good];
            B0 = one(T),
            Te = one(T),
            nring = 2,
        )
    end
end
end # module
