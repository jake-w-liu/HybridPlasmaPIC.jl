using HybridPlasmaPIC, Test

# Exact Bernstein certificates on the four origin/edge triangles. Nonpositive
# Bernstein coefficients imply |R(-a+ib)|² <= 1 throughout each triangle, not
# just at sampled frequencies. The positive-coordinate polygon is downward
# closed, so independent upper bounds on damping and frequency are sufficient.
const _RK4_STABILITY_VERTICES = (
    (0 // 1, 14 // 5),
    (1 // 1, 5 // 2),
    (37 // 20, 37 // 20),
    (5 // 2, 1 // 1),
    (557 // 200, 0 // 1),
)

@testset "RK4 mixed-rate polygon has an exact stability certificate" begin
    Q = Rational{BigInt}
    C = Complex{Q}
    function multiply(p, q)
        r = Dict{Tuple{Int,Int},C}()
        for ((i, j), x) in p, ((k, l), y) in q
            key = (i + k, j + l)
            r[key] = get(r, key, zero(C)) + x * y
        end
        return r
    end
    degree = 16
    # Certify a slightly larger polygon as well, leaving room for rounding of
    # coefficients and the logarithmic timestep evaluation.
    expansion = 1 + Q(1, 2^20)
    for edge = 1:4
        a, b = _RK4_STABILITY_VERTICES[edge:edge+1]
        nx, ny = a[2] - b[2], b[1] - a[1]
        intercept = nx * a[1] + ny * a[2]
        @test nx > 0 && ny > 0
        @test all(v -> nx * v[1] + ny * v[2] <= intercept, _RK4_STABILITY_VERTICES)
        z = Dict((1, 0) => expansion * C(-a[1], a[2]), (0, 1) => expansion * C(-b[1], b[2]))
        power = Dict((0, 0) => one(C))
        polynomial = copy(power)
        for k = 1:4
            power = multiply(power, z)
            for (ij, coefficient) in power
                polynomial[ij] = get(polynomial, ij, zero(C)) + coefficient / factorial(k)
            end
        end
        norm_squared = multiply(polynomial, Dict(ij => conj(x) for (ij, x) in polynomial))
        norm_squared[(0, 0)] -= 1
        @test all(isreal, values(norm_squared))
        for i = 0:degree, j = 0:degree-i
            coefficient = sum(
                real(x) * Q(
                    binomial(i, k) * binomial(j, l),
                    binomial(degree, k + l) * binomial(k + l, k),
                ) for ((k, l), x) in norm_squared if k <= i && l <= j
            )
            @test coefficient <= 0
        end
    end
end

@testset "recommended_dt controls combined Hall waves and damping" begin
    amplification(z) = abs(1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24)
    for T in (Float32, Float64)
        g = FourierGrid((33,), (T(2π),))
        K = π / Float64(g.dx[1])
        wave_bound = K * (sqrt(K^2 + 4) + K) / 2
        η = (2.785 / 2.8) * wave_bound / K^2
        dt = recommended_dt(g; η, safety = 0.8)
        k = 16 * 2π / Float64(g.L[1])
        @test amplification(complex(-η * k^2 * dt, k^2 * dt)) <= 1

        # Check the complete rectangle of independently bounded modal rates,
        # including hyperresistive and mixed damping, all supported dimensions,
        # anisotropic grids, safety=1, and particle/magnetic subcycling.
        for D = 1:3, safety in (0.2, 0.8, 1.0), NB in (1, 3)
            grid = FourierGrid(ntuple(_ -> 4, D), ntuple(d -> T(2π * d), D))
            corner = sqrt(sum(dx -> (π / Float64(dx))^2, grid.dx))
            ω = corner * (sqrt(corner^2 + 4) + corner) / 2
            for ratio in (1e-12, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e12)
                for partition in (0.0, 0.4, 1.0)
                    η = ratio * ω * partition / corner^2
                    ηH = ratio * ω * (1 - partition) / corner^4
                    dt = recommended_dt(grid; η, ηH, safety, NB)
                    ν = η * corner^2 + ηH * corner^4
                    for x in (0.0, 0.2, 0.7, 1.0), y in (0.0, 0.3, 0.8, 1.0)
                        @test amplification(complex(-x * ν, y * ω) * (dt / NB)) <=
                              1 + 32eps(Float64)
                    end
                end
            end
        end
    end
end

@testset "recommended_dt damps the production Hall–Ohmic Fourier mode" begin
    for T in (Float32, Float64), evolving in (false, true), safety in (0.8, 1.0)
        g = FourierGrid((33,), (T(2π),))
        K = π / Float64(g.dx[1])
        wave_bound = K * (sqrt(K^2 + 4) + K) / 2
        η = (2.785 / 2.8) * wave_bound / K^2
        closure = evolving ? TransportElectrons(T(5 / 3)) : IsothermalElectrons(one(T))
        st = HybridStepper(g, HybridModel(closure; η), CIC(), 0)
        f = st.fields
        fill!(st.fn, one(T))
        fill!(f.ninv, one(T))
        fill!(f.pe, one(T))
        fill!(f.B[1], one(T))
        x = T.(0:32) .* g.dx[1]
        @. f.B[2] = T(0.001) * cos(T(16) * x)
        @. f.B[3] = T(0.001) * sin(T(16) * x)
        initial = sum(abs2, f.B[2]) + sum(abs2, f.B[3])
        dt = T(recommended_dt(g; η, safety))
        if evolving
            HybridPlasmaPIC._rk4_Bp!(st, dt, 1)
        else
            HybridPlasmaPIC._rk4_B!(st, dt)
        end
        measured = sqrt((sum(abs2, f.B[2]) + sum(abs2, f.B[3])) / initial)
        z = complex(-Float64(T(η)) * 16^2 * dt, Float64(16^2 * dt))
        expected = abs(1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24)
        @test isapprox(measured, expected; rtol = 128eps(T))
        @test measured <= 1 + 128eps(T)
    end
end
