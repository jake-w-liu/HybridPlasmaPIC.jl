module ParallelWaveHistoryTests
using HybridPlasmaPIC, LinearAlgebra, Test
include(joinpath(@__DIR__, "..", "validation", "parallel_wave_history.jl"))

@testset "Cold CIC eigenmode initialization and adjoint force" begin
    amplitude = 1e-4
    for n in (8, 16, 32), nppc in (4, 9, 16), mode in (1, 2), branch in (:whistler, :ioncyclotron)
        state, ps, continuum, discrete =
            _parallel_wave_state(:hybrid, branch, mode, n, nppc, 0.0, amplitude)
        grid = state.g
        U =
            (
                mode_amplitude(state.fields.ui[2], grid, (mode,)) +
                im * mode_amplitude(state.fields.ui[3], grid, (mode,))
            ) / n
        b =
            (
                mode_amplitude(state.fields.B[2], grid, (mode,)) +
                im * mode_amplitude(state.fields.B[3], grid, (mode,))
            ) / n
        @test isapprox(b, amplitude; rtol = 1e-12)
        @test isapprox(U / b, discrete.velocity_ratio; rtol = 1e-12)
        @test maximum(abs, state.fields.n .- 1) < 1e-12
        # On a midpoint lattice, mean(f*(1-f))=1/6+1/(12*nppc^2).
        theta = 2pi * mode / n
        q = 1 - 2(1 - cos(theta)) * (1 / 6 + 1 / (12nppc^2))
        @test isapprox(discrete.coupling, q; rtol = 1e-13)
        @test isapprox(discrete.rate * U, im * (q - 1) * U + im * mode * q * b; rtol = 1e-12)
        @test isapprox(discrete.rate * b, im * mode * U + im * mode^2 * b; rtol = 1e-12)
        # Gather actual production E and B and inspect the initial Lorentz force.
        # The derivative uses no timestep or Boris update approximation.
        E = ntuple(_ -> zeros(length(ps.weight)), 3)
        B = ntuple(_ -> zeros(length(ps.weight)), 3)
        gather_vector!(E, state.fields.E, ps, grid, CIC())
        gather_vector!(B, state.fields.B, ps, grid, CIC())
        for j in eachindex(ps.weight)
            acceleration = complex(
                E[2][j] + ps.v[3][j] * B[1][j] - ps.v[1][j] * B[3][j],
                E[3][j] + ps.v[1][j] * B[2][j] - ps.v[2][j] * B[1][j],
            )
            expected = discrete.rate * complex(ps.v[2][j], ps.v[3][j])
            @test isapprox(acceleration, expected; atol = 1e-12amplitude, rtol = 1e-11)
        end
        # Independent closed-form continuum branches select the same root signs.
        omega =
            branch == :whistler ? (sqrt(mode^4 + 4mode^2) + mode^2) / 2 :
            (mode^2 - sqrt(mode^4 + 4mode^2)) / 2
        @test isapprox(imag(continuum.rate), omega; rtol = 1e-13)
    end
end

@testset "Production histories and configuration failures" begin
    for solver in (:hall, :hybrid), branch in (:whistler, :ioncyclotron)
        h = _parallel_wave_history(;
            solver,
            branch,
            mode = 1,
            n = 8,
            nppc = 4,
            dt = 0.01,
            duration = 0.08,
        )
        @test length(h.samples) == length(h.times) == 9
        @test h.times[1] == 0
        @test isapprox(h.times[end], 0.08; rtol = 1e-14)
        @test all(isfinite, h.samples)
        @test h.samples[end] != h.samples[1]
        @test h.density_error < 1e-10
    end
    for kwargs in (
        (; solver = :unknown),
        (; branch = :unknown),
        (; mode = 0),
        (; mode = 8),
        (; dt = 0.0),
        (; dt = NaN),
        (; duration = Inf),
        (; duration = 0.085),
        (; nppc = 0),
        (; subcycles = 0),
        (; eta = -1.0),
        (; eta = 0.01),
        (; amplitude = 0.0),
    )
        config = merge((; solver = :hybrid, branch = :whistler, mode = 1), kwargs)
        @test_throws ArgumentError _parallel_wave_history(; config...)
    end
end
end # module
