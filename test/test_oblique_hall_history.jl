module ObliqueHallHistoryTests
using HybridPlasmaPIC, LinearAlgebra, Test
include(joinpath(@__DIR__, "..", "validation", "oblique_hall_history.jl"))
include(joinpath(@__DIR__, "..", "validation", "mode_history.jl"))

# Independently derived in (density,ux,uy,uz,by,bz), unlike the oracle's
# (ux,uy,uz,bx,by,bz,pressure) frequency matrix.
function reduced_generator(k, angle, pressure)
    c, s = cosd(angle), sind(angle)
    A = zeros(ComplexF64, 6, 6)
    A[1, 2] = -im * k
    A[2, 1], A[2, 5] = -im * k * pressure, -im * k * s
    A[3, 5], A[4, 6] = im * k * c, im * k * c
    A[5, 2], A[5, 3], A[5, 6] = -im * k * s, im * k * c, -k^2 * c
    A[6, 4], A[6, 5] = im * k * c, k^2 * c
    return A
end

@testset "Independent oblique matrix, branch identities and normalization" begin
    for k in (0.7, 1.0, 2.0),
        angle in (0.0, 15.0, 30.0, 60.0, 75.0),
        (Ti, Te) in ((0.0, 0.0), (0.0, 0.2), (0.3, 0.2)),
        direction in (-1, 1)

        pressure = Ti + Te
        branches =
            pressure == 0 ? (:whistler, :ioncyclotron) : (:whistler, :ioncyclotron, :acoustic)
        for branch in branches
            ref = @inferred _oblique_hall_reference(k, angle, Ti, Te, branch, direction)
            v = ref.vector[[1, 2, 3, 4, 6, 7]]
            A = reduced_generator(k, angle, pressure)
            @test norm(A * v - ref.rate * v) <= 1e-11norm(A)
            @test ref.omega * direction > 0
            @test isapprox(norm(ref.vector), 1; atol = 1e-14)
            @test abs(ref.vector[5]) < 1e-12
            @test abs(imag(ref.vector[ref.component])) < 1e-14
            @test real(ref.vector[ref.component]) > 0
            # Independent scalar dispersion polynomial, including acoustic limit.
            x = ref.omega^2
            first =
                (x - k^2 * cosd(angle)^2) *
                (x^2 - k^2 * (1 + pressure) * x + k^4 * pressure * cosd(angle)^2)
            second = k^4 * cosd(angle)^2 * x * (x - k^2 * pressure)
            @test abs(first - second) <= 1e-10max(1, abs(first), abs(second))
            @test all(row -> row.overlap >= 0.8, ref.trace)
            @test all(row -> row.angle == 0 || row.assignment_gap > 1e-3, ref.trace)
            @test maximum(row.angle for row in ref.trace) == angle
        end
    end
    for branch in (:whistler, :ioncyclotron, :acoustic)
        ion = _oblique_hall_reference(1.0, 45.0, 0.5, 0.0, branch, 1)
        electron = _oblique_hall_reference(1.0, 45.0, 0.0, 0.5, branch, 1)
        @test ion.omega == electron.omega
        @test ion.vector == electron.vector
    end
    for branch in (:whistler, :ioncyclotron)
        ref = _oblique_hall_reference(1.0, 0.0, 0.0, 0.0, branch, 1)
        expected = branch == :whistler ? (sqrt(5) + 1) / 2 : (sqrt(5) - 1) / 2
        @test isapprox(ref.omega, expected; rtol = 1e-14)
    end
    # Warm acoustic/electromagnetic degeneracy must not choose an arbitrary basis.
    for speed in ((sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2)
        @test_throws ArgumentError _oblique_hall_reference(1.0, 30.0, speed^2, 0.0, :whistler, 1)
    end
    @test_throws ArgumentError _oblique_hall_reference(1.0, 30.0, 0.0, 0.0, :acoustic, 1)
    @test_throws ArgumentError _oblique_assignment(ones(ComplexF64, 7, 2), ones(ComplexF64, 7, 2))
end

function rhs_coefficient(k, angle, Ti, Te, column, amplitude)
    grid = FourierGrid((16,), (2pi / k,))
    state = HallMHDState(grid, HallMHDModel(IsothermalElectrons(Te); Ti))
    state.fields.n .= 1
    state.fields.B[1] .= cosd(angle)
    state.fields.B[2] .= sind(angle)
    fields = (state.fields.n, state.fields.ui..., state.fields.B[2], state.fields.B[3])
    for j = 1:16
        fields[column][j] += amplitude * cos(k * (j - 1) * grid.dx[1] + 0.37)
    end
    rhs = hall_mhd_rhs!(state)
    arrays = (rhs.dn, rhs.du..., rhs.dB[2], rhs.dB[3])
    return ComplexF64[2mode_amplitude(a, grid, (1,)) / 16 for a in arrays]
end

@testset "Central derivative of actual production Hall RHS" begin
    for k in (1.0, 2.0),
        angle in (23.0, 61.0),
        (Ti, Te) in ((0.0, 0.0), (0.0, 0.5), (0.3, 0.2)),
        column = 1:6

        A = reduced_generator(k, angle, Ti + Te)
        for amplitude in (1e-3, 5e-4)
            plus = rhs_coefficient(k, angle, Ti, Te, column, amplitude)
            minus = rhs_coefficient(k, angle, Ti, Te, column, -amplitude)
            measured = (plus - minus) / (2amplitude)
            expected = A[:, column] * cis(0.37)
            @test norm(measured - expected) <= 2amplitude^2 * max(1, norm(expected)) + 1e-10
        end
    end
end

@testset "Actual initialized fields, production histories and errors" begin
    for angle in (30.0, 60.0),
        (Ti, Te) in ((0.0, 0.0), (0.3, 0.2)),
        direction in (-1, 1),
        branch in
        (Ti + Te == 0 ? (:whistler, :ioncyclotron) : (:whistler, :ioncyclotron, :acoustic))

        state, ref = _oblique_hall_state(16, 1, angle, Ti, Te, branch, direction, 1e-5)
        initial = _oblique_mode_vector(state, 1)
        @test norm(initial - 1e-5ref.vector) < 1e-15
        rhs = hall_mhd_rhs!(state)
        fields = (rhs.dn, rhs.du..., rhs.dB...)
        measured = ComplexF64[2mode_amplitude(f, state.g, (1,)) / 16 for f in fields]
        @test norm(measured - ref.rate * initial) / norm(initial) < 1e-8
        history = _oblique_hall_history(; branch, angle, Ti, Te, direction, duration = 0.08)
        @test history.samples[:, 1] == initial
        @test length(history.times) == size(history.samples, 2) == 9
        @test history.samples[:, end] != initial
        @test isapprox(history.times[end], 0.08; atol = 1e-15)
        @test history.mass_error < 1e-14
        @test history.longitudinal_error < 1e-14
        @test history.floor_count == 0
        fit = _fit_mode_history(vec(history.samples[ref.component, :]), history.dt)
        @test isapprox(fit.omega, ref.omega; rtol = 1e-6)
        @test abs(fit.gamma) < 1e-7
        @test _oblique_history_errors(history).history_error < 1e-7
    end
    basic = (; branch = :whistler, angle = 30.0, Ti = 0.0, Te = 0.0, duration = 0.08)
    for override in (
        (; n = 4),
        (; mode = 0),
        (; mode = 8),
        (; dt = 0.0),
        (; dt = Inf),
        (; duration = Inf),
        (; duration = 0.085),
        (; amplitude = 0.0),
        (; amplitude = 0.002),
        (; angle = 90.0),
        (; angle = NaN),
        (; Ti = -1.0),
        (; Te = Inf),
        (; direction = 0),
        (; branch = :unknown),
    )
        @test_throws ArgumentError _oblique_hall_history(; merge(basic, override)...)
    end
    mktempdir() do dir
        path = joinpath(dir, "raw.csv")
        history = _oblique_hall_history(; basic..., output = path)
        lines = readlines(path)
        @test length(lines) == 10
        @test length(split(first(lines), ',')) == 15
        @test parse.(Float64, split(last(lines), ','))[1] == history.times[end]
        @test_throws ArgumentError _oblique_hall_history(; basic..., output = path)
        @test readlines(path) == lines
    end
    state, ref = _oblique_hall_state(16, 1, 30.0, 0.0, 0.0, :whistler, 1, 1e-5)
    state.fields.n .= 0.0
    io = IOBuffer()
    @test_throws ArgumentError _sample_oblique_hall!(state, ref, 1, 0.01, 8, io)
    @test length(split(chomp(String(take!(io))), '\n')) == 2
    static = _fit_mode_history(fill(1.0 + 0im, 16), 0.01)
    @test static.periods == 0
end
end

module ObliqueHallCaseTests
using Test, TOML
include(joinpath(@__DIR__, "..", "validation", "35_production_oblique_hall", "run.jl"))
@testset "Frozen oblique matrix and visible run failures" begin
    profiles = _oblique_hall_profiles()
    orders = _oblique_order_profiles()
    @test length(profiles) == 60
    @test length(unique(p.label for p in profiles)) == 60
    @test length(orders) == 12
    @test all(p -> p.duration == 80 && p.amplitude == 1e-5 && p.mode == 1, profiles)
    @test Set(p.angle for p in profiles) == Set((30.0, 60.0))
    @test all(p -> p.branch != :acoustic || p.Ti + p.Te > 0, profiles)
    @test all(p -> p.n == (p.refinement == :grid ? 32 : 16), profiles)
    @test all(p -> p.dt == (p.refinement == :time ? 0.005 : 0.01), profiles)
    for p in profiles
        matched = filter(
            q ->
                (q.angle, q.Ti, q.Te, q.branch, q.direction) ==
                (p.angle, p.Ti, p.Te, p.branch, p.direction),
            profiles,
        )
        @test Set(q.refinement for q in matched) == Set((:baseline, :grid, :time))
    end
    for interval in (0.0, 100.0)
        mktempdir() do dir
            p = merge(
                first(profiles),
                (; label = "failure_interval$(interval)", dt = interval, duration = 800.0),
            )
            @test_throws ArgumentError _run_oblique_hall_case(dir, [p], orders)
            receipt = TOML.parsefile(joinpath(dir, "35_production_oblique_hall_provenance.toml"))
            @test receipt["status"] == "error"
            @test receipt["failure_type"] == "ArgumentError"
            @test !isempty(receipt["failure_cause"])
            @test isempty(receipt["completed_profiles"])
            @test receipt["active_profile"] == p.label
            if interval > 0
                raw = joinpath(dir, "35_production_oblique_hall_$(p.label).csv")
                @test length(readlines(raw)) >= 2
                @test haskey(receipt["artifact_sha256"], basename(raw))
            end
            before = read(joinpath(dir, "35_production_oblique_hall_provenance.toml"))
            @test_throws ArgumentError _run_oblique_hall_case(dir, [p], orders)
            @test read(joinpath(dir, "35_production_oblique_hall_provenance.toml")) == before
        end
    end
end
end
