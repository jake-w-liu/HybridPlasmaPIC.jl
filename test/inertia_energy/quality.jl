module InertiaEnergyQualityTests
include("support.jl")
import JET
public_inertia(ps, B, n, c, g, d) = energy_budget(ps, B, n, c, g; de2 = d)
@testset "Inertia and shared exact square sum concrete compiler calls" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        J = ntuple(_ -> ones(T, g.n), 3)
        for fields in (J, map(a -> view(a, ntuple(_ -> Colon(), D)...), J))
            for (f, args) in (
                (H._field_square_sum_exact, (fields,)),
                (H._electron_inertia_energy, (fields, T(0.25), g.dx)),
            )
                f(args...)
                report = JET.report_call(f, Tuple{Core.Typeof.(args)...})
                if !isempty(JET.get_reports(report))
                    show(stdout, MIME"text/plain"(), report)
                    println()
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end
@testset "Call-local inertia arithmetic and independent-grid public concurrency" begin
    for T in (Float32, Float64), profile in (:ordinary, :overflow, :underflow)
        e =
            profile === :ordinary ? 0 :
            (T === Float32 ? 80 : 600) * (profile === :overflow ? 1 : -1)
        ce =
            profile === :ordinary ? 0 :
            (T === Float32 ? 100 : 800) * (profile === :overflow ? -1 : 1)
        J = (fill(ldexp(one(T), e), 4), zeros(T, 4), zeros(T, 4))
        de2 = ldexp(one(T), ce)
        dx = (one(T),)
        expected = H._electron_inertia_energy(J, de2, dx)
        before = map(copy, J)
        values = Vector{T}(undef, 64)
        Threads.@threads for i in eachindex(values)
            values[i] = H._electron_inertia_energy(J, de2, dx)
        end
        @test all(isequal(v, expected) for v in values)
        @test all(isequal.(J, before))
        # FourierGrid owns mutable FFT scratch: each independent public call gets
        # its own grid, while the physical fields and particles are shared read-only.
        grids = [FourierGrid((4,), (T(2pi),)) for _ = 1:32]
        B = (zeros(T, 4), zeros(T, 4), ldexp(one(T), e) .* T[1, 0, -1, 0])
        n = ones(T, 4)
        ps = ParticleSet{1,T}(0)
        c = IsothermalElectrons(zero(T))
        expected = energy_budget(ps, B, n, c, first(grids); de2 = de2)
        before = deepcopy((B, n, ps.x, ps.v, ps.weight))
        outputs = Vector{typeof(expected)}(undef, length(grids))
        Threads.@threads for i in eachindex(grids)
            outputs[i] = energy_budget(ps, B, n, c, grids[i]; de2 = de2)
        end
        for output in outputs
            @test isequal(output, expected)
        end
        @test isequal((B, n, ps.x, ps.v, ps.weight), before)
    end
end
@testset "Public inertia budgets retain concrete inference" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        B = ntuple(_ -> ones(T, g.n), 3)
        n = ones(T, g.n)
        c = IsothermalElectrons(T(0.5))
        ps = ParticleSet{D,T}(0)
        for population in (ps, [ps])
            args = (population, B, n, c, g, T(0.25))
            public_inertia(args...)
            report = JET.report_call(public_inertia, Tuple{Core.Typeof.(args)...})
            if !isempty(JET.get_reports(report))
                show(stdout, MIME"text/plain"(), report)
                println()
            end
            @test isempty(JET.get_reports(report))
        end
    end
end

end
