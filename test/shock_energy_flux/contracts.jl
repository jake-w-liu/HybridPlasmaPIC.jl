module ShockFluxContractsTests
using HybridPlasmaPIC, Test, Random
import JET
include("legacy.jl")

const H = HybridPlasmaPIC

function population(::Type{T}, n = 8) where {T}
    ps = ParticleSet{1,T}(n)
    ps.x[1] .= T.(0:n-1) ./ 2
    ps.weight .= T(0.75)
    ps.v[1] .= T(1.25)
    ps.v[2] .= T(-0.5)
    ps.v[3] .= T(0.25)
    ps
end

function view_population(ps)
    part(a) = view(vcat(a, a), 1:2:length(a)*2)
    # Values need not match the original population; both owners see the same
    # independent strided storage and are checked for mutation.
    ParticleSet{1,eltype(ps.weight)}(
        map(part, ps.x),
        map(part, ps.v),
        part(ps.weight),
        part(ps.id),
        part(ps.tag),
        ps.q,
        ps.m,
    )
end

function failure(f, args...)
    try
        f(args...)
        return nothing
    catch err
        return (typeof(err), sprint(showerror, err))
    end
end

println((version = VERSION, threads = Threads.nthreads(), JET = pkgversion(JET)))
flush(stdout)
@testset "Shock flux compatibility, storage and validation ownership" begin
    for T in (Float16, Float32, Float64, BigFloat)
        sh = PerpShock(8, T(7))
        ps = population(T)
        for p in (ps, view_population(ps), population(T, 0))
            before = deepcopy((p.x, p.v, p.weight, sh))
            @test isequal(H._ion_energy_flux(sh, p), Legacy._ion_energy_flux(sh, p))
            @test isequal(
                H.boundary_energy_flux(sh; ps = p),
                Legacy.boundary_energy_flux(sh; ps = p),
            )
            @test isequal((p.x, p.v, p.weight), before[1:3])
            @test all(
                isequal(getfield(sh, key), getfield(before[4], key)) for
                key in fieldnames(typeof(sh)) if key != :s
            )
            @test isequal(sh.s.H, before[4].s.H)
        end
        for gamma in (one(T), T(Inf), T(NaN)), pressure in (zero(T), one(T)), u in (zero(T), one(T))
            sh.γe = gamma
            sh.pe .= pressure
            sh.ux .= u
            @test isequal(H.boundary_energy_flux(sh), Legacy.boundary_energy_flux(sh))
        end
        sh = PerpShock(8, T(7))
        for field in (:n, :ux, :uy, :pe, :Ey, :Bz), value in (T(Inf), T(-Inf), T(NaN))
            changed = deepcopy(sh)
            getfield(changed, field)[1] = value
            @test isequal(H.boundary_energy_flux(changed), Legacy.boundary_energy_flux(changed))
        end
        # Compare the first observed exception, including text, with the owner.
        for invalid in (:position, :vx, :vy, :vz, :weight, :all)
            p = population(T)
            invalid in (:position, :all) && (p.x[1][1] = T(NaN))
            invalid in (:vx, :all) && (p.v[1][1] = T(Inf))
            invalid in (:vy, :all) && (p.v[2][1] = T(NaN))
            invalid in (:vz, :all) && (p.v[3][1] = T(-Inf))
            invalid in (:weight, :all) && (p.weight[1] = T(NaN))
            old = failure(Legacy._ion_energy_flux, sh, p)
            new = failure(H._ion_energy_flux, sh, p)
            @test old !== nothing
            @test isequal(old, new)
        end
        for norm in (zero(T), T(-1), T(Inf), T(NaN))
            altered = deepcopy(sh)
            altered.s.H[1] = norm
            @test isequal(H._ion_energy_flux(altered, ps), Legacy._ion_energy_flux(altered, ps))
        end
    end
end

@testset "Shock flux concrete compiler analysis" begin
    for T in (Float32, Float64)
        sh = PerpShock(8, T(7))
        ps = population(T)
        for p in (ps, view_population(ps))
            calls = (
                (H._ion_energy_flux, (sh, p)),
                (H.boundary_energy_flux, (sh,)),
                (H._shock_bulk_energy_flux, (T(1.2), T(-1.1), T(0.3))),
                (H._shock_electron_enthalpy_flux, (T(1.5), T(0.3), T(-1.1))),
                (H._shock_flux_factor_quantum, (T(1.25),)),
                (H._shock_flux_quantum, (sh, p)),
            )
            for (f, args) in calls
                report = JET.report_call(f, Tuple{Core.Typeof.(args)...})
                if !isempty(JET.get_reports(report))
                    show(stdout, MIME"text/plain"(), report)
                    println()
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
    for (f, args) in (
        (H._shock_flux_float32_product_error, (1.0f0,)),
        (H._shock_flux_float32_sum_error, (1.0f0,)),
    )
        report = JET.report_call(f, Tuple{Core.Typeof.(args)...})
        if !isempty(JET.get_reports(report))
            show(stdout, MIME"text/plain"(), report)
            println()
        end
        @test isempty(JET.get_reports(report))
    end
end

@testset "Concurrent shared-input ordinary and complete shock flux" begin
    for T in (Float32, Float64), cold in (false, true)
        sh = PerpShock(8, T(7))
        ps = population(T)
        if cold
            huge = ldexp(one(T), T === Float32 ? 70 : 600)
            ps.v[1][1] = huge
            ps.v[1][2] = -huge
            ps.x[1][1] = ps.x[1][2] = zero(T)
        end
        before = deepcopy((sh.s.H, ps.x, ps.v, ps.weight))
        expected = H._ion_energy_flux(sh, ps)
        output = Vector{Vector{T}}(undef, 64)
        Threads.@threads for i in eachindex(output)
            output[i] = H._ion_energy_flux(sh, ps)
        end
        for value in output
            @test isequal(value, expected)
        end
        @test isequal((sh.s.H, ps.x, ps.v, ps.weight), before)
        @test all(output[i] !== output[j] for i in eachindex(output) for j = 1:i-1)
    end
end
println("SHOCK_FLUX_CONTRACTS_COMPLETE")
flush(stdout)

end
