module ShockTotalQualityTests
using HybridPlasmaPIC, Test
import JET
const H = HybridPlasmaPIC
public_particles(sh, ps) = boundary_energy_flux(sh; ps = ps)
@testset "Complete boundary total concrete compiler calls" begin
    for T in (Float32, Float64)
        sh = PerpShock(8, T(7))
        ps = ParticleSet{1,T}(2)
        ps.x[1] .= T[0, 7]
        ps.v[1] .= one(T)
        part(a) = view(vcat(a, a), 1:2:2length(a))
        strided = ParticleSet{1,T}(
            map(part, ps.x),
            map(part, ps.v),
            part(ps.weight),
            part(ps.id),
            part(ps.tag),
            ps.q,
            ps.m,
        )
        for p in (ps, strided)
            calls = (
                (public_particles, (sh, p)),
                (H._shock_ion_flux_state, (sh, p)),
                (H._shock_total_flux, (sh, p, 1, T(1), T(2), T(3), 0.0, nothing)),
                (H._shock_total_flux, (sh, nothing, 1, T(1), T(2), T(3), 0.0, nothing)),
                (H._shock_zero_boundary_ion_flux, (sh, p, 1)),
                (H._shock_exact_boundary_ion_flux, (sh, p, 1)),
                (H._shock_exact_bulk_flux, (T(1), T(2), T(3))),
                (H._shock_exact_enthalpy_flux, (T(1.5), T(2), T(3))),
                (H._shock_exact_particle_flux, (T(1), T(2), T(3), T(4))),
            )
            for (f, args) in calls
                f(args...)
                report =
                    JET.report_call(f, Tuple{Core.Typeof.(args)...}; stacktrace_types_limit = 0)
                if !isempty(JET.get_reports(report))
                    show(stdout, MIME"text/plain"(), report)
                    println()
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end
end
