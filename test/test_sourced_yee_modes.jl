using HybridPlasmaPIC, Test
# Independent modal ODE: e'=c²*kappa*b-j, b'=-kappa*e.
# With zero initial fields and constant j, e=-j*sin(c*kappa*t)/(c*kappa),
# b=j*(1-cos(c*kappa*t))/(c²*kappa). No package stencil computes this reference.
function forced_error(n, steps, mode, periodic; continuum = false)
    g = YeeGrid((n,), (1.0,); boundaries = (periodic ? :periodic : :pec,))
    m = MaxwellVacuum(g; c = 1.25)
    h = 0.25 / steps
    j = 0.375
    @assert h < maxwell_cfl(m)
    k = (periodic ? 2pi : pi) * mode
    kappa = continuum ? k : 2sin(k * g.dx[1] / 2) / g.dx[1]
    current = map(a -> zero(a), m.fields.E)
    x = yee_coordinates(g, :electric, 2)[1]
    current[2] .= h * j * sin.(k * x)
    if !periodic
        current[2][1] = current[2][end] = 0.0
    end
    for _ = 1:steps
        step_maxwell!(m, h, current)
    end
    e = -j * sin(m.c * kappa * 0.25) / (m.c * kappa)
    b = j * (1 - cos(m.c * kappa * 0.25)) / (m.c^2 * kappa)
    ey = e * sin.(k * x)
    bz = b * cos.(k * yee_coordinates(g, :magnetic, 3)[1])
    return max(maximum(abs, m.fields.E[2] - ey), maximum(abs, m.fields.B[3] - bz))
end
@testset "Prescribed-current modal time and space order" begin
    for periodic in (false, true), mode in (1, 2, 3)
        errors = [forced_error(16, steps, mode, periodic) for steps in (64, 128, 256)]
        orders = log2.(errors[1:end-1] ./ errors[2:end])
        println((gate = :time, periodic = periodic, mode = mode, errors = errors, orders = orders))
        @test all(>(1.95), orders)
        errors = [forced_error(n, 32n, mode, periodic; continuum = true) for n in (12, 24, 48)]
        orders = log2.(errors[1:end-1] ./ errors[2:end])
        println((gate = :space, periodic = periodic, mode = mode, errors = errors, orders = orders))
        @test all(>(1.8), orders)
    end
end
