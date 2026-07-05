# Audit: is recommended_dt's kmax = pi/min(dx) the actual stiffest whistler mode
# in 2D? Frozen-moment B-subcycle (exactly what _rk4_B! integrates), uniform B0
# along the in-plane diagonal, n=1, u_i=0. The electron-MHD whistler for this
# operator is omega = (k.b)|k| (normalized), maximized at the k-space CORNER
# |k|^2 = sum_d (pi/dx_d)^2 -- a factor D larger than (pi/min dx)^2 on a square grid.
using HybridPlasmaPIC, Random
const H = HybridPlasmaPIC
T = Float64
n = 32
g = FourierGrid((n, n), (T(2pi), T(2pi)))
model = HybridModel(IsothermalElectrons(0.0))
st = HybridStepper(g, model, CIC(), 4)

fill!(st.fn, one(T))
for c = 1:3
    fill!(st.fui[c], zero(T))
end
f = st.fields
H._ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, model.closure, T(model.nfloor), f.floor_count, g)

function runcase(dt, nsteps; label = "")
    rng = MersenneTwister(7)
    b0 = (1 / sqrt(2), 1 / sqrt(2), 0.0)
    for c = 1:3
        fill!(f.B[c], b0[c])
        f.B[c] .+= 1e-8 .* randn(rng, size(f.B[c]))
    end
    H.project_divfree!(f.B, g)
    dev0 = maximum(abs, f.B[3])
    ok = true
    for i = 1:nsteps
        H._rk4_B!(st, dt)
        if !all(all(isfinite, f.B[c]) for c = 1:3)
            println("$label: NON-FINITE at substep $i")
            ok = false
            break
        end
    end
    dev1 = ok ? maximum(abs, f.B[3]) : NaN
    println("$label dt=$(round(dt, sigdigits=4))  dev0=$(round(dev0, sigdigits=3)) -> dev=$(dev1)")
end

kmax1d = pi / g.dx[1]
K = kmax1d
omega_formula = 0.5 * (sqrt(K^4 + 4K^2) + K^2)
# grid's largest representable mode along each axis is m = n/2 (Nyquist) or n/2-1;
# use the corner mode just inside Nyquist: m=(15,15) => k=(15,15), |k|=15*sqrt(2)
kc = 15.0 * sqrt(2)
omega_corner = kc^2          # k parallel to B0 diagonal: omega = (k.b)|k| = |k|^2
println("omega used by recommended_dt = $omega_formula ; true corner-mode omega ≈ $omega_corner")

dt_rec = recommended_dt(g; NB = 1, integrator = :rk4)
println("recommended dt = $dt_rec ; omega_corner*dt = $(omega_corner*dt_rec)  (RK4 limit 2.83)")
runcase(dt_rec, 2000; label = "recommended_dt   ")
dt_safe = 0.8 * 2.8 / omega_corner
runcase(dt_safe, 2000; label = "corner-corrected ")
