# Audit: recommended_dt ignores eta/etaH. With hyper-resistivity etaH the field
# RHS has a REAL decay rate etaH*k^4 that exceeds the whistler omega ~ k^2 by
# etaH*k^2; at recommended_dt the RK4 subcycle then sits far outside the RK4
# real-axis stability interval (~2.79) and blows up. 1D, frozen moments n=1, ui=0.
using HybridPlasmaPIC, Random
const H = HybridPlasmaPIC
T = Float64
n = 64
g = FourierGrid((n,), (T(2pi),))
etaH = 0.01
model = HybridModel(IsothermalElectrons(0.0); ηH = etaH)
st = HybridStepper(g, model, CIC(), 4)
fill!(st.fn, one(T))
for c = 1:3
    fill!(st.fui[c], zero(T))
end
f = st.fields
H._ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, model.closure, T(model.nfloor), f.floor_count, g)

dt = recommended_dt(g; NB = 1, integrator = :rk4)
kmax = pi / g.dx[1]
println("dt_rec = $dt ; etaH*k^4*dt at k=31: $(etaH*31.0^4*dt)  (RK4 real-axis limit 2.79)")

rng = MersenneTwister(3)
fill!(f.B[1], one(T))
f.B[2] .= 1e-8 .* randn(rng, n)
f.B[3] .= 1e-8 .* randn(rng, n)
H.project_divfree!(f.B, g)
for i = 1:200
    H._rk4_B!(st, dt)
    if !all(all(isfinite, f.B[c]) for c = 1:3)
        println("NON-FINITE at substep $i with etaH=$etaH at recommended_dt")
        exit()
    end
end
println("stayed finite; max|B2| = $(maximum(abs, f.B[2]))")
