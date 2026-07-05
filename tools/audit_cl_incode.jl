# In-code confirmation: does _cl_subcycle_B! (as implemented, restart+average each
# particle step) amplify a whistler mode at omega*h well below the claimed bound 2?
# 1D, B0 = x_hat, frozen n=1, u_i=0: electron-MHD whistler omega = k^2 (normalized).
# Track the Fourier amplitude of mode m over repeated subcycled steps and compare the
# measured per-step growth with the scalar-recursion prediction.
using HybridPlasmaPIC, Random
const H = HybridPlasmaPIC
T = Float64
n = 32
g = FourierGrid((n,), (T(2pi),))
model = HybridModel(IsothermalElectrons(0.0))
st = CAMCLStepper(g, model, CIC(), 4)

fill!(st.fn, one(T))
for c = 1:3
    fill!(st.fui[c], zero(T))
end
f = st.fields
H._ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, model.closure, T(model.nfloor), f.floor_count, g)

# scalar prediction (same recursion as _cl_subcycle_B!)
function cl_mult(lam_dt::Complex{Float64}, NB::Int)
    h = lam_dt / NB
    A = 1.0 + 0im
    C = 1.0 + 0im
    C += (h / 2) * A
    for k = 0:NB-1
        A += h * C
        k < NB - 1 && (C += h * A)
    end
    C += (h / 2) * A
    return 0.5 * (A + C)
end

m_mode = 15                       # k = 15: stiffest non-Nyquist mode, omega = 225
omega = Float64(m_mode^2)  # all other modes have smaller omega*h
NB = 4
for zh in (1.0, 1.6)
    dt = zh * NB / omega          # omega * (dt/NB) = zh
    rng = MersenneTwister(11)
    fill!(f.B[1], one(T))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    k = Float64(m_mode)
    # seed a right-circular whistler eigenmode at mode m
    f.B[2] .= 1e-8 .* cos.(k .* x)
    f.B[3] .= 1e-8 .* sin.(k .* x)
    a0 = abs(mode_amplitude(f.B[2], g, (m_mode,)))
    nst = 120
    for i = 1:nst
        H._cl_subcycle_B!(st, T(dt), NB)
    end
    a1 = abs(mode_amplitude(f.B[2], g, (m_mode,)))
    grow = (a1 / a0)^(1 / nst)
    pred = abs(cl_mult(-im * omega * dt, NB))
    println(
        "omega*h=$zh NB=$NB: measured per-step growth = $(round(grow, sigdigits=6)); scalar prediction |m| = $(round(pred, sigdigits=6))",
    )
end
