# Temporal convergence of HybridStepper step! (NB=1 vs NB=4) and CAM-CL step_camcl!
# (NB=2), on a smooth deterministic magnetized 1D problem. Error at fixed t_final
# against a fine-dt reference; expect observed order ~2 in all cases.
using HybridPlasmaPIC, Random
const H = HybridPlasmaPIC
T = Float64

function setup(n, L, nppc)
    g = FourierGrid((n,), (T(L),))
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    k = 2pi / L
    for p = 1:N
        xp = ps.x[1][p]
        ps.v[1][p] = 0.03 * sin(k * xp)
        ps.v[2][p] = 0.05 * cos(k * xp)
        ps.v[3][p] = 0.0
    end
    return g, ps, k
end

function runcase(kind::Symbol, dt, nsteps, NB)
    g, ps, k = setup(16, 2pi, 100)
    model = HybridModel(IsothermalElectrons(0.5))
    x = [(i - 1) * g.dx[1] for i = 1:16]
    if kind === :hybrid
        st = HybridStepper(g, model, CIC(), nparticles(ps))
        fill!(st.fields.B[1], 1.0)
        st.fields.B[2] .= 0.02 .* cos.(k .* x)
        init!(st, ps)
        for _ = 1:nsteps
            step!(st, ps, dt; NB = NB)
        end
        return st.fields.B, ps
    else
        st = CAMCLStepper(g, model, CIC(), nparticles(ps))
        fill!(st.fields.B[1], 1.0)
        st.fields.B[2] .= 0.02 .* cos.(k .* x)
        init_camcl!(st, ps)
        for _ = 1:nsteps
            step_camcl!(st, ps, dt; NB = NB)
        end
        return st.fields.B, ps
    end
end

# compare only INTEGER-level quantities (B at n, x at n): particle v lives at the
# half level t_f - dt/2, which differs between runs by O(dt) purely by convention.
function err(Ba, psa, Bb, psb)
    e = 0.0
    for c = 1:3
        e += sum(abs2, Ba[c] .- Bb[c])
    end
    e += sum(abs2, psa.x[1] .- psb.x[1]) / length(psa.x[1])
    sqrt(e)
end

tf = 0.32
for (kind, NB) in ((:hybrid, 1), (:hybrid, 4), (:camcl, 2), (:camcl, 4))
    Bref, psref = runcase(kind, tf / 512, 512, max(NB, kind === :camcl ? 2 : 1))
    prev = NaN
    print("$kind NB=$NB: ")
    for m in (8, 16, 32, 64)
        B, ps = runcase(kind, tf / m, m, NB)
        e = err(B, ps, Bref, psref)
        if !isnan(prev)
            print(" order=", round(log2(prev / e), digits = 2))
        end
        prev = e
    end
    println()
end
