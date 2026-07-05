# Exact per-particle-step multiplier of _cl_subcycle_B! for the scalar mode
# db/dt = lambda*b (lambda = -i*omega whistler, lambda = -nu resistive/hyper-resistive),
# reproducing the implemented algorithm: half-step bootstrap, NB staggered kicks,
# final half-step sync, average. Also RK4 for comparison.
function cl_mult(lam_dt::Complex{Float64}, NB::Int)
    h = lam_dt / NB
    A = 1.0 + 0im
    C = 1.0 + 0im
    C += (h / 2) * A                 # bootstrap F(A)
    for k = 0:NB-1
        A += h * C                   # F(C)
        if k < NB - 1
            C += h * A               # F(A)
        end
    end
    C += (h / 2) * A                 # final sync F(A)
    return 0.5 * (A + C)
end

rk4_R(z::Complex{Float64}) = 1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24

println("== whistler (imaginary axis): |m| per step, z = omega*h per substep ==")
for NB in (2, 4, 8)
    bad = 0.0
    zbad = 0.0
    for z = 0.005:0.005:2.0
        m = abs(cl_mult(-im * z * NB, NB))
        if m > 1 + 1e-12 && m > bad
            bad = m
            zbad = z
        end
    end
    if bad > 0
        # growth over 1000 particle steps at worst z
        println("NB=$NB: max |m| = $bad at omega*h = $zbad ; |m|^1000 = $(bad^1000)")
    else
        println("NB=$NB: |m| <= 1 for all omega*h <= 2")
    end
end
println("\n== whistler at the recommended_dt operating point omega*h = 0.8*2 = 1.6 ==")
for NB in (2, 4, 8)
    m = abs(cl_mult(-im * 1.6 * NB, NB))
    println("NB=$NB: |m| = $m ; per-omega_ci-time growth over 1000 steps = $(m^1000)")
end

println("\n== resistive damping (real axis): m vs exact exp(-nu*dt), z = nu*h ==")
for NB in (2, 4, 8)
    for zh in (0.5, 1.0, 1.5, 2.0)
        m = cl_mult(complex(-zh * NB, 0.0), NB)
        ex = exp(-zh * NB)
        println("NB=$NB nu*h=$zh: m = $(real(m))  exact = $ex  |m|>1: $(abs(m)>1)")
    end
end

println("\n== RK4 real-axis limit (for eta/etaH in the RK4 subcycle) ==")
for z in (2.0, 2.5, 2.785, 2.8, 3.0)
    println("z=nu*h=$z: |R(-z)| = $(abs(rk4_R(complex(-z,0.0))))")
end
