# WKB ray tracing through a stratified hybrid plasma: trace a fast
# magnetosonic wave packet through a 1-D density wave and check the exact
# local invariant k(x) = ω / c(x), c² = v_A²(x) + c_s²(x).
#
#   julia --project=. examples/ray_tracing.jl

using HybridPlasmaPIC, Printf

L = 20.0                                        # box length [d_i]
nfun = x -> 1.0 + 0.3 * sin(2π * x / L)         # density wave [n0]
Te, γe = 0.08, 2.0                              # polytropic electrons

med = AnalyticRayMedium((x, y, z) -> nfun(x), (x, y, z) -> (0.0, 0.0, 1.0); Te, γe)
ray = trace_ray(med, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 0.01, nsteps = 1500)

@printf(
    "fast-branch ray:  ω = %.4f Ω_ci   status = %s   steps = %d\n",
    ray.ω,
    ray.status,
    length(ray.t) - 1
)
@printf("max |D| residual along ray: %.2e (integration accuracy monitor)\n", maximum(ray.residual))

worst = 0.0
for m = 1:50:length(ray.t)
    nx = nfun(ray.x[1, m])
    c = sqrt(1 / nx + γe * Te * nx^(γe - 1))    # local fast speed [v_A]
    global worst = max(worst, abs(ray.k[1, m] - ray.ω / c) / (ray.ω / c))
end
@printf("k(x) = ω/c(x) satisfied to %.2e (relative, along the whole ray)\n", worst)

# same medium sampled on a grid (as a snapshot of a simulation would be)
g = SpectralOperators.FourierGrid((256,), (L,))
xg = [(i - 1) * L / 256 for i = 1:256]
gmed = GridRayMedium(g, nfun.(xg), (zeros(256), zeros(256), fill(1.0, 256)); Te, γe)
gray = trace_ray(gmed, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 0.01, nsteps = 1500)
@printf(
    "grid-medium ray endpoint differs from analytic by %.2e d_i\n",
    abs(gray.x[1, end] - ray.x[1, end])
)
