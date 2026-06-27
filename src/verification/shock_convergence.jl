# shock_converge.jl — Phase 11 (1-D) collisionless perpendicular-shock
# Mach-sweep extension + numerical convergence study.
#
# Thin research drivers on top of the verified reflecting-wall PerpShock model
# (`run_perp_shock`, SHK-002). Two functions:
#
#   • mach_sweep(; MAs=[1.2,2.0,4.0,6.0], kwargs...)
#       — covers the full 1-D M_A = 1.2, 2, 4, 6 set: one run_perp_shock result
#         per upstream Alfvén Mach number, each augmented with its `MA`. The
#         supercritical (M_A ≥ 2) runs form a real, flux-frozen shock with
#         downstream compression 2 < n2 < 4.
#
#   • convergence_study(; MA=4.0, params...)
#       — runs the SAME physical shock at two grid resolutions (N) AND two
#         particles-per-cell (nppc), reporting the downstream compression n2 for
#         each. The shock statistics are converged when the n2 values agree
#         across both the Δx (= Lx/N) and the ppc refinements — i.e. the result
#         is insensitive to the numerical resolution, as a well-resolved kinetic
#         shock should be.
#
# Both are pure post-`run_perp_shock` orchestration (no new physics): they exist
# to make the Mach sweep + convergence demonstration a single, testable call.

"""
    mach_sweep(; MAs=[1.2, 2.0, 4.0, 6.0], N=512, Lx=120.0, Te=0.125, γe=5/3,
                 vthi=0.35, η=0.02, nppc=32, nsteps=600, seed=1)
        -> Vector{NamedTuple}

Run [`run_perp_shock`](@ref) once for each upstream Alfvén Mach number in `MAs`,
forwarding the shared keyword arguments, and return the per-`MA` diagnostic
NamedTuples (each augmented with its input `MA`). The default `MAs` covers the
full 1-D `M_A = 1.2, 2, 4, 6` perpendicular-shock set.

Each element is the full `run_perp_shock` result
`(; MA, n2, Bz2, Vs, X_rh, frozen_ratio, reflected_fraction, M_real, xf)`. The
supercritical runs (`M_A ≥ 2`) form a real, flux-frozen shock: the frozen-in
ratio `(Bz2/B0)/n2 ≈ 1` and the downstream compression sits in `2 < n2 < 4`
(above the weak limit, below the strong-shock fluid ceiling). `nsteps`/`nppc`/`N`
are kept modest by default so the whole sweep runs in a few seconds.
"""
function mach_sweep(;
    MAs = [1.2, 2.0, 4.0, 6.0],
    N::Integer = 512,
    Lx::Real = 120.0,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    nppc::Integer = 32,
    nsteps::Integer = 600,
    seed::Integer = 1,
)
    out = NamedTuple[]
    for MA in MAs
        r = run_perp_shock(;
            MA = MA,
            N = Int(N),
            Lx = Lx,
            Te = Te,
            γe = γe,
            vthi = vthi,
            η = η,
            nppc = Int(nppc),
            nsteps = Int(nsteps),
            seed = Int(seed),
        )
        push!(out, (; MA = Float64(MA), r...))
    end
    return out
end

"""
    convergence_study(; MA=4.0, Ns=(256, 512), ppcs=(32, 64), Lx=120.0,
                        Te=0.125, γe=5/3, vthi=0.35, η=0.02, nsteps=500, seed=1)
        -> (; MA, base, grid, ppc, n2_base, n2_fine_N, n2_fine_ppc,
              rel_grid, rel_ppc, rel_max, converged)

Numerical convergence study for the 1-D perpendicular hybrid shock at Alfvén
Mach number `MA`. Runs [`run_perp_shock`](@ref) at:

  • a **base** resolution `(N=Ns[1], nppc=ppcs[1])`,
  • a **grid-refined** run `(N=Ns[2], nppc=ppcs[1])` — halves Δx = Lx/N,
  • a **particle-refined** run `(N=Ns[1], nppc=ppcs[2])` — doubles ppc,

all with otherwise identical physical parameters, and reports the downstream
compression `n2` of each. The shock statistics are *converged* (insensitive to
the discretization) when both relative changes
`rel_grid = |n2_fine_N − n2_base|/n2_base` and
`rel_ppc  = |n2_fine_ppc − n2_base|/n2_base` are small; `converged` is `true`
when `rel_max = max(rel_grid, rel_ppc) ≤ tol` (the `tol` keyword, default 0.12 — the
caller decides, e.g. ~0.12). Each `*_run` field carries the full
`run_perp_shock` NamedTuple so the caller can inspect `frozen_ratio`, `Vs`, etc.
"""
function convergence_study(;
    MA::Real = 4.0,
    Ns = (256, 512),
    ppcs = (32, 64),
    Lx::Real = 120.0,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    nsteps::Integer = 500,
    seed::Integer = 1,
    tol::Real = 0.12,
)
    length(Ns) == 2 || throw(ArgumentError("Ns must hold exactly two resolutions"))
    length(ppcs) == 2 || throw(ArgumentError("ppcs must hold exactly two ppc values"))
    Nlo, Nhi = Int(Ns[1]), Int(Ns[2])
    plo, phi = Int(ppcs[1]), Int(ppcs[2])

    runshock(N, nppc) = run_perp_shock(;
        MA = MA,
        N = N,
        Lx = Lx,
        Te = Te,
        γe = γe,
        vthi = vthi,
        η = η,
        nppc = nppc,
        nsteps = Int(nsteps),
        seed = Int(seed),
    )

    base = runshock(Nlo, plo)        # reference: coarse grid, few particles
    grid = runshock(Nhi, plo)        # refine Δx (more nodes)
    ppc = runshock(Nlo, phi)         # refine particle sampling (more ppc)

    n2_base = base.n2
    n2_fine_N = grid.n2
    n2_fine_ppc = ppc.n2

    rel_grid = abs(n2_fine_N - n2_base) / abs(n2_base)
    rel_ppc = abs(n2_fine_ppc - n2_base) / abs(n2_base)
    rel_max = max(rel_grid, rel_ppc)

    return (;
        MA = Float64(MA),
        base = base,
        grid = grid,
        ppc = ppc,
        n2_base = n2_base,
        n2_fine_N = n2_fine_N,
        n2_fine_ppc = n2_fine_ppc,
        rel_grid = rel_grid,
        rel_ppc = rel_ppc,
        rel_max = rel_max,
        converged = rel_max <= tol,
    )
end
