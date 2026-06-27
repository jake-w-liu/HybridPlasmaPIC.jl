# PhaseSpace.jl — velocity & phase-space histograms (from diagnostics.jl)

"""
    velocity_histogram(ps, comp; nbins=64, vmin, vmax) -> (centers, counts)

Weighted 1-D velocity distribution of component `comp` (1,2,3). `sum(counts)` ≈
total weight inside [vmin,vmax].
"""
function velocity_histogram(
    ps::ParticleSet{D,T},
    comp::Int;
    nbins::Int = 64,
    vmin = nothing,
    vmax = nothing,
) where {D,T}
    nbins > 0 || throw(ArgumentError("nbins must be positive"))
    v = ps.v[comp]
    w = ps.weight
    lo = vmin === nothing ? minimum(v) : T(vmin)
    hi = vmax === nothing ? maximum(v) : T(vmax)
    hi <= lo && (hi = lo + one(T))           # degenerate range ⇒ avoid /0 (NaN→crash)
    dv = (hi - lo) / nbins
    counts = zeros(T, nbins)
    @inbounds for p in eachindex(v)
        b = floor(Int, (v[p] - lo) / dv) + 1
        b == nbins + 1 && v[p] <= hi && (b = nbins)   # include the top edge (v==hi)
        1 <= b <= nbins && (counts[b] += w[p])
    end
    centers = [lo + (i - T(0.5)) * dv for i = 1:nbins]
    return centers, counts
end

"""
    phase_space_histogram(ps, sdim, vcomp; nx=64, nv=64, ...) -> (xc, vc, counts)

Weighted 2-D (x_{sdim}, v_{vcomp}) phase-space histogram.
"""
function phase_space_histogram(
    ps::ParticleSet{D,T},
    sdim::Int,
    vcomp::Int;
    nx::Int = 64,
    nv::Int = 64,
    xmin = nothing,
    xmax = nothing,
    vmin = nothing,
    vmax = nothing,
) where {D,T}
    nx > 0 || throw(ArgumentError("nx must be positive"))
    nv > 0 || throw(ArgumentError("nv must be positive"))
    x = ps.x[sdim]
    v = ps.v[vcomp]
    w = ps.weight
    xlo = xmin === nothing ? minimum(x) : T(xmin)
    xhi = xmax === nothing ? maximum(x) : T(xmax)
    vlo = vmin === nothing ? minimum(v) : T(vmin)
    vhi = vmax === nothing ? maximum(v) : T(vmax)
    xhi <= xlo && (xhi = xlo + one(T))       # degenerate range ⇒ avoid /0 (NaN→crash)
    vhi <= vlo && (vhi = vlo + one(T))
    dx = (xhi - xlo) / nx
    dv = (vhi - vlo) / nv
    counts = zeros(T, nx, nv)
    @inbounds for p in eachindex(x)
        i = floor(Int, (x[p] - xlo) / dx) + 1
        j = floor(Int, (v[p] - vlo) / dv) + 1
        i == nx + 1 && x[p] <= xhi && (i = nx)        # include top edges
        j == nv + 1 && v[p] <= vhi && (j = nv)
        (1 <= i <= nx && 1 <= j <= nv) && (counts[i, j] += w[p])
    end
    xc = [xlo + (i - T(0.5)) * dx for i = 1:nx]
    vc = [vlo + (j - T(0.5)) * dv for j = 1:nv]
    return xc, vc, counts
end

# ---------------------------------------------------------------- spectra
