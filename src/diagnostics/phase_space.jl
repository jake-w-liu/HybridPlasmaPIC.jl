# PhaseSpace.jl — velocity & phase-space histograms (from diagnostics.jl)

@inline function _require_finite_hist_value(name::AbstractString, value)
    isfinite(value) || throw(ArgumentError("$(name) must be finite"))
    return value
end

function _check_velocity_component(comp::Int)
    1 <= comp <= 3 || throw(ArgumentError("velocity component must be in 1:3, got $comp"))
    return nothing
end

function _check_spatial_component(sdim::Int, ::Val{D}) where {D}
    1 <= sdim <= D || throw(ArgumentError("spatial dimension must be in 1:$D, got $sdim"))
    return nothing
end

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
    _check_velocity_component(comp)
    v = ps.v[comp]
    w = ps.weight
    if isempty(v) && (vmin === nothing || vmax === nothing)
        throw(ArgumentError("velocity bounds must be provided for empty particle sets"))
    end
    lo = _require_finite_hist_value("velocity lower bound", vmin === nothing ? minimum(v) : T(vmin))
    hi = _require_finite_hist_value("velocity upper bound", vmax === nothing ? maximum(v) : T(vmax))
    hi <= lo && (hi = lo + one(T))           # degenerate range ⇒ avoid /0 (NaN→crash)
    dv = (hi - lo) / nbins
    counts = zeros(T, nbins)
    @inbounds for p in eachindex(v)
        vp = _require_finite_hist_value("velocity[$p]", v[p])
        _require_finite_hist_value("weight[$p]", w[p])
        b = floor(Int, (vp - lo) / dv) + 1
        b == nbins + 1 && vp <= hi && (b = nbins)     # include the top edge (v==hi)
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
    _check_spatial_component(sdim, Val(D))
    _check_velocity_component(vcomp)
    x = ps.x[sdim]
    v = ps.v[vcomp]
    w = ps.weight
    if isempty(x) && (xmin === nothing || xmax === nothing || vmin === nothing || vmax === nothing)
        throw(ArgumentError("position and velocity bounds must be provided for empty particle sets"))
    end
    xlo =
        _require_finite_hist_value("position lower bound", xmin === nothing ? minimum(x) : T(xmin))
    xhi =
        _require_finite_hist_value("position upper bound", xmax === nothing ? maximum(x) : T(xmax))
    vlo =
        _require_finite_hist_value("velocity lower bound", vmin === nothing ? minimum(v) : T(vmin))
    vhi =
        _require_finite_hist_value("velocity upper bound", vmax === nothing ? maximum(v) : T(vmax))
    xhi <= xlo && (xhi = xlo + one(T))       # degenerate range ⇒ avoid /0 (NaN→crash)
    vhi <= vlo && (vhi = vlo + one(T))
    dx = (xhi - xlo) / nx
    dv = (vhi - vlo) / nv
    counts = zeros(T, nx, nv)
    @inbounds for p in eachindex(x)
        xp = _require_finite_hist_value("position[$p]", x[p])
        vp = _require_finite_hist_value("velocity[$p]", v[p])
        _require_finite_hist_value("weight[$p]", w[p])
        i = floor(Int, (xp - xlo) / dx) + 1
        j = floor(Int, (vp - vlo) / dv) + 1
        i == nx + 1 && xp <= xhi && (i = nx)          # include top edges
        j == nv + 1 && vp <= vhi && (j = nv)
        (1 <= i <= nx && 1 <= j <= nv) && (counts[i, j] += w[p])
    end
    xc = [xlo + (i - T(0.5)) * dx for i = 1:nx]
    vc = [vlo + (j - T(0.5)) * dv for j = 1:nv]
    return xc, vc, counts
end

# ---------------------------------------------------------------- spectra
