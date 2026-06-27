# spacecraft.jl — synthetic spacecraft diagnostics for 1-D shock runs.
#
# A virtual probe samples a periodic grid field at an arbitrary (possibly
# moving) position via CIC linear interpolation, exactly the transpose-consistent
# gather used by deposit.jl: nodes sit at (i−1)·dx, i=1..n, periodic, and a point
# at fractional cell position s = x/dx interpolates linearly between nodes
# floor(s) and floor(s)+1 with weights (1−frac, frac).
#
# On top of the probe this file provides the standard shock-physics frame
# transforms (shock-frame velocity boost, de Hoffmann–Teller frame velocity) and
# a reflected-ion classifier, the building blocks for analysing collisionless
# shock crossings.
#
# Conventions (documented for the classifier): the upstream plasma lies at +x,
# and the shock front itself propagates in the +x direction at speed Vs. In the
# shock rest frame a particle velocity vx maps to (vx − Vs). A particle counts as
# "reflected" when it sits upstream of the front (x > x_shock) AND in the shock
# frame it is travelling back toward the upstream region (vx − Vs > 0).

# ---------------------------------------------------------------- gather

"""
    gather_at(field::Array{T,1}, g::FourierGrid{1,T}, xpos::Real) -> T

CIC (linear) interpolation of the periodic 1-D grid `field` at physical position
`xpos`. Grid node `i` (1-based) is located at `(i-1)·dx`; the value is the linear
blend of the two bracketing nodes with periodic wraparound, matching the gather
stencil in `deposit.jl`.
"""
function gather_at(field::Array{T,1}, g::FourierGrid{1,T}, xpos::Real) where {T}
    n = g.n[1]
    length(field) == n || throw(DimensionMismatch("field length must match the 1-D grid size"))
    dx = g.dx[1]
    s = T(xpos) / dx                 # fractional cell position
    i0 = floor(Int, s)               # 0-based base node
    f = s - i0                        # ∈ [0,1)
    ia = mod(i0, n) + 1              # 1-based, wrapped
    ib = mod(i0 + 1, n) + 1
    @inbounds return (one(T) - f) * field[ia] + f * field[ib]
end

# ---------------------------------------------------------------- probe

"""
    SyntheticProbe(x0; T=Float64)

A virtual spacecraft at position `x0` (stored as `Float64`) recording a time
series of `(t, value)` samples. Use [`sample!`](@ref) to record the local field
value and [`advance!`](@ref) to move the probe.
"""
mutable struct SyntheticProbe{T<:AbstractFloat}
    x::Float64
    t::Vector{Float64}
    val::Vector{T}
end

SyntheticProbe(x0::Real; T::Type{<:AbstractFloat} = Float64) =
    SyntheticProbe{T}(Float64(x0), Float64[], T[])

"""
    sample!(probe, field, g, t)

Record `(t, gather_at(field, g, probe.x))` into the probe's time series and
return the sampled value.
"""
function sample!(
    probe::SyntheticProbe{T},
    field::Array{T,1},
    g::FourierGrid{1,T},
    t::Real,
) where {T}
    v = gather_at(field, g, probe.x)
    push!(probe.t, Float64(t))
    push!(probe.val, v)
    return v
end

"""
    advance!(probe, vx, dt)

Move the probe by `vx·dt` along x (a moving spacecraft). Returns the new
position.
"""
function advance!(probe::SyntheticProbe, vx::Real, dt::Real)
    probe.x += Float64(vx) * Float64(dt)
    return probe.x
end

# ---------------------------------------------------------------- frame transforms

"""
    shock_frame(vx, Vs)

Boost an x-velocity into a frame moving at speed `Vs` along x: returns `vx - Vs`.
"""
shock_frame(vx::Real, Vs::Real) = vx - Vs

"""
    dehoffmann_teller_velocity(u::NTuple{3}, B::NTuple{3}) -> NTuple{3}

de Hoffmann–Teller frame velocity `V_HT = B × (u × B) / |B|²`, i.e. the
component of the bulk flow `u` perpendicular to `B`. In the dHT frame the
transformed perpendicular flow vanishes, so the residual flow `u - V_HT` is
parallel to `B` and the motional electric field `E = -(u - V_HT) × B` vanishes.

If `|B|² = 0` the frame is undefined; this returns the zero vector.
"""
function dehoffmann_teller_velocity(u::NTuple{3,<:Real}, B::NTuple{3,<:Real})
    Bx, By, Bz = B
    B2 = Bx * Bx + By * By + Bz * Bz
    if B2 == 0
        return (zero(B2), zero(B2), zero(B2))
    end
    ux, uy, uz = u
    # w = u × B
    wx = uy * Bz - uz * By
    wy = uz * Bx - ux * Bz
    wz = ux * By - uy * Bx
    # B × w
    vx = By * wz - Bz * wy
    vy = Bz * wx - Bx * wz
    vz = Bx * wy - By * wx
    return (vx / B2, vy / B2, vz / B2)
end

# ---------------------------------------------------------------- reflected ions

"""
    classify_reflected(ps::ParticleSet{1}, x_shock, frame_Vs) -> Vector{Bool}

Flag each particle as reflected at a 1-D shock. With the upstream plasma at +x
and the shock front moving in +x at `frame_Vs`, a particle is reflected when it
is upstream of the front (`x > x_shock`) AND, in the shock frame, it travels back
toward upstream (`vx - frame_Vs > 0`).
"""
function classify_reflected(ps::ParticleSet{1,T}, x_shock::Real, frame_Vs::Real) where {T}
    x = ps.x[1]
    vx = ps.v[1]
    N = nparticles(ps)
    out = Vector{Bool}(undef, N)
    xs = T(x_shock)
    Vs = T(frame_Vs)
    @inbounds for p = 1:N
        out[p] = (x[p] > xs) && (shock_frame(vx[p], Vs) > zero(T))
    end
    return out
end

# ---------------------------------------------------------------- four-spacecraft timing

"""
    crossing_time(ts, vals, level) -> Float64

First linear-interpolated time at which the time series `vals` crosses `level`
(either direction). Returns `NaN` if no crossing exists. The building block for
multi-spacecraft timing of a boundary crossing.
"""
function _require_finite_real_sequence(name::AbstractString, xs)
    @inbounds for i in eachindex(xs)
        x = xs[i]
        x isa Real || throw(ArgumentError("$(name) must contain real values"))
        isfinite(x) || throw(ArgumentError("$(name) must contain only finite values"))
    end
    return xs
end

@inline function _require_finite_point3(name::AbstractString, p, ::Type{T}) where {T<:AbstractFloat}
    return (
        _require_finite_real("$(name)[1]", p[1], T),
        _require_finite_real("$(name)[2]", p[2], T),
        _require_finite_real("$(name)[3]", p[3], T),
    )
end

function crossing_time(ts::AbstractVector, vals::AbstractVector, level::Real)
    length(ts) == length(vals) || throw(DimensionMismatch("ts and vals must have the same length"))
    _require_finite_real_sequence("ts", ts)
    _require_finite_real_sequence("vals", vals)
    level = _require_finite_real("level", level, Float64)
    n = length(vals)
    @inbounds for i = 1:n-1
        a = vals[i] - level
        b = vals[i+1] - level
        a == 0 && return Float64(ts[i])
        if a * b < 0
            f = a / (a - b)
            return Float64(ts[i]) + f * (Float64(ts[i+1]) - Float64(ts[i]))
        end
    end
    return n > 0 && vals[end] == level ? Float64(ts[end]) : NaN
end

"""
    four_spacecraft_timing(positions::NTuple{4,NTuple{3}}, times::NTuple{4})
        -> (; normal, speed, slowness)

Standard four-spacecraft timing analysis (e.g. Cluster/MMS). Given four probe
positions `rᵢ` and the times `tᵢ` at which a planar boundary `n̂·r = V t` crosses
each, solve `m·(rᵢ−r₀) = tᵢ−t₀` (i=1,2,3) for the slowness vector `m = n̂/V`, then
return the unit normal `n̂ = m/|m|`, the phase speed `V = 1/|m|`, and `m`.

The crossing times must be non-degenerate (the four probes not coplanar with the
boundary motion), else the 3×3 system is singular and the result is non-finite.
"""
@inline _det3(a11, a12, a13, a21, a22, a23, a31, a32, a33) =
    a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31)

function four_spacecraft_timing(positions::NTuple{4,NTuple{3,<:Real}}, times::NTuple{4,<:Real})
    pos = ntuple(i -> _require_finite_point3("positions[$i]", positions[i], Float64), 4)
    ts = ntuple(i -> _require_finite_real("times[$i]", times[i], Float64), 4)
    r0 = pos[1]
    t0 = ts[1]
    r1 = pos[2]
    r2 = pos[3]
    r3 = pos[4]
    a11 = Float64(r1[1]) - Float64(r0[1])
    a12 = Float64(r1[2]) - Float64(r0[2])
    a13 = Float64(r1[3]) - Float64(r0[3])
    a21 = Float64(r2[1]) - Float64(r0[1])
    a22 = Float64(r2[2]) - Float64(r0[2])
    a23 = Float64(r2[3]) - Float64(r0[3])
    a31 = Float64(r3[1]) - Float64(r0[1])
    a32 = Float64(r3[2]) - Float64(r0[2])
    a33 = Float64(r3[3]) - Float64(r0[3])
    b1 = ts[2] - t0
    b2 = ts[3] - t0
    b3 = ts[4] - t0
    det = _det3(a11, a12, a13, a21, a22, a23, a31, a32, a33)
    if det == 0.0
        return (; normal = (NaN, NaN, NaN), speed = NaN, slowness = (NaN, NaN, NaN))
    end
    mx = _det3(b1, a12, a13, b2, a22, a23, b3, a32, a33) / det
    my = _det3(a11, b1, a13, a21, b2, a23, a31, b3, a33) / det
    mz = _det3(a11, a12, b1, a21, a22, b2, a31, a32, b3) / det
    m2 = mx * mx + my * my + mz * mz
    if !(m2 > 0.0 && isfinite(m2))
        return (; normal = (NaN, NaN, NaN), speed = NaN, slowness = (NaN, NaN, NaN))
    end
    sp = 1.0 / sqrt(m2)
    n̂ = (mx * sp, my * sp, mz * sp)
    return (; normal = n̂, speed = sp, slowness = (mx, my, mz))
end

"""
    four_spacecraft_traces(; MA=3.0, probes, level=nothing, ...)
        -> (; traces, times, crossings, normal, speed)

Run a 3-D perpendicular shock while recording the perpendicular field `B_z` at
four virtual spacecraft `probes` (a `NTuple{4,NTuple{3}}` of (x,y,z)
positions), then recover the shock normal and speed by
[`four_spacecraft_timing`](@ref) of the `B_z` half-rise crossings. The physical
controls `Te`, `γe`, `vthi`, `η`, `db_turb`, and `field_method` match
[`run_perp_shock3d`](@ref). For the planar perpendicular shock the recovered
normal is ≈ x̂. `level` defaults to the midpoint between B0 and the downstream
peak.
"""
function four_spacecraft_traces(;
    MA::Real = 3.0,
    nx::Integer = 48,
    ny::Integer = 8,
    nz::Integer = 8,
    Lx::Real = 70.0,
    Ly::Real = 10.0,
    Lz::Real = 10.0,
    nppc::Integer = 8,
    nsteps::Integer = 420,
    dt::Real = 0.03,
    seed::Integer = 1,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    db_turb::Real = 0.0,
    field_method::Symbol = :rk4,
    probes::NTuple{4,NTuple{3,<:Real}} = (
        (5.0, 1.0, 1.0),
        (8.0, 8.0, 1.0),
        (8.0, 1.0, 8.0),
        (11.0, 4.0, 4.0),
    ),
    level::Union{Nothing,Real} = nothing,
)
    nsteps >= 1 || throw(ArgumentError("nsteps must be positive"))
    T = Float64
    _require_valid_positive_shock_ma(MA, T)
    B0 = one(T)
    probesT = ntuple(i -> _require_finite_point3("probes[$i]", probes[i], T), 4)
    sh, ps = _load_shock3d(; MA, nx, ny, nz, Lx, Ly, Lz, Te, γe, vthi, η, nppc, seed, db_turb)

    traces = ntuple(_ -> Float64[], 4)
    times = Float64[]
    dx = sh.sbp.dx
    for st = 1:nsteps
        step_shock3d!(sh, ps, T(dt); NB = 2, field_method = field_method)
        push!(times, st * T(dt))
        for q = 1:4
            xp, yp, zp = probesT[q]
            v = _gather3d(sh.B[3], xp, yp, zp, dx, sh.dy, sh.dz, nx, ny, nz)
            push!(traces[q], v)
        end
    end
    lvl = level === nothing ? (B0 + maximum(maximum, traces)) / 2 : _require_finite_real("level", level, T)
    crossings = ntuple(q -> crossing_time(times, traces[q], lvl), 4)
    res =
        all(isfinite, crossings) ? four_spacecraft_timing(probesT, crossings) :
        (; normal = (NaN, NaN, NaN), speed = NaN, slowness = (NaN, NaN, NaN))
    return (; traces, times, crossings, normal = res.normal, speed = res.speed)
end
