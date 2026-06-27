# Boris.jl — Boris pusher (from particles.jl)

@inline function boris_kick(
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    qm::T,
    dt::T,
) where {T}
    h = qm * dt / 2
    vmx = vx + h * Ex
    vmy = vy + h * Ey
    vmz = vz + h * Ez
    tx = h * Bx
    ty = h * By
    tz = h * Bz
    t2 = tx * tx + ty * ty + tz * tz
    f = 2 / (1 + t2)
    sx = f * tx
    sy = f * ty
    sz = f * tz
    # v' = v⁻ + v⁻ × t
    vpx = vmx + (vmy * tz - vmz * ty)
    vpy = vmy + (vmz * tx - vmx * tz)
    vpz = vmz + (vmx * ty - vmy * tx)
    # v⁺ = v⁻ + v' × s
    vnx = vmx + (vpy * sz - vpz * sy)
    vny = vmy + (vpz * sx - vpx * sz)
    vnz = vmz + (vpx * sy - vpy * sx)
    return (vnx + h * Ex, vny + h * Ey, vnz + h * Ez)
end

"""
    push_uniform!(ps, E::NTuple{3}, B::NTuple{3}, dt)

Advance every particle one Boris leapfrog step in a spatially uniform field:
v^{n-1/2} → v^{n+1/2}, then x^{n} → x^{n+1} = x^n + dt·v^{n+1/2} (only the D
spatial coordinates move). Returns `ps`.
"""
function push_uniform!(ps::ParticleSet{D,T}, E::NTuple{3}, B::NTuple{3}, dt::Real) where {D,T}
    qm = ps.q / ps.m
    dtT = T(dt)
    Ex, Ey, Ez = T(E[1]), T(E[2]), T(E[3])
    Bx, By, Bz = T(B[1]), T(B[2]), T(B[3])
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = boris_kick(vx[p], vy[p], vz[p], Ex, Ey, Ez, Bx, By, Bz, qm, dtT)
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
    end
    @inbounds for d = 1:D
        xd = ps.x[d]
        vd = ps.v[d]
        for p in eachindex(xd)
            xd[p] += dtT * vd[p]
        end
    end
    return ps
end

function _check_gathered_field_lengths(ps::ParticleSet{D}, E, B, xmid) where {D}
    N = nparticles(ps)
    all(length(E[c]) == N for c = 1:3) ||
        throw(DimensionMismatch("E component length must match particle count"))
    all(length(B[c]) == N for c = 1:3) ||
        throw(DimensionMismatch("B component length must match particle count"))
    if xmid !== nothing
        length(xmid) == D ||
            throw(DimensionMismatch("xmid must have one component per spatial dimension"))
        all(length(xmid[d]) == N for d = 1:D) ||
            throw(DimensionMismatch("xmid component length must match particle count"))
    end
    return nothing
end

"""
    push_gathered!(ps, E, B, dt; xmid=nothing)

Advance particles one Boris leapfrog step using per-particle gathered electric
and magnetic fields. `E` and `B` are 3-tuples of particle-length vectors. If
`xmid` is provided, it receives `x^n + (dt/2) v^{n+1/2}` before positions are
advanced to `x^{n+1}`.
"""
function push_gathered!(
    ps::ParticleSet{D,T},
    E::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    B::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    dt::Real;
    xmid = nothing,
) where {D,T}
    _check_gathered_field_lengths(ps, E, B, xmid)
    qm = ps.q / ps.m
    dtT = T(dt)
    h = dtT / T(2)
    vx, vy, vz = ps.v
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = boris_kick(
            vx[p],
            vy[p],
            vz[p],
            T(Ex[p]),
            T(Ey[p]),
            T(Ez[p]),
            T(Bx[p]),
            T(By[p]),
            T(Bz[p]),
            qm,
            dtT,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            xmid !== nothing && (xmid[d][p] = ps.x[d][p] + h * ps.v[d][p])
            ps.x[d][p] += dtT * ps.v[d][p]
        end
    end
    return ps
end

function _push_uniform_broadcast!(
    ps::ParticleSet{D,T},
    E::NTuple{3},
    B::NTuple{3},
    dt::Real,
) where {D,T}
    return push_uniform!(ps, E, B, dt)
end

# ---------------------------------------------------------------- boundaries
