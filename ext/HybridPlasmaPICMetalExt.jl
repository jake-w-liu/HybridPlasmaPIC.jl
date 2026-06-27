module HybridPlasmaPICMetalExt

import HybridPlasmaPIC
import Metal

HybridPlasmaPIC.extension_dependency_module(::Val{:metal}) = Metal
HybridPlasmaPIC.extension_device_array_type(::Val{:metal}) = Metal.MtlArray

_metal_functional() =
    try
        !isdefined(Metal, :functional) || getproperty(Metal, :functional)()
    catch err
        false
    end

function HybridPlasmaPIC.backend_memory_status(::Val{:metal})
    return HybridPlasmaPIC.BackendMemoryStatus(
        :metal,
        _metal_functional(),
        false;
        note = "Metal.jl does not expose memory-pool counters",
    )
end

HybridPlasmaPIC.reclaim_backend_memory!(::Val{:metal}) = false

@inline function _shape_code(::HybridPlasmaPIC.NGP)
    return Int32(1)
end

@inline function _shape_code(::HybridPlasmaPIC.CIC)
    return Int32(2)
end

@inline function _shape_code(::HybridPlasmaPIC.TSC)
    return Int32(3)
end

@inline function _shape_width(code::Int32)
    return code == Int32(1) ? Int32(1) : (code == Int32(2) ? Int32(2) : Int32(3))
end

@inline function _round_nearest_even_int32(s::Float32)
    floored = floor(s)
    base = unsafe_trunc(Int32, floored)
    frac = s - floored
    half = 0.5f0
    if frac > half
        return base + Int32(1)
    elseif frac < half
        return base
    else
        return (base & Int32(1)) == Int32(0) ? base : base + Int32(1)
    end
end

@inline function _floor_int32(s::Float32)
    return unsafe_trunc(Int32, floor(s))
end

@inline function _stencil_base(code::Int32, s::Float32)
    if code == Int32(1)
        return _round_nearest_even_int32(s)
    elseif code == Int32(2)
        return _floor_int32(s)
    else
        return _round_nearest_even_int32(s) - Int32(1)
    end
end

@inline function _stencil_weight(code::Int32, s::Float32, base::Int32, offset::Int32)
    if code == Int32(1)
        return 1.0f0
    elseif code == Int32(2)
        f = s - base
        return offset == Int32(1) ? 1.0f0 - f : f
    else
        center = base + Int32(1)
        δ = s - center
        half = 0.5f0
        if offset == Int32(1)
            return half * (half - δ)^2
        elseif offset == Int32(2)
            return 0.75f0 - δ^2
        else
            return half * (half + δ)^2
        end
    end
end

@inline function _periodic_index(base::Int32, offset::Int32, n::Int32)
    r = rem(base + offset - Int32(1), n)
    r < Int32(0) && (r += n)
    return r + Int32(1)
end

function _gather1_kernel!(out, field, x1, n1::Int32, dx1::Float32, code::Int32, N::Int32)
    p = Int32(Metal.thread_position_in_grid().x)
    p > N && return nothing
    s1 = x1[p] / dx1
    base1 = _stencil_base(code, s1)
    width = _shape_width(code)
    acc = 0.0f0
    o1 = Int32(1)
    while o1 <= width
        w1 = _stencil_weight(code, s1, base1, o1)
        i1 = _periodic_index(base1, o1, n1)
        acc += field[i1] * w1
        o1 += Int32(1)
    end
    out[p] = acc
    return nothing
end

function _gather2_kernel!(
    out,
    field,
    x1,
    x2,
    n1::Int32,
    n2::Int32,
    dx1::Float32,
    dx2::Float32,
    code::Int32,
    N::Int32,
)
    p = Int32(Metal.thread_position_in_grid().x)
    p > N && return nothing
    s1 = x1[p] / dx1
    s2 = x2[p] / dx2
    base1 = _stencil_base(code, s1)
    base2 = _stencil_base(code, s2)
    width = _shape_width(code)
    acc = 0.0f0
    o2 = Int32(1)
    while o2 <= width
        w2 = _stencil_weight(code, s2, base2, o2)
        i2 = _periodic_index(base2, o2, n2)
        o1 = Int32(1)
        while o1 <= width
            w1 = _stencil_weight(code, s1, base1, o1)
            i1 = _periodic_index(base1, o1, n1)
            acc += field[i1, i2] * w1 * w2
            o1 += Int32(1)
        end
        o2 += Int32(1)
    end
    out[p] = acc
    return nothing
end

function _gather3_kernel!(
    out,
    field,
    x1,
    x2,
    x3,
    n1::Int32,
    n2::Int32,
    n3::Int32,
    dx1::Float32,
    dx2::Float32,
    dx3::Float32,
    code::Int32,
    N::Int32,
)
    p = Int32(Metal.thread_position_in_grid().x)
    p > N && return nothing
    s1 = x1[p] / dx1
    s2 = x2[p] / dx2
    s3 = x3[p] / dx3
    base1 = _stencil_base(code, s1)
    base2 = _stencil_base(code, s2)
    base3 = _stencil_base(code, s3)
    width = _shape_width(code)
    acc = 0.0f0
    o3 = Int32(1)
    while o3 <= width
        w3 = _stencil_weight(code, s3, base3, o3)
        i3 = _periodic_index(base3, o3, n3)
        o2 = Int32(1)
        while o2 <= width
            w2 = _stencil_weight(code, s2, base2, o2)
            i2 = _periodic_index(base2, o2, n2)
            o1 = Int32(1)
            while o1 <= width
                w1 = _stencil_weight(code, s1, base1, o1)
                i1 = _periodic_index(base1, o1, n1)
                acc += field[i1, i2, i3] * w1 * w2 * w3
                o1 += Int32(1)
            end
            o2 += Int32(1)
        end
        o3 += Int32(1)
    end
    out[p] = acc
    return nothing
end

@inline function _axis_weight_for_cell(code::Int32, s::Float32, cell::Int32, n::Int32)
    base = _stencil_base(code, s)
    width = _shape_width(code)
    weight = 0.0f0
    offset = Int32(1)
    while offset <= width
        if _periodic_index(base, offset, n) == cell
            weight += _stencil_weight(code, s, base, offset)
        end
        offset += Int32(1)
    end
    return weight
end

function _deposit1_kernel!(out, vals, x1, n1::Int32, dx1::Float32, code::Int32, Np::Int32)
    i1 = Int32(Metal.thread_position_in_grid().x)
    i1 > n1 && return nothing
    acc = 0.0f0
    p = Int32(1)
    while p <= Np
        s1 = x1[p] / dx1
        w1 = _axis_weight_for_cell(code, s1, i1, n1)
        acc += vals[p] * w1
        p += Int32(1)
    end
    out[i1] = acc
    return nothing
end

function _deposit2_kernel!(
    out,
    vals,
    x1,
    x2,
    n1::Int32,
    n2::Int32,
    dx1::Float32,
    dx2::Float32,
    code::Int32,
    Np::Int32,
    Ncells::Int32,
)
    q = Int32(Metal.thread_position_in_grid().x)
    q > Ncells && return nothing
    q0 = q - Int32(1)
    i1 = rem(q0, n1) + Int32(1)
    i2 = div(q0, n1) + Int32(1)
    acc = 0.0f0
    p = Int32(1)
    while p <= Np
        s1 = x1[p] / dx1
        s2 = x2[p] / dx2
        w1 = _axis_weight_for_cell(code, s1, i1, n1)
        if w1 != 0.0f0
            w2 = _axis_weight_for_cell(code, s2, i2, n2)
            acc += vals[p] * w1 * w2
        end
        p += Int32(1)
    end
    out[i1, i2] = acc
    return nothing
end

function _deposit3_kernel!(
    out,
    vals,
    x1,
    x2,
    x3,
    n1::Int32,
    n2::Int32,
    n3::Int32,
    dx1::Float32,
    dx2::Float32,
    dx3::Float32,
    code::Int32,
    Np::Int32,
    Ncells::Int32,
)
    q = Int32(Metal.thread_position_in_grid().x)
    q > Ncells && return nothing
    q0 = q - Int32(1)
    plane = n1 * n2
    i1 = rem(q0, n1) + Int32(1)
    i2 = rem(div(q0, n1), n2) + Int32(1)
    i3 = div(q0, plane) + Int32(1)
    acc = 0.0f0
    p = Int32(1)
    while p <= Np
        s1 = x1[p] / dx1
        s2 = x2[p] / dx2
        s3 = x3[p] / dx3
        w1 = _axis_weight_for_cell(code, s1, i1, n1)
        if w1 != 0.0f0
            w2 = _axis_weight_for_cell(code, s2, i2, n2)
            if w2 != 0.0f0
                w3 = _axis_weight_for_cell(code, s3, i3, n3)
                acc += vals[p] * w1 * w2 * w3
            end
        end
        p += Int32(1)
    end
    out[i1, i2, i3] = acc
    return nothing
end

function _metal_launch_config(N::Int)
    threads = min(max(N, 1), 256)
    groups = cld(N, threads)
    return threads, groups
end

function _check_metal_gather_args(out, field, ps, g)
    N = HybridPlasmaPIC.nparticles(ps)
    length(out) == N || throw(DimensionMismatch("out length must match particle count"))
    size(field) == g.n || throw(DimensionMismatch("field size must match grid size"))
    N <= typemax(Int32) ||
        throw(ArgumentError("Metal gather supports at most $(typemax(Int32)) particles"))
    all(n -> n <= typemax(Int32), g.n) ||
        throw(ArgumentError("Metal gather supports grid axes up to $(typemax(Int32)) cells"))
    return N
end

function _check_metal_deposit_args(out, vals, ps, g)
    Np = HybridPlasmaPIC.nparticles(ps)
    length(vals) == Np || throw(DimensionMismatch("vals length must match particle count"))
    size(out) == g.n || throw(DimensionMismatch("out size must match grid size"))
    Ncells = length(out)
    Np <= typemax(Int32) ||
        throw(ArgumentError("Metal deposition supports at most $(typemax(Int32)) particles"))
    Ncells <= typemax(Int32) ||
        throw(ArgumentError("Metal deposition supports at most $(typemax(Int32)) grid cells"))
    all(n -> n <= typemax(Int32), g.n) ||
        throw(ArgumentError("Metal deposition supports grid axes up to $(typemax(Int32)) cells"))
    return Np, Ncells
end

function HybridPlasmaPIC.deposit_scalar!(
    out::Metal.MtlArray{T,1},
    ps::HybridPlasmaPIC.ParticleSet{1,T,X,V,W,I,G},
    vals::Metal.MtlArray{T,1},
    g::HybridPlasmaPIC.FourierGrid{1,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    Np, Ncells = _check_metal_deposit_args(out, vals, ps, g)
    threads, groups = _metal_launch_config(Ncells)
    Metal.@metal threads = threads groups = groups _deposit1_kernel!(
        out,
        vals,
        ps.x[1],
        Int32(g.n[1]),
        g.dx[1],
        _shape_code(shape),
        Int32(Np),
    )
    return out
end

function HybridPlasmaPIC.deposit_scalar!(
    out::Metal.MtlArray{T,2},
    ps::HybridPlasmaPIC.ParticleSet{2,T,X,V,W,I,G},
    vals::Metal.MtlArray{T,1},
    g::HybridPlasmaPIC.FourierGrid{2,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    Np, Ncells = _check_metal_deposit_args(out, vals, ps, g)
    threads, groups = _metal_launch_config(Ncells)
    Metal.@metal threads = threads groups = groups _deposit2_kernel!(
        out,
        vals,
        ps.x[1],
        ps.x[2],
        Int32(g.n[1]),
        Int32(g.n[2]),
        g.dx[1],
        g.dx[2],
        _shape_code(shape),
        Int32(Np),
        Int32(Ncells),
    )
    return out
end

function HybridPlasmaPIC.deposit_scalar!(
    out::Metal.MtlArray{T,3},
    ps::HybridPlasmaPIC.ParticleSet{3,T,X,V,W,I,G},
    vals::Metal.MtlArray{T,1},
    g::HybridPlasmaPIC.FourierGrid{3,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    Np, Ncells = _check_metal_deposit_args(out, vals, ps, g)
    threads, groups = _metal_launch_config(Ncells)
    Metal.@metal threads = threads groups = groups _deposit3_kernel!(
        out,
        vals,
        ps.x[1],
        ps.x[2],
        ps.x[3],
        Int32(g.n[1]),
        Int32(g.n[2]),
        Int32(g.n[3]),
        g.dx[1],
        g.dx[2],
        g.dx[3],
        _shape_code(shape),
        Int32(Np),
        Int32(Ncells),
    )
    return out
end

function HybridPlasmaPIC.gather_scalar!(
    out::Metal.MtlArray{T,1},
    field::Metal.MtlArray{T,1},
    ps::HybridPlasmaPIC.ParticleSet{1,T,X,V,W,I,G},
    g::HybridPlasmaPIC.FourierGrid{1,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    N = _check_metal_gather_args(out, field, ps, g)
    N == 0 && return out
    threads, groups = _metal_launch_config(N)
    Metal.@metal threads = threads groups = groups _gather1_kernel!(
        out,
        field,
        ps.x[1],
        Int32(g.n[1]),
        g.dx[1],
        _shape_code(shape),
        Int32(N),
    )
    return out
end

function HybridPlasmaPIC.gather_scalar!(
    out::Metal.MtlArray{T,1},
    field::Metal.MtlArray{T,2},
    ps::HybridPlasmaPIC.ParticleSet{2,T,X,V,W,I,G},
    g::HybridPlasmaPIC.FourierGrid{2,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    N = _check_metal_gather_args(out, field, ps, g)
    N == 0 && return out
    threads, groups = _metal_launch_config(N)
    Metal.@metal threads = threads groups = groups _gather2_kernel!(
        out,
        field,
        ps.x[1],
        ps.x[2],
        Int32(g.n[1]),
        Int32(g.n[2]),
        g.dx[1],
        g.dx[2],
        _shape_code(shape),
        Int32(N),
    )
    return out
end

function HybridPlasmaPIC.gather_scalar!(
    out::Metal.MtlArray{T,1},
    field::Metal.MtlArray{T,3},
    ps::HybridPlasmaPIC.ParticleSet{3,T,X,V,W,I,G},
    g::HybridPlasmaPIC.FourierGrid{3,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    N = _check_metal_gather_args(out, field, ps, g)
    N == 0 && return out
    threads, groups = _metal_launch_config(N)
    Metal.@metal threads = threads groups = groups _gather3_kernel!(
        out,
        field,
        ps.x[1],
        ps.x[2],
        ps.x[3],
        Int32(g.n[1]),
        Int32(g.n[2]),
        Int32(g.n[3]),
        g.dx[1],
        g.dx[2],
        g.dx[3],
        _shape_code(shape),
        Int32(N),
    )
    return out
end

function HybridPlasmaPIC.gather_vector!(
    out::Tuple{<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1}},
    field::Tuple{<:Metal.MtlArray{T,D},<:Metal.MtlArray{T,D},<:Metal.MtlArray{T,D}},
    ps::HybridPlasmaPIC.ParticleSet{D,T,X,V,W,I,G},
    g::HybridPlasmaPIC.FourierGrid{D,T},
    shape::HybridPlasmaPIC.ShapeFunction,
) where {
    D,
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.gather_scalar!(out[1], field[1], ps, g, shape)
    HybridPlasmaPIC.gather_scalar!(out[2], field[2], ps, g, shape)
    HybridPlasmaPIC.gather_scalar!(out[3], field[3], ps, g, shape)
    return out
end

function HybridPlasmaPIC.push_uniform!(
    ps::HybridPlasmaPIC.ParticleSet{D,T,X,V,W,I,G},
    E::NTuple{3},
    B::NTuple{3},
    dt::Real,
) where {
    D,
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    return HybridPlasmaPIC._push_uniform_broadcast!(ps, E, B, dt)
end

function HybridPlasmaPIC.push_gathered!(
    ps::HybridPlasmaPIC.ParticleSet{D,T,X,V,W,I,G},
    E::Tuple{<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1}},
    B::Tuple{<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1},<:Metal.MtlArray{T,1}},
    dt::Real;
    xmid = nothing,
) where {
    D,
    T,
    X<:Metal.MtlArray{T,1},
    V<:Metal.MtlArray{T,1},
    W<:Metal.MtlArray{T,1},
    I<:Metal.MtlArray{UInt64,1},
    G<:Metal.MtlArray{UInt32,1},
}
    HybridPlasmaPIC.validate_particle_backend_eltype(Val(:metal), T)
    HybridPlasmaPIC._check_gathered_field_lengths(ps, E, B, xmid)
    if xmid !== nothing
        length(xmid) == D ||
            throw(DimensionMismatch("xmid must have one component per spatial dimension"))
        all(xmid[d] isa Metal.MtlArray{T,1} for d = 1:D) ||
            throw(ArgumentError("Metal gathered-field push requires Metal xmid storage"))
    end

    qm = ps.q / ps.m
    dtT = T(dt)
    hq = qm * dtT / T(2)
    hx = dtT / T(2)
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    vx, vy, vz = ps.v

    vmx = vx .+ hq .* Ex
    vmy = vy .+ hq .* Ey
    vmz = vz .+ hq .* Ez
    tx = hq .* Bx
    ty = hq .* By
    tz = hq .* Bz
    t2 = tx .* tx .+ ty .* ty .+ tz .* tz
    f = T(2) ./ (one(T) .+ t2)
    sx = f .* tx
    sy = f .* ty
    sz = f .* tz
    vpx = vmx .+ (vmy .* tz .- vmz .* ty)
    vpy = vmy .+ (vmz .* tx .- vmx .* tz)
    vpz = vmz .+ (vmx .* ty .- vmy .* tx)
    vnx = vmx .+ (vpy .* sz .- vpz .* sy)
    vny = vmy .+ (vpz .* sx .- vpx .* sz)
    vnz = vmz .+ (vpx .* sy .- vpy .* sx)
    vx .= vnx .+ hq .* Ex
    vy .= vny .+ hq .* Ey
    vz .= vnz .+ hq .* Ez

    for d = 1:D
        xmid !== nothing && (xmid[d] .= ps.x[d] .+ hx .* ps.v[d])
        ps.x[d] .+= dtT .* ps.v[d]
    end
    return ps
end

function HybridPlasmaPIC.validate_particle_backend_eltype(::Val{:metal}, ::Type{T}) where {T}
    T === Float32 ||
        throw(ArgumentError("Metal particle backend storage requires Float32; got $(T)"))
    return nothing
end

function HybridPlasmaPIC.disallow_scalar_indexing!(::Val{:metal})
    if isdefined(Metal, :allowscalar)
        getproperty(Metal, :allowscalar)(false)
    end
    return nothing
end

end # module
