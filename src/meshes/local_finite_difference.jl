# local_finite_difference.jl (§5.3 meshes/) — the SBP one-sided/central FD
# stencils live in SpectralOperators.jl (sbp_deriv!). Yee Maxwell differences
# below act between edges/faces, rather than on the collocated spectral fields.

function _yee_check_shape(fields, g::YeeGrid{D,T}, electric::Bool) where {D,T}
    length(fields) == 3 || throw(DimensionMismatch("a Yee vector has three components"))
    for c = 1:3
        fields[c] isa StridedArray{T,D} || throw(
            ArgumentError("Yee operators require CPU strided arrays matching the grid precision"),
        )
        size(fields[c]) == _yee_shape(g, electric, c) ||
            throw(DimensionMismatch("Yee component $c has the wrong shape"))
        Base.require_one_based_indexing(fields[c])
    end
    return nothing
end

@inline function _yee_difference_indices(
    I::CartesianIndex{D},
    axis::Int,
    g::YeeGrid{D},
    forward::Bool,
) where {D}
    axis > D && return nothing
    index = I[axis]
    if forward
        next = index == g.n[axis] && g.periodic[axis] ? 1 : index + 1
        J = CartesianIndex(ntuple(d -> d == axis ? next : I[d], D))
        return J, I
    end
    # Even magnetic extension at the PEC wall; tangential E stays constrained.
    !g.periodic[axis] && (index == 1 || index == g.n[axis] + 1) && return nothing
    previous = index == 1 ? g.n[axis] : index - 1
    J = CartesianIndex(ntuple(d -> d == axis ? previous : I[d], D))
    return I, J
end

@inline function _yee_difference(a, I, axis, g, forward)
    indices = _yee_difference_indices(I, axis, g, forward)
    indices === nothing && return zero(eltype(a))
    positive, negative = indices
    return (a[positive] - a[negative]) / g.dx[axis]
end

function _yee_difference_exact(a, I, axis, g, forward)
    indices = _yee_difference_indices(I, axis, g, forward)
    indices === nothing && return zero(Rational{BigInt})
    positive, negative = indices
    return (_exact_rational(a[positive]) - _exact_rational(a[negative])) /
           _exact_rational(g.dx[axis])
end

function _yee_curl_exact(fields, I, c, g, electric)
    j, k = mod1(c + 1, 3), mod1(c + 2, 3)
    return _yee_difference_exact(fields[k], I, j, g, electric) -
           _yee_difference_exact(fields[j], I, k, g, electric)
end

# ---------------------------------------------------------------- row engine
#
# Component c of a curl-type stencil evaluates `op(old, A, B)` at every
# entry: A is the difference of `sa` (component c+2) along axis c+1 and B the
# difference of `sb` (component c+1) along axis c+2. A forward difference maps
# edges to faces (electric source); a backward difference maps faces to edges
# and is absent on PEC walls, where tangential E is constrained. Rows run along
# axis 1 with fixed row offsets, so ordinary loops have unit stride and no
# per-entry index logic; periodic wraps and absent axis-1 differences at row
# ends are evaluated as separate entries. Absent differences (suppressed axes
# and walls) are passed as `nothing`, never as an arithmetic zero.

@inline _yee_row(a, J::CartesianIndex) = LinearIndices(a)[1, J] - 1
@inline _yee_row(a, J::CartesianIndex, axis::Int, index::Int) =
    LinearIndices(a)[1, CartesianIndex(Base.setindex(Tuple(J), index, axis - 1))] - 1

# Linear offsets of the (+, -) samples of one difference in row J, relative to
# the axis-1 index, and whether the difference exists in that row.
@inline function _yee_term_offsets(
    s,
    J::CartesianIndex,
    ::Val{axis},
    g::YeeGrid{D},
    ::Val{forward},
) where {axis,forward,D}
    axis > D && return (0, 0, false)
    if axis == 1
        row = _yee_row(s, J)
        return forward ? (row + 1, row, true) : (row, row - 1, true)
    end
    index, n = J[axis-1], g.n[axis]
    if forward
        next = index == n && g.periodic[axis] ? 1 : index + 1
        return (_yee_row(s, J, axis, next), _yee_row(s, J), true)
    end
    !g.periodic[axis] && (index == 1 || index == n + 1) && return (0, 0, false)
    return (_yee_row(s, J), _yee_row(s, J, axis, index == 1 ? n : index - 1), true)
end

@inline _yee_load(::Nothing, k) = nothing
@inline _yee_load(a, k) = @inbounds a[k]
@inline _yee_pair(s, plus, minus, i, ::Val{true}) = @inbounds (s[plus+i], s[minus+i])
@inline _yee_pair(s, plus, minus, i, ::Val{false}) = nothing

@inline _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, hasA::Val, hasB::Val) =
    op(_yee_load(old, base + i), _yee_pair(sa, ap, am, i, hasA), _yee_pair(sb, bp, bm, i, hasB))

# `fallback === nothing` is the vectorizable native pass; it returns whether
# every result is finite. Otherwise each nonfinite native result is replaced
# by `fallback(I)`, the same discrete expression evaluated exactly. `ivdep` is
# valid: callers reject aliasing between output and sources, and an in-place
# stage reads and writes only the same entry in each iteration. The only
# reduction is the integer finiteness flag, so no float is reassociated.
function _yee_segment!(
    op::F,
    fallback::G,
    out,
    old,
    sa,
    sb,
    base,
    ap,
    am,
    bp,
    bm,
    range,
    hasA::Val,
    hasB::Val,
) where {F,G}
    ok = true
    if fallback === nothing
        @inbounds @simd ivdep for i in range
            value = _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, hasA, hasB)
            out[base+i] = value
            ok &= isfinite(value)
        end
    else
        for i in range
            value = _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, hasA, hasB)
            isfinite(value) || (value = fallback(CartesianIndices(out)[base+i]))
            out[base+i] = value
        end
    end
    return ok
end

@inline function _yee_segment!(
    op,
    fallback,
    out,
    old,
    sa,
    sb,
    base,
    ap,
    am,
    hasA::Bool,
    bp,
    bm,
    hasB::Bool,
    range,
)
    if hasA
        hasB ?
        _yee_segment!(
            op,
            fallback,
            out,
            old,
            sa,
            sb,
            base,
            ap,
            am,
            bp,
            bm,
            range,
            Val(true),
            Val(true),
        ) :
        _yee_segment!(
            op,
            fallback,
            out,
            old,
            sa,
            sb,
            base,
            ap,
            am,
            bp,
            bm,
            range,
            Val(true),
            Val(false),
        )
    else
        hasB ?
        _yee_segment!(
            op,
            fallback,
            out,
            old,
            sa,
            sb,
            base,
            ap,
            am,
            bp,
            bm,
            range,
            Val(false),
            Val(true),
        ) :
        _yee_segment!(
            op,
            fallback,
            out,
            old,
            sa,
            sb,
            base,
            ap,
            am,
            bp,
            bm,
            range,
            Val(false),
            Val(false),
        )
    end
end

# One row-end entry (periodic wrap or absent wall difference), without loop setup.
@inline function _yee_entry!(
    op,
    fallback,
    out,
    old,
    sa,
    sb,
    base,
    ap,
    am,
    hasA::Bool,
    bp,
    bm,
    hasB::Bool,
    i,
)
    value = if hasA
        hasB ? _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, Val(true), Val(true)) :
        _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, Val(true), Val(false))
    else
        hasB ? _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, Val(false), Val(true)) :
        _yee_value(op, old, sa, sb, base, ap, am, bp, bm, i, Val(false), Val(false))
    end
    if fallback !== nothing && !isfinite(value)
        value = fallback(CartesianIndices(out)[base+i])
    end
    @inbounds out[base+i] = value
    return isfinite(value)
end

function _yee_rows!(
    op::F,
    fallback::G,
    out::AbstractArray{T,D},
    old,
    sa,
    sb,
    g::YeeGrid{D},
    ja::Val{JA},
    jb::Val{JB},
    direction::Val{forward},
) where {F,G,T,D,JA,JB,forward}
    n1 = size(out, 1)
    ok = true
    for J in CartesianIndices(Base.tail(size(out)))
        base = _yee_row(out, J)
        ap, am, hasA = _yee_term_offsets(sa, J, ja, g, direction)
        bp, bm, hasB = _yee_term_offsets(sb, J, jb, g, direction)
        first, last = 1, n1
        if JA == 1 || JB == 1
            # Only one difference runs along axis 1; at most two row-end entries.
            if forward
                if g.periodic[1]
                    last = n1 - 1
                    ok &= if JA == 1
                        _yee_entry!(
                            op,
                            fallback,
                            out,
                            old,
                            sa,
                            sb,
                            base,
                            ap - n1,
                            am,
                            hasA,
                            bp,
                            bm,
                            hasB,
                            n1,
                        )
                    else
                        _yee_entry!(
                            op,
                            fallback,
                            out,
                            old,
                            sa,
                            sb,
                            base,
                            ap,
                            am,
                            hasA,
                            bp - n1,
                            bm,
                            hasB,
                            n1,
                        )
                    end
                end
            elseif g.periodic[1]
                first = 2
                ok &= if JA == 1
                    _yee_entry!(
                        op,
                        fallback,
                        out,
                        old,
                        sa,
                        sb,
                        base,
                        ap,
                        am + n1,
                        hasA,
                        bp,
                        bm,
                        hasB,
                        1,
                    )
                else
                    _yee_entry!(
                        op,
                        fallback,
                        out,
                        old,
                        sa,
                        sb,
                        base,
                        ap,
                        am,
                        hasA,
                        bp,
                        bm + n1,
                        hasB,
                        1,
                    )
                end
            else
                first, last = 2, n1 - 1
                for i in (1, n1)
                    ok &= if JA == 1
                        _yee_entry!(op, fallback, out, old, sa, sb, base, 0, 0, false, bp, bm, hasB, i)
                    else
                        _yee_entry!(op, fallback, out, old, sa, sb, base, ap, am, hasA, 0, 0, false, i)
                    end
                end
            end
        end
        ok &= _yee_segment!(
            op,
            fallback,
            out,
            old,
            sa,
            sb,
            base,
            ap,
            am,
            hasA,
            bp,
            bm,
            hasB,
            first:last,
        )
    end
    return ok
end

@inline _yee_curl_combination(::Nothing, ::Nothing) = nothing
@inline _yee_curl_combination(a, ::Nothing) = a
@inline _yee_curl_combination(::Nothing, b) = -b
@inline _yee_curl_combination(a, b) = a - b

struct _YeeCurlOp{T}
    dxa::T
    dxb::T
end

@inline _yee_quotient(::Nothing, dx) = nothing
@inline _yee_quotient(pair::Tuple, dx) = (pair[1] - pair[2]) / dx

@inline function (op::_YeeCurlOp{T})(::Nothing, a, b) where {T}
    value = _yee_curl_combination(_yee_quotient(a, op.dxa), _yee_quotient(b, op.dxb))
    return value === nothing ? zero(T) : value
end

@inline function _yee_curl_component!(
    out,
    fields,
    g::YeeGrid{D,T},
    direction::Val{forward},
    ::Val{c},
    fallback,
) where {D,T,forward,c}
    ja, jb = mod1(c + 1, 3), mod1(c + 2, 3)
    op = _YeeCurlOp{T}(ja <= D ? g.dx[ja] : one(T), jb <= D ? g.dx[jb] : one(T))
    return _yee_rows!(
        op,
        fallback,
        out[c],
        nothing,
        fields[jb],
        fields[ja],
        g,
        Val(ja),
        Val(jb),
        direction,
    )
end

function _yee_curl!(out, fields, g::YeeGrid{D,T}, direction::Val{forward}) where {D,T,forward}
    ok = _yee_curl_component!(out, fields, g, direction, Val(1), nothing)
    ok &= _yee_curl_component!(out, fields, g, direction, Val(2), nothing)
    ok &= _yee_curl_component!(out, fields, g, direction, Val(3), nothing)
    ok && return out
    # Intermediate overflow: evaluate the nonfinite entries exactly and round
    # once; an exact value outside the target range is a signed infinity.
    for (c, component) in ((1, Val(1)), (2, Val(2)), (3, Val(3)))
        exact = I -> _from_exact(T, _yee_curl_exact(fields, I, c, g, forward))
        _yee_curl_component!(out, fields, g, direction, component, exact)
    end
    return out
end

"""
    yee_curl!(out, fields, g; field=:electric)

Compatible Yee curl: electric edges to magnetic faces or magnetic faces to
electric edges, `(f[i+1]-f[i])/dx` differences evaluated in the element type.
Validates shapes, finite input and input/output aliasing before writing. An
entry whose native evaluation overflows is evaluated exactly and rounded once,
so a nonfinite output means the exact curl is outside the floating-point range.
Magnetic input must satisfy homogeneous PEC normal-field constraints.
"""
function yee_curl!(out, fields, g::YeeGrid; field::Symbol = :electric)
    electric = _yee_electric(field)
    _yee_check_shape(fields, g, electric)
    _yee_check_shape(out, g, !electric)
    for a in out, b in fields
        Base.mightalias(a, b) && throw(ArgumentError("Yee curl input and output must not alias"))
    end
    for c = 1:3, k = 1:c-1
        Base.mightalias(out[c], out[k]) && throw(ArgumentError("Yee curl outputs must not alias"))
    end
    for a in fields
        _yee_all_finite(a) || throw(ArgumentError("Yee curl input must be finite"))
    end
    return _yee_curl!(out, fields, g, Val(electric))
end

"""
    yee_divergence!(out, fields, g; field=:magnetic)

Magnetic divergence is cell centered with shape `g.n`. Electric divergence is
on interior nodes: `n[d]-1` on PEC axes and `n[d]` on periodic axes. Boundary
surface charge is outside this vacuum bulk diagnostic. These divergences
annihilate their compatible curls to floating-point roundoff. Input must be
finite; an entry whose native evaluation overflows is evaluated exactly and
rounded once.
"""
function yee_divergence!(out, fields, g::YeeGrid{D,T}; field::Symbol = :magnetic) where {D,T}
    electric = _yee_electric(field)
    _yee_check_shape(fields, g, electric)
    out isa StridedArray{T,D} ||
        throw(ArgumentError("Yee divergence output must match the grid precision and CPU storage"))
    sizes = ntuple(d -> g.n[d] - Int(electric && !g.periodic[d]), D)
    size(out) == sizes || throw(DimensionMismatch("Yee divergence output has the wrong shape"))
    Base.require_one_based_indexing(out)
    any(a -> Base.mightalias(out, a), fields) &&
        throw(ArgumentError("Yee divergence input and output must not alias"))
    for a in fields
        _yee_all_finite(a) || throw(ArgumentError("Yee divergence input must be finite"))
    end
    ok = true
    node(I) = CartesianIndex(ntuple(d -> I[d] + Int(electric && !g.periodic[d]), D))
    for I in CartesianIndices(out)
        J = node(I)
        value = zero(T)
        for d = 1:D
            value += _yee_difference(fields[d], J, d, g, !electric)
        end
        out[I] = value
        ok &= isfinite(value)
    end
    ok && return out
    for I in CartesianIndices(out)
        isfinite(out[I]) && continue
        J = node(I)
        out[I] = _from_exact(T, sum(d -> _yee_difference_exact(fields[d], J, d, g, !electric), 1:D))
    end
    return out
end

"""
    yee_electric_charge!(out, E, g::YeeGrid{1})

Write extensive electric charge on every 1D node, including PEC surface charge,
in epsilon0=1 units with unit transverse area. For PEC the result is
`[Ex[1], diff(Ex)..., -Ex[end]]` (length `n+1`); periodic nodes have
`Ex[i]-Ex[prev(i)]` (length `n`). This is an incidence, without division by dx.
The sum telescopes to zero in exact arithmetic, including the two PEC surfaces.

Inputs must have valid electric shapes, finite data and homogeneous tangential
PEC values; output must be a matching CPU strided vector and must not alias any
input. Checks precede writes. Each entry is one correctly rounded difference or
negation, with signed infinity if it is outside the target range. The existing
`yee_divergence!` retains its bulk-only meaning. For a prescribed current step
the exact charge update is the negative of this same incidence applied to the
current impulse.
"""
function yee_electric_charge!(out, E, g::YeeGrid{1,T}) where {T}
    _yee_validate_electric_data(E, g)
    out isa StridedVector{T} ||
        throw(ArgumentError("Yee charge output must match the grid precision and CPU storage"))
    length(out) == g.n[1] + Int(!g.periodic[1]) ||
        throw(DimensionMismatch("Yee charge output has the wrong shape"))
    Base.require_one_based_indexing(out)
    any(a -> Base.mightalias(out, a), E) &&
        throw(ArgumentError("Yee charge input and output must not alias"))
    x = E[1]
    if g.periodic[1]
        for i in eachindex(out)
            out[i] = x[i] - x[i == 1 ? length(x) : i - 1]
        end
    else
        out[1] = x[1]
        for i = 2:length(x)
            out[i] = x[i] - x[i-1]
        end
        out[end] = -x[end]
    end
    return out
end
