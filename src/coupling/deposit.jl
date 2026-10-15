# Deposit.jl — scalar deposition (from deposit.jl)

function _validate_deposition_inputs(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    vals::AbstractVector,
    g::FourierGrid{D,T},
) where {D,T}
    length(vals) == nparticles(ps) || throw(DimensionMismatch("vals length ≠ particle count"))
    size(out) == g.n ||
        throw(DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"))
    (_mightalias_particle_storage(out, ps) || _particle_arrays_mightalias(out, vals)) &&
        throw(ArgumentError("deposition output must not alias particle storage or values"))
    _validate_particle_positions(ps)
    return nothing
end

# Symbolic Real coefficients (for example π) need not have a rational
# representation. They retain legacy arithmetic instead of entering the exact
# replay; supported represented operands must also be finite.
@inline _deposition_exact_operand(value) = false
@inline _deposition_exact_operand(value::Union{AbstractFloat,Integer,Rational}) = isfinite(value)

@inline _deposition_range_lost(x::Union{Float16,Float32,Float64}) = !isfinite(x) || issubnormal(x)
@inline _deposition_range_lost(x) = !isfinite(x)
@inline _deposition_product_lost(a, b, value) =
    _deposition_range_lost(value) || (iszero(value) && !iszero(a) && !iszero(b))

# Generic deposition and normalized moments share this arithmetic. Checking is
# specialized away for generic deposition, whose rounded raw-sum contract stays
# unchanged. Moment callers retain the original factors for an exact cold replay.
function _deposit_scalar_kernel!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    factors::Tuple,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    ::Val{checked},
    particles = 1:nparticles(ps),
) where {D,T,checked}
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    eligible = true
    lost = false
    @inbounds for p in particles
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        firstvals = first(factors)
        val = firstvals[firstindex(firstvals)+p-1]
        checked && (eligible &= _deposition_exact_operand(val))
        for values in Base.tail(factors)
            factor = values[firstindex(values)+p-1]
            product = val * factor
            if checked
                eligible &= _deposition_exact_operand(factor)
                lost |= _deposition_product_lost(val, factor, product)
            end
            val = product
        end
        for c in stamp
            o = Tuple(c)
            w = val
            for d = 1:D
                factor = st[d][2][o[d]]
                product = w * factor
                if checked
                    eligible &= _deposition_exact_operand(factor)
                    lost |= _deposition_product_lost(w, factor, product)
                end
                w = product
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            value = out[idx...] + w
            out[idx...] = value
            if checked
                stored = out[idx...]
                lost |=
                    _deposition_range_lost(value) ||
                    _deposition_range_lost(stored) ||
                    (iszero(stored) && !iszero(value))
            end
        end
    end
    return (; eligible, lost)
end

"""
    deposit_scalar!(out, ps, vals, g, shape)

Accumulate `out_g = Σ_p vals[p] · S_g(x_p)` for the given shape (periodic wrap).
`out` is zeroed first. Partition of unity ⇒ `Σ_g out_g = Σ_p vals[p]`.
The output must not overlap `vals` or any particle storage.
"""
function deposit_scalar!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    vals::AbstractVector,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    _validate_deposition_inputs(out, ps, vals, g)
    _deposit_scalar_kernel!(out, ps, (vals,), g, shape, Val(false))
    return out
end
