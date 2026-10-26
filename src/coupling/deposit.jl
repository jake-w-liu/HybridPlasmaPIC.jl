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

# Normwise value-range gate for the checked moment path: finite inputs must
# never yield NaN/Inf, and a nonzero factor product that has already underflowed
# far below the grid scale is replayed exactly rather than accumulated.
@inline _deposition_value_lost(val::Union{Float16,Float32,Float64}) =
    !isfinite(val) ||
    (!iszero(val) && abs(val) < floatmin(typeof(val)) * exp2(precision(typeof(val)) + 2))
@inline _deposition_value_lost(val) = !isfinite(val)

# Periodic cell index for a stencil point. `raw = base + o - 1` lies in
# [-1, n+1] for NGP/CIC/TSC stencils on positions s ∈ [0, n], so one conditional
# correction suffices on ordinary grids; mod remains the fallback for tiny grids
# where a single add/subtract can still leave the index out of range.
@inline function _wrapped_cell(raw::Int, n::Int)
    raw < 0 && (raw += n)
    raw >= n && (raw -= n)
    (raw < 0 || raw >= n) && (raw = mod(raw, n))
    return raw + 1
end

# Fused deposition: every output in `outs` accumulates `Σ_p val_p · S_g(x_p)`
# sharing one stencil computation per particle. `factorss[k]` is the factor
# tuple of output k. Per output the arithmetic order matches the unfused
# kernel exactly — val = first factor; val *= each further factor in order;
# w = val; w *= S_d for d = 1..D; out[idx] += w — so each output is bitwise
# identical to a separate pass. Eligibility (exact replay admission) is tracked
# per output as before; the range gate is normwise: a nonfinite factor product
# or one already far below working precision is flagged per particle and output,
# and each accumulated cell is checked finite once at the end.
# Factor products are formed left to right exactly as in the scalar kernel.
@inline _moment_factor_value(values, p) = @inbounds values[firstindex(values)+p-1]
@inline _moment_factor_product(factors::Tuple, p) =
    _moment_factor_product(Base.tail(factors), p, _moment_factor_value(first(factors), p))
@inline _moment_factor_product(::Tuple{}, p, val) = val
@inline _moment_factor_product(factors::Tuple, p, val) =
    _moment_factor_product(Base.tail(factors), p, val * _moment_factor_value(first(factors), p))
@inline _moment_factors_exact(factors::Tuple, p) =
    all(values -> _deposition_exact_operand(_moment_factor_value(values, p)), factors)
@inline _moment_factors_nonzero(factors::Tuple, p) =
    all(values -> !iszero(_moment_factor_value(values, p)), factors)

function _deposit_moments_kernel!(
    outs::NTuple{K,Array{T,D}},
    ps::ParticleSet{D,T},
    factorss::Tuple,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    ::Val{checked},
    particles = 1:nparticles(ps),
) where {D,T,K,checked}
    length(factorss) == K || throw(DimensionMismatch("one factor tuple is required per output"))
    for out in outs
        fill!(out, zero(T))
    end
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    emask = (one(UInt64) << K) - one(UInt64)
    lmask = zero(UInt64)
    @inbounds for p in particles
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        # Tuple maps unroll over heterogeneous factor tuples; values share type T.
        vals = map(factors -> _moment_factor_product(factors, p), factorss)::NTuple{K,T}
        if checked
            exact = map(factors -> _moment_factors_exact(factors, p), factorss)
            nonzero = map(factors -> _moment_factors_nonzero(factors, p), factorss)
            # Normwise underflow gate: `smin` is the smallest nonzero stencil
            # weight product this particle deposits (factors >= 0, partition of
            # unity => smin > 0). If |val| * smin falls below floatmin, the
            # smallest contribution underflows to zero or subnormal; a field at
            # that scale loses it silently, so replay is required.
            smin = one(T)
            for d = 1:D
                smin *= minimum(w -> iszero(w) ? one(T) : w, st[d][2])
            end
            for k = 1:K
                emask &= ~(UInt64(!exact[k]) << (k - 1))
                lost_k =
                    _deposition_value_lost(vals[k]) ||
                    (iszero(vals[k]) && nonzero[k]) ||
                    (!iszero(vals[k]) && abs(vals[k]) * smin < floatmin(T))
                lmask |= UInt64(lost_k) << (k - 1)
            end
        end
        for c in stamp
            o = Tuple(c)
            idx = ntuple(d -> _wrapped_cell(st[d][1] + o[d] - 1, n[d]), D)
            sfactors = ntuple(d -> st[d][2][o[d]], D)
            checked && (all(_deposition_exact_operand, sfactors) || (emask = zero(UInt64)))
            for k = 1:K
                w = vals[k]
                for d = 1:D
                    w = w * sfactors[d]
                end
                outs[k][idx...] += w
            end
        end
    end
    if checked
        for k = 1:K
            out = outs[k]
            for I in eachindex(out)
                lmask |= UInt64(!isfinite(out[I])) << (k - 1)
            end
        end
    end
    return emask, lmask
end

# Per-output replay flags decoded from the kernel's packed masks. Literal
# indices keep heterogeneous factor/state tuples fully inferable and
# allocation-free.
@inline _moment_state(masks::Tuple{UInt64,UInt64}, k::Int) =
    (; eligible = Bool((masks[1] >> (k - 1)) & 1), lost = Bool((masks[2] >> (k - 1)) & 1))

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
        allnz = !iszero(val)
        for values in Base.tail(factors)
            factor = values[firstindex(values)+p-1]
            if checked
                eligible &= _deposition_exact_operand(factor)
                allnz &= !iszero(factor)
            end
            val = val * factor
        end
        if checked
            # Same normwise underflow gate as the fused kernel: replay when the
            # smallest nonzero stencil product would underflow |val| * smin
            # below floatmin, or the factor product flushed nonzero to zero.
            smin = one(T)
            for d = 1:D
                smin *= minimum(w -> iszero(w) ? one(T) : w, st[d][2])
            end
            lost |=
                _deposition_value_lost(val) ||
                (iszero(val) && allnz) ||
                (!iszero(val) && abs(val) * smin < floatmin(T))
        end
        for c in stamp
            o = Tuple(c)
            w = val
            for d = 1:D
                sfactor = st[d][2][o[d]]
                checked && (eligible &= _deposition_exact_operand(sfactor))
                w = w * sfactor
            end
            idx = ntuple(d -> _wrapped_cell(st[d][1] + o[d] - 1, n[d]), D)
            out[idx...] += w
        end
    end
    if checked
        for I in eachindex(out)
            lost |= !isfinite(out[I])
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
