# Moments.jl — density/momentum/current/pressure moments (from deposit.jl)

@inline function _deposit_weighted_scalar!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    weight::AbstractVector{T},
    value::AbstractVector{T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    np = nparticles(ps)
    length(weight) == np ||
        throw(DimensionMismatch("weight length $(length(weight)) must equal particle count $np"))
    length(value) == np ||
        throw(DimensionMismatch("value length $(length(value)) must equal particle count $np"))
    size(out) == g.n ||
        throw(DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"))
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, wi) in enumerate(eachindex(weight))
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        val = weight[wi] * value[wi]
        for c in stamp
            o = Tuple(c)
            w = val
            for d = 1:D
                w *= st[d][2][o[d]]
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            out[idx...] += w
        end
    end
    return out
end

@inline function _deposit_weighted_product!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    weight::AbstractVector{T},
    left::AbstractVector{T},
    right::AbstractVector{T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    np = nparticles(ps)
    length(weight) == np ||
        throw(DimensionMismatch("weight length $(length(weight)) must equal particle count $np"))
    length(left) == np ||
        throw(DimensionMismatch("left length $(length(left)) must equal particle count $np"))
    length(right) == np ||
        throw(DimensionMismatch("right length $(length(right)) must equal particle count $np"))
    size(out) == g.n ||
        throw(DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"))
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, wi) in enumerate(eachindex(weight))
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        val = weight[wi] * left[wi] * right[wi]
        for c in stamp
            o = Tuple(c)
            w = val
            for d = 1:D
                w *= st[d][2][o[d]]
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            out[idx...] += w
        end
    end
    return out
end

# Accumulate one weighted covariance component with the parallel/merge form of
# Welford's recurrence. Using W*a/(W+a) * Δx*Δy avoids raw vᵢvⱼ moments and
# their catastrophic Inf-Inf cancellation for a hot numerical offset with a
# representable thermal spread. The exponent-scaled coefficient and fallback
# product also preserve contributions whose direct intermediates would
# underflow or overflow.
@inline function _deposit_weighted_covariance!(
    out::Array{T,D},
    mean_left::Array{T,D},
    mean_right::Array{T,D},
    weight_sum::Array{T,D},
    ps::ParticleSet{D,T},
    left::AbstractVector{T},
    right::AbstractVector{T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    same_component::Bool,
) where {D,T}
    fill!(out, zero(T))
    fill!(mean_left, zero(T))
    same_component || fill!(mean_right, zero(T))
    fill!(weight_sum, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for p in eachindex(ps.weight)
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        for c in stamp
            o = Tuple(c)
            cell_weight = ps.weight[p]
            for d = 1:D
                cell_weight *= st[d][2][o[d]]
            end
            iszero(cell_weight) && continue
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            old_weight = weight_sum[idx...]
            if iszero(old_weight)
                mean_left[idx...] = left[p]
                same_component || (mean_right[idx...] = right[p])
                weight_sum[idx...] = cell_weight
                continue
            end
            new_weight = old_weight + cell_weight
            (isfinite(new_weight) && new_weight > zero(T)) || throw(
                ArgumentError("pressure_tensor!: deposited particle weight exceeds numeric range"),
            )
            left_delta = left[p] - mean_left[idx...]
            right_delta = same_component ? left_delta : right[p] - mean_right[idx...]
            coefficient = _finite_muldiv(old_weight, cell_weight, new_weight)
            direct_contribution = zero(T)
            use_direct = isfinite(left_delta) && isfinite(right_delta) && coefficient > zero(T)
            if use_direct
                # Attach the small coefficient to the larger deviation first:
                # this retains the exact ordinary-range recurrence while
                # reducing both underflow and overflow risk.
                if abs(left_delta) >= abs(right_delta)
                    first_term = coefficient * left_delta
                    direct_contribution = first_term * right_delta
                    use_direct &= !iszero(first_term) || iszero(left_delta)
                else
                    first_term = coefficient * right_delta
                    direct_contribution = first_term * left_delta
                    use_direct &= !iszero(first_term) || iszero(right_delta)
                end
                use_direct &= isfinite(direct_contribution)
            end
            if use_direct
                out[idx...] += direct_contribution
            else
                out[idx...] += _weighted_covariance_contribution(
                    old_weight,
                    cell_weight,
                    new_weight,
                    left[p],
                    mean_left[idx...],
                    right[p],
                    same_component ? mean_left[idx...] : mean_right[idx...],
                )
            end
            mean_left[idx...] = _convex_weighted_mean(
                mean_left[idx...],
                left[p],
                old_weight,
                cell_weight,
                new_weight,
            )
            if !same_component
                mean_right[idx...] = _convex_weighted_mean(
                    mean_right[idx...],
                    right[p],
                    old_weight,
                    cell_weight,
                    new_weight,
                )
            end
            weight_sum[idx...] = new_weight
        end
    end
    return out
end

@inline function _centered_frexp(value::T, center::T) where {T}
    delta = value - center
    if isfinite(delta)
        return frexp(delta)
    end
    # Opposite near-range values can have an infinite direct difference even
    # though a sufficiently small weight makes the covariance representable.
    scale = max(abs(value), abs(center))
    isfinite(scale) || return T(NaN), 0
    iszero(scale) && return zero(T), 0
    normalized_delta = value / scale - center / scale
    df, de = frexp(normalized_delta)
    sf, se = frexp(scale)
    mf, me = frexp(df * sf)
    return mf, de + se + me
end

@inline function _weighted_covariance_contribution(
    old_weight::T,
    cell_weight::T,
    new_weight::T,
    left::T,
    left_center::T,
    right::T,
    right_center::T,
) where {T}
    lf, le = _centered_frexp(left, left_center)
    rf, re = _centered_frexp(right, right_center)
    (iszero(lf) || iszero(rf)) && return zero(T)
    of, oe = frexp(old_weight)
    wf, we = frexp(cell_weight)
    nf, ne = frexp(new_weight)
    return ldexp((of * wf * lf * rf) / nf, oe + we + le + re - ne)
end

@inline function _pressure_floor_adjustment(
    mass::T,
    left_momentum::T,
    right_momentum::T,
    density::T,
    density_floor::T,
) where {T}
    gap = density_floor - density
    (iszero(mass) || iszero(left_momentum) || iszero(right_momentum) || iszero(gap)) &&
        return zero(T)
    mf, me = frexp(mass)
    lf, le = frexp(left_momentum)
    rf, re = frexp(right_momentum)
    gf, ge = frexp(gap)
    df, de = frexp(density)
    ff, fe = frexp(density_floor)
    value, adjustment = frexp((mf * lf * rf * gf) / (df * ff))
    return ldexp(value, me + le + re + ge - de - fe + adjustment)
end

function _pressure_requires_extended_precision(ps::ParticleSet)
    accumulation_limit = floatmax(eltype(ps.weight)) / max(nparticles(ps), 1)
    @inbounds for p in eachindex(ps.weight)
        w = ps.weight[p]
        for c = 1:3
            v = ps.v[c][p]
            moment = w * v
            (!isfinite(moment) || (iszero(moment) && !iszero(w) && !iszero(v))) && return true
        end
        for (i, j) in _PT_PAIRS
            vi = ps.v[i][p]
            vj = ps.v[j][p]
            second = w * vi * vj
            (
                !isfinite(second) ||
                (iszero(second) && !iszero(w) && !iszero(vi) && !iszero(vj)) ||
                abs(second) > accumulation_limit
            ) && return true
        end
    end
    return false
end

@inline function _pressure_extended_precision(::Type{T}) where {T<:AbstractFloat}
    exponent_span = exponent(floatmax(T)) - exponent(nextfloat(zero(T)))
    # A deposited covariance term contains up to three finite T factors. Four
    # exponent spans plus significand/headroom bits make their sums and the
    # subsequent mean correction exact for every possible particle count.
    return 4exponent_span + 4precision(T) + 128
end

@inline _pressure_extended_precision(::Type{BigFloat}) = max(precision(BigFloat), 8728)

function _deposit_pressure_extended!(P, nbuf, mom, ps, g, shape, density_floor)
    T = eltype(ps.weight)
    precision_bits = _pressure_extended_precision(T)
    setprecision(BigFloat, precision_bits) do
        weight_sum = zeros(BigFloat, g.n)
        momentum_sum = ntuple(_ -> zeros(BigFloat, g.n), 3)
        second_sum = ntuple(_ -> zeros(BigFloat, g.n), 6)
        n = g.n
        D = ndims(weight_sum)
        stamp = CartesianIndices(ntuple(_ -> width(shape), D))
        @inbounds for p in eachindex(ps.weight)
            st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
            velocity = ntuple(component -> BigFloat(ps.v[component][p]), 3)
            for c in stamp
                o = Tuple(c)
                cell_weight = ps.weight[p]
                for d = 1:D
                    cell_weight *= st[d][2][o[d]]
                end
                iszero(cell_weight) && continue
                idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
                weight = BigFloat(cell_weight)
                weight_sum[idx...] += weight
                for component = 1:3
                    momentum_sum[component][idx...] += weight * velocity[component]
                end
                for (component, (i, j)) in enumerate(_PT_PAIRS)
                    second_sum[component][idx...] += weight * velocity[i] * velocity[j]
                end
            end
        end

        cell_volume = prod(BigFloat(dx) for dx in g.dx)
        mass = BigFloat(ps.m)
        floor_big = BigFloat(density_floor)
        for I in eachindex(weight_sum)
            density_big = weight_sum[I] / cell_volume
            nbuf[I] = T(density_big)
            for component = 1:3
                mom[component][I] = T(momentum_sum[component][I] / cell_volume)
            end
            denominator = max(density_big, floor_big)
            for (component, (i, j)) in enumerate(_PT_PAIRS)
                left_momentum = momentum_sum[i][I] / cell_volume
                right_momentum = momentum_sum[j][I] / cell_volume
                pressure =
                    mass * (
                        second_sum[component][I] / cell_volume -
                        left_momentum * right_momentum / denominator
                    )
                P[component][I] = T(pressure)
            end
        end
    end
    return P
end

function _validate_pressure_scratch_aliasing(P, ps, work, nbuf, mom)
    grids = (P..., nbuf, mom...)
    for i in eachindex(grids)
        _mightalias_particle_storage(grids[i], ps) &&
            throw(ArgumentError("pressure_tensor!: grid buffers must not alias particle storage"))
        for j = (i+1):length(grids)
            Base.mightalias(grids[i], grids[j]) && throw(
                ArgumentError("pressure_tensor!: grid output and scratch buffers must not alias"),
            )
        end
    end
    if work !== nothing
        _mightalias_particle_storage(work, ps) &&
            throw(ArgumentError("pressure_tensor!: work must not alias particle storage"))
        for grid in grids
            Base.mightalias(work, grid) && throw(
                ArgumentError(
                    "pressure_tensor!: work must not alias grid output or scratch buffers",
                ),
            )
        end
    end
    return nothing
end

function density!(
    nout::Array{T,D},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    deposit_scalar!(nout, ps, ps.weight, g, shape)
    nout ./= prod(g.dx)
    return nout
end

"Momentum density (n u)_c = (Σ_p w_p v_{c,p} S_g(x_p)) / ΔV for c = 1,2,3."
function momentum!(
    mom::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::Union{Nothing,AbstractVector{T}} = nothing,
) where {D,T}
    ΔV = prod(g.dx)
    if work === nothing
        for c = 1:3
            _deposit_weighted_scalar!(mom[c], ps, ps.weight, ps.v[c], g, shape)
            mom[c] ./= ΔV
        end
    else
        length(work) == nparticles(ps) || throw(
            DimensionMismatch(
                "work length $(length(work)) must equal particle count $(nparticles(ps))",
            ),
        )
        for c = 1:3
            @. work = ps.weight * ps.v[c]
            deposit_scalar!(mom[c], ps, work, g, shape)
            mom[c] ./= ΔV
        end
    end
    return mom
end

"Ion current density J_c = q · (Σ_p w_p v_{c,p} S_g(x_p))/ΔV."
function current!(
    J::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::Union{Nothing,AbstractVector{T}} = nothing,
) where {D,T}
    momentum!(J, ps, g, shape; work)
    for c = 1:3
        J[c] .*= ps.q
    end
    return J
end

# (i,j) index pairs of the 6 independent symmetric-tensor components, in order
# (xx, yy, zz, xy, xz, yz).
const _PT_PAIRS = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))

"""
    pressure_tensor!(P, ps, g, shape; nfloor=1e-6, work=..., nbuf=..., mom=...)

Deposit the ion pressure tensor `P_ij = m·Σ_p w_p (v_i−U_i)(v_j−U_j) S / ΔV`
(`U` the local bulk velocity), as 6 components in the order (xx,yy,zz,xy,xz,yz).
Computed as the centered second moment: second-moment density − ρ U_i U_j.
"""
function pressure_tensor!(
    P::NTuple{6,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    nfloor = 1e-6,
    work::Union{Nothing,AbstractVector{T}} = nothing,
    nbuf::Union{Nothing,Array{T,D}} = nothing,
    mom::Union{Nothing,NTuple{3,<:Array{T,D}}} = nothing,
) where {D,T}
    nf = _require_finite_positive_real("nfloor", nfloor, T)
    Np = nparticles(ps)
    nbuf === nothing && (nbuf = similar(P[1]))
    mom === nothing && (mom = ntuple(_ -> similar(P[1]), 3))
    for c = 1:6
        size(P[c]) == g.n ||
            throw(DimensionMismatch("P[$c] size $(size(P[c])) does not match grid size $(g.n)"))
    end

    size(nbuf) == g.n ||
        throw(DimensionMismatch("nbuf size $(size(nbuf)) does not match grid size $(g.n)"))
    for c = 1:3
        size(mom[c]) == g.n ||
            throw(DimensionMismatch("mom[$c] size $(size(mom[c])) does not match grid size $(g.n)"))
    end

    work === nothing ||
        length(work) == Np ||
        throw(DimensionMismatch("work length $(length(work)) must equal particle count $Np"))
    _validate_pressure_scratch_aliasing(P, ps, work, nbuf, mom)

    ΔV = prod(g.dx)
    mq = ps.m
    stable_weights = all(w -> isfinite(w) && w >= zero(T), ps.weight)
    if stable_weights && _pressure_requires_extended_precision(ps)
        _deposit_pressure_extended!(P, nbuf, mom, ps, g, shape, nf)
    elseif stable_weights
        # The six independent online accumulations reuse `mom[1:2]` and
        # `nbuf` as scratch. Recompute the documented density/momentum outputs
        # once afterward; total stencil passes remain the same as the former
        # raw-second-moment implementation.
        for (idx, (i, j)) in enumerate(_PT_PAIRS)
            _deposit_weighted_covariance!(
                P[idx],
                mom[1],
                mom[2],
                nbuf,
                ps,
                ps.v[i],
                ps.v[j],
                g,
                shape,
                i == j,
            )
        end
        density!(nbuf, ps, g, shape)
        momentum!(mom, ps, g, shape; work)
        for (idx, (i, j)) in enumerate(_PT_PAIRS)
            Pij = P[idx]
            mi = mom[i]
            mj = mom[j]
            @inbounds for I in eachindex(Pij)
                pressure = _finite_muldiv(mq, Pij[I], ΔV)
                nv = nbuf[I]
                if zero(T) < nv < nf
                    # Preserve the existing density-floor definition:
                    # Sᵢⱼ - mᵢmⱼ/nfloor = covariance +
                    # mᵢmⱼ(1/n - 1/nfloor).
                    pressure += _pressure_floor_adjustment(mq, mi[I], mj[I], nv, nf)
                end
                Pij[I] = pressure
            end
        end
    else
        # Preserve legacy arithmetic for non-physical signed/non-finite weights;
        # the stable recurrence requires a non-negative measure.
        density!(nbuf, ps, g, shape)
        momentum!(mom, ps, g, shape; work)
        for (idx, (i, j)) in enumerate(_PT_PAIRS)
            vi = ps.v[i]
            vj = ps.v[j]
            if work === nothing
                _deposit_weighted_product!(P[idx], ps, ps.weight, vi, vj, g, shape)
            else
                @inbounds @. work = ps.weight * vi * vj
                deposit_scalar!(P[idx], ps, work, g, shape)
            end
            Pij = P[idx]
            mi = mom[i]
            mj = mom[j]
            @inbounds for I in eachindex(Pij)
                nv = max(nbuf[I], nf)
                Pij[I] = mq * (Pij[I] / ΔV - mi[I] * mj[I] / nv)
            end
        end
    end
    return P
end

"""
    temperature_components(P, n; nfloor=1e-6)

Per-cell directional temperatures `(T_x, T_y, T_z) = (P_xx, P_yy, P_zz)/n` from a
pressure tensor produced by [`pressure_tensor!`](@ref).
"""
function temperature_components(P::NTuple{6,<:Array{T,D}}, n::Array{T,D}; nfloor = 1e-6) where {T,D}
    nf = T(nfloor)
    isfinite(nf) && nf > zero(T) || throw(ArgumentError("nfloor must be finite and positive"))
    return ntuple(c -> P[c] ./ max.(n, nf), 3)
end
