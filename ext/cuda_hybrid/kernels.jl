struct _CUDARNOps <: _CH._ScalarOps end
@inline _CH._scalar_add(::_CUDARNOps, a::Float64, b::Float64) = CUDA.add_rn(a, b)
@inline _CH._scalar_sub(::_CUDARNOps, a::Float64, b::Float64) = CUDA.sub_rn(a, b)
@inline _CH._scalar_mul(::_CUDARNOps, a::Float64, b::Float64) = CUDA.mul_rn(a, b)
@inline _CH._scalar_div(::_CUDARNOps, a::Float64, b::Float64) = CUDA.div_rn(a, b)
for op in (:_scalar_add, :_scalar_sub, :_scalar_mul, :_scalar_div)
    @eval begin
        @inline _CH.$op(o::_CUDARNOps, a::Float64, b::Int) = _CH.$op(o, a, Float64(b))
        @inline _CH.$op(o::_CUDARNOps, a::Int, b::Float64) = _CH.$op(o, Float64(a), b)
    end
end

@inline _ch_index() =
    (Int(CUDA.blockIdx().x) - 1) * Int(CUDA.blockDim().x) + Int(CUDA.threadIdx().x)

@inline function _ch_flag!(status, category::Int, index::Int, component::Int)
    if CUDA.atomic_cas!(pointer(status, 1), Int64(0), Int64(category)) == 0
        @inbounds status[2] = index
        @inbounds status[3] = component
    end
    return nothing
end

function _ch_particle_values!(status, x, v, weight, L, count)
    p = _ch_index()
    if p <= count
        @inbounds begin
            _ch_supported(x[p]) && 0 <= x[p] < L || _ch_flag!(status, 1, p, 1)
            _ch_supported(weight[p]) && weight[p] >= 0 || _ch_flag!(status, 1, p, 5)
            for c = 1:3
                _ch_supported(v[c][p]) || _ch_flag!(status, 1, p, c + 1)
            end
        end
    end
    return nothing
end

function _ch_field_values!(status, arrays, density, pressure, inverse, nfloor, check_density, count)
    i = _ch_index()
    if i <= count
        @inbounds for c in eachindex(arrays)
            _ch_supported(arrays[c][i]) || _ch_flag!(status, 2, i, c)
        end
        @inbounds if check_density
            density[i] >= nfloor || _ch_flag!(status, 2, i, 1)
            pressure[i] >= 0 || _ch_flag!(status, 2, i, 14)
            inverse[i] > 0 || _ch_flag!(status, 2, i, 16)
        end
    end
    return nothing
end

function _ch_vector_values!(status, vector, count)
    i = _ch_index()
    if i <= count
        @inbounds for c in eachindex(vector)
            _ch_supported(vector[c][i]) || _ch_flag!(status, 2, i, c)
        end
    end
    return nothing
end

@inline function _ch_stencil(x, L, dx, n)
    s = _CH._periodic_cell_position(x, L, dx)
    base = floor(s)
    if !(isfinite(s) && isfinite(base) && 0 <= base <= n)
        return ((0, (0.0, 0.0)),), false
    end
    integer_base = unsafe_trunc(Int, base)
    return ((integer_base, _CH._cic_weights(s, integer_base, _CUDARNOps())),), true
end

function _ch_gather!(status, out, field, x, L, dx, n, count)
    p = _ch_index()
    if p <= count
        @inbounds st, valid = _ch_stencil(x[p], L, dx, n)
        if !valid
            _ch_flag!(status, 3, p, 0)
            return nothing
        end
        value = _CH._gather_vector_stencil(
            st,
            (n,),
            CartesianIndices((2,)),
            Float64,
            Val(:periodic),
            Val(false),
            _CUDARNOps(),
        ) do idx
            @inbounds (field[1][idx[1]], field[2][idx[1]], field[3][idx[1]])
        end
        @inbounds for c = 1:3
            _ch_supported(value[c]) || _ch_flag!(status, 3, p, c)
            out[c][p] = value[c]
        end
    end
    return nothing
end

function _ch_scales!(status, out, v, E, B, coefficient, mode, count)
    p = _ch_index()
    if p <= count
        @inbounds begin
            u = (v[1][p], v[2][p], v[3][p])
            e = (E[1][p], E[2][p], E[3][p])
            b = (B[1][p], B[2][p], B[3][p])
            r =
                mode isa Val{:boris} ? _CH._boris_scale_record(u, e, b, coefficient, _CUDARNOps()) :
                _CH._euler_scale_record(u, e, b, coefficient, _CUDARNOps())
            all(isfinite, (r.vmax, r.emax, r.bmax, r.scale)) || _ch_flag!(status, 4, p, 0)
            out[1][p], out[2][p], out[3][p] = r.vmax, r.emax, r.bmax
            out[4][p] = _CH._batch_nonzero_min(Inf, r.scale)
        end
    end
    return nothing
end

function _ch_euler!(status, out, v, E, B, hq, mode, count)
    p = _ch_index()
    if p <= count
        @inbounds begin
            u = (v[1][p], v[2][p], v[3][p])
            e = (E[1][p], E[2][p], E[3][p])
            b = (B[1][p], B[2][p], B[3][p])
            value = first(_CH._euler_kick_parts(u..., e..., b..., hq, mode, _CUDARNOps()))
            for c = 1:3
                _ch_supported(value[c]) || _ch_flag!(status, 4, p, c)
                out[c][p] = value[c]
            end
        end
    end
    return nothing
end

function _ch_boris_drift!(status, xnew, vnew, xmid, x, v, E, B, qm, h, L, count)
    p = _ch_index()
    if p <= count
        @inbounds begin
            value = first(
                _CH._boris_kick_parts(
                    v[1][p],
                    v[2][p],
                    v[3][p],
                    E[1][p],
                    E[2][p],
                    E[3][p],
                    B[1][p],
                    B[2][p],
                    B[3][p],
                    qm,
                    h,
                    Val(false),
                    _CUDARNOps(),
                ),
            )
            for c = 1:3
                _ch_supported(value[c]) || _ch_flag!(status, 4, p, c)
                vnew[c][p] = value[c]
            end
            endpoint = _CH._particle_coordinate_native(x[p], value[1], h, _CUDARNOps())
            midpoint = _CH._particle_coordinate_native(x[p], value[1], h / 2, _CUDARNOps())
            if !(isfinite(endpoint) && isfinite(midpoint))
                _ch_flag!(status, 4, p, 4)
                return nothing
            end
            endpoint = _CH._periodic_coordinate(endpoint, 0.0, L, L)
            midpoint = _CH._periodic_coordinate(midpoint, 0.0, L, L)
            _ch_supported(endpoint) && 0 <= endpoint < L || _ch_flag!(status, 4, p, 4)
            _ch_supported(midpoint) && 0 <= midpoint < L || _ch_flag!(status, 4, p, 5)
            xnew[p], xmid[p] = endpoint, midpoint
        end
    end
    return nothing
end

function _ch_deposit!(status, raw, x, v, weight, L, dx, n, count)
    p = _ch_index()
    if p <= count
        @inbounds st, valid = _ch_stencil(x[p], L, dx, n)
        if !valid
            _ch_flag!(status, 3, p, 0)
            return nothing
        end
        factors = ((weight,), (weight, v[1]), (weight, v[2]), (weight, v[3]))
        values = map(f -> _CH._moment_factor_product_ops(f, p, _CUDARNOps()), factors)
        smin = _CH._moment_stencil_minimum(st, Float64, _CUDARNOps())
        for k = 1:4
            r = _CH._moment_particle_replay(
                values[k],
                _CH._moment_factors_exact(factors[k], p),
                _CH._moment_factors_nonzero(factors[k], p),
                smin,
                _CUDARNOps(),
            )
            if !r.eligible || r.lost
                _ch_flag!(status, 5, p, k)
                return nothing
            end
        end
        @inbounds for o = 1:2
            sfactors = (st[1][2][o],)
            if !all(_CH._deposition_exact_operand, sfactors)
                _ch_flag!(status, 5, p, 0)
                return nothing
            end
            node = _CH._wrapped_cell(st[1][1] + o - 1, n)
            for k = 1:4
                value = _CH._moment_stencil_value(values[k], sfactors, _CUDARNOps())
                CUDA.atomic_add!(pointer(raw[k], node), value)
            end
        end
    end
    return nothing
end

function _ch_accumulate!(status, out, raw, q, dx, count)
    i = _ch_index()
    if i <= count
        @inbounds for k = 1:4
            value, lost = _CH._moment_candidate(raw[k][i], 1.0, dx, Val(:postscale), _CUDARNOps())
            part = _CH._charge_moment_parts(out[k][i], q, value, _CUDARNOps())
            if lost ||
               _CH._deposition_product_lost(q, value, part.increment) ||
               !_ch_supported(part.value)
                _ch_flag!(status, 5, i, k)
            end
            out[k][i] = part.value
        end
    end
    return nothing
end

function _ch_normalize!(status, n, u, nfloor, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            if !(_ch_supported(n[i]) && n[i] >= nfloor)
                _ch_flag!(status, 5, i, 1)
                return nothing
            end
            inv = _CH._hybrid_density_inverse(n[i], Float64, _CUDARNOps())
            _ch_supported(inv) && inv > 0 || _ch_flag!(status, 5, i, 5)
            for c = 1:3
                value = _CH._ion_moment_value(u[c][i], inv, _CUDARNOps())
                _CH._deposition_product_lost(u[c][i], inv, value) && _ch_flag!(status, 5, i, c + 1)
                _ch_supported(value) || _ch_flag!(status, 5, i, c + 1)
                u[c][i] = value
            end
        end
    end
    return nothing
end

function _ch_pressure!(status, pe, ninv, n, Te, nfloor, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            if !(_ch_supported(n[i]) && n[i] >= nfloor)
                _ch_flag!(status, 6, i, 1)
                return nothing
            end
            pressure = _CH._isothermal_pressure_value(Te, n[i], _CUDARNOps())
            inv = _CH._hybrid_density_inverse(n[i], Float64, _CUDARNOps())
            if _CH._deposition_product_lost(Te, n[i], pressure) ||
               !_ch_supported(pressure) ||
               !_ch_supported(inv) ||
               inv <= 0
                _ch_flag!(status, 6, i, 1)
            end
            pe[i], ninv[i] = pressure, inv
        end
    end
    return nothing
end

function _ch_ohm!(status, E, u, B, J, gp, ninv, eta, etaH, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            b = (B[1][i], B[2][i], B[3][i])
            v = (u[1][i], u[2][i], u[3][i])
            j = (J[1][i], J[2][i], J[3][i])
            inv = ninv[i]
            r1 = _CH._ohm_scalar_native(
                v[2],
                v[3],
                b[2],
                b[3],
                j[2],
                j[3],
                inv,
                eta,
                j[1],
                etaH,
                0.0,
                gp[1][i],
                _CUDARNOps(),
            )
            r2 = _CH._ohm_scalar_native(
                v[3],
                v[1],
                b[3],
                b[1],
                j[3],
                j[1],
                inv,
                eta,
                j[2],
                etaH,
                0.0,
                nothing,
                _CUDARNOps(),
            )
            r3 = _CH._ohm_scalar_native(
                v[1],
                v[2],
                b[1],
                b[2],
                j[1],
                j[2],
                inv,
                eta,
                j[3],
                etaH,
                0.0,
                nothing,
                _CUDARNOps(),
            )
            for (c, r) in enumerate((r1, r2, r3))
                r.needs_completion && _ch_flag!(status, 6, i, c)
                _ch_supported(r.value) || _ch_flag!(status, 6, i, c)
                E[c][i] = r.value
            end
        end
    end
    return nothing
end

function _ch_negate!(status, vector, count)
    i = _ch_index()
    if i <= count
        @inbounds for c = 1:3
            vector[c][i] = -vector[c][i] # Exact unary sign of canonical _faraday_negate!.
        end
    end
    return nothing
end

function _ch_affine!(status, out, initial, rates, op, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            out[1][i] = initial[1][i] # D1 Bx rate is identically zero; retain its stored bits.
            for c = 2:3
                r =
                    length(rates) == 1 ?
                    _CH._rk4_affine_native(op, initial[c][i], rates[1][c][i], _CUDARNOps()) :
                    _CH._rk4_affine_native(
                        op,
                        initial[c][i],
                        rates[1][c][i],
                        rates[2][c][i],
                        rates[3][c][i],
                        rates[4][c][i],
                        _CUDARNOps(),
                    )
                r.resolved && _ch_supported(r.value) || _ch_flag!(status, 7, i, c)
                out[c][i] = r.value
            end
        end
    end
    return nothing
end

function _ch_fft_copy!(status, out, field, scale, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            value = field[i]
            scaled = isone(scale) ? value : CUDA.div_rn(value, scale)
            if !isfinite(value) || _CH.SpectralOperators._fft_input_loses_nonzero(value, scaled)
                _ch_flag!(status, 8, i, 1)
            end
            out[i] = ComplexF64(scaled, 0.0)
        end
    end
    return nothing
end

function _ch_fft_mode!(status, out, hat, ik, sign, accumulate, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            value =
                accumulate ? _CH.SpectralOperators._ik_accum_value(out[i], sign, hat[i], ik[i]) :
                _CH.SpectralOperators._ik_store_value(hat[i], ik[i])
            isfinite(value) || _ch_flag!(status, 8, i, 2)
            out[i] = value
        end
    end
    return nothing
end

function _ch_fft_store!(status, out, transformed, normalization, scale, count)
    i = _ch_index()
    if i <= count
        @inbounds begin
            # Match ScaledPlan: normalize the complex inverse first, then take real.
            isfinite(transformed[i]) || _ch_flag!(status, 8, i, 3)
            value = real(transformed[i] * normalization)
            if !_ch_supported(value) ||
               _CH.SpectralOperators._fft_input_loses_nonzero(real(transformed[i]), value)
                _ch_flag!(status, 8, i, 3)
            end
            restored = isone(scale) ? value : CUDA.mul_rn(value, scale)
            _ch_supported(restored) || _ch_flag!(status, 8, i, 3)
            out[i] = restored
        end
    end
    return nothing
end
