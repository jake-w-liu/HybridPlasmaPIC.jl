# Flux-density/covector transforms share one represented metric at each node.
function _shear_flux_to_physical!(out, beta, mesh)
    for I in CartesianIndices(out[1])
        s = mesh._s[I[2]]
        out[1][I] = _shear_linear(beta[1][I], s, beta[2][I])
        out[2][I] = beta[2][I]
        out[3][I] = beta[3][I]
    end
    return out
end
function _shear_covariant!(out, field, mesh)
    for I in CartesianIndices(out[1])
        out[1][I] = field[1][I]
        out[2][I] = _shear_linear(field[2][I], mesh._s[I[2]], field[1][I])
        out[3][I] = field[3][I]
    end
    return out
end

struct ShearedHybridFieldReport{T}
    divergence::T
    tolerance::T
    source_norm::T
    source_risk::Bool
    range_resolved::Bool
end

function _sheared_divergence_receipt(
    st::ShearedPeriodicHybridPIC,
    beta::NTuple{3,Array{T,3}},
) where {T}
    g, hats = st.g, st.work.hats
    magnitude = maximum(a -> maximum(abs, a), beta)
    isfinite(magnitude) || throw(ArgumentError("sheared magnetic flux must be finite"))
    safe = SpectralOperators._fft_safe_magnitude(T, length(beta[1]), one(T))
    scale = SpectralOperators._fft_input_scale(magnitude, safe)
    for c = 1:3
        SpectralOperators._copy_fft_input!(hats[c], beta[c], scale)
        g.plan * hats[c]
        all(isfinite, hats[c]) || throw(ArgumentError("sheared divergence FFT is not finite"))
    end
    Q = _exact_rational
    residual = reference = zero(Rational{BigInt})
    for I in CartesianIndices(hats[1])
        re = im = bound = zero(Rational{BigInt})
        for c = 1:3
            k = Q(g.kvec[c][I[c]])
            z = hats[c][I]
            re += k * Q(real(z))
            im += k * Q(imag(z))
            bound += abs(k) * max(abs(Q(real(z))), abs(Q(imag(z))))
        end
        residual = max(residual, abs(re), abs(im))
        reference = max(reference, bound)
    end
    factor = Q(scale) / length(beta[1])
    residual *= factor
    reference *= factor
    tolerance = Q(st.div_atol) + Q(st.div_rtol) * reference
    residual <= tolerance ||
        throw(ArgumentError("sheared magnetic flux failed its fresh divergence gate"))
    values = map(x -> _from_exact(T, x), (residual, tolerance, reference))
    return ShearedHybridFieldReport{T}(values..., false, all(isfinite, values))
end

function _sheared_complete_fields!(st, beta::NTuple{3,Array{T,3}}) where {T}
    f, w, g = st.work.fields, st.work, st.g
    _shear_flux_to_physical!(f.B, beta, st.mesh)
    _shear_covariant!(w.Bcov, f.B, st.mesh)
    curl!(w.curl, w.Bcov, g)
    all(a -> all(isfinite, a), w.curl) ||
        throw(ArgumentError("sheared logical current is not finite"))
    _shear_flux_to_physical!(f.J, w.curl, st.mesh)
    electron_pressure!(f.pe, f.n, st.closure)
    all(x -> isfinite(x) && (iszero(st.closure.Te) || x > 0), f.pe) || throw(
        ArgumentError(
            "sheared isothermal pressure must materialize without losing positive pressure",
        ),
    )
    gradient!(w.gradient, f.pe, g)
    all(a -> all(isfinite, a), w.gradient) ||
        throw(ArgumentError("sheared pressure gradient is not finite"))
    for I in CartesianIndices(f.pe)
        f.gradp[1][I] = w.gradient[1][I]
        f.gradp[2][I] = _shear_linear(w.gradient[2][I], -st.mesh._s[I[2]], w.gradient[1][I])
        f.gradp[3][I] = w.gradient[3][I]
    end
    foreach(a -> fill!(a, zero(T)), w.lapJ)
    risk = _ohm_field_rows!(f.E, f.ui, f.B, f.J, w.lapJ, f.gradp, f.ninv, zero(T), zero(T))
    all(a -> all(isfinite, a), f.E) || throw(ArgumentError("sheared electric field is not finite"))
    # A finite completed local field is not an original-product differentiated
    # source certificate. All accepted phases are Faraday ready under this rule.
    risk && throw(
        ArgumentError(
            "sheared original-product mapped Faraday recovery is unsupported; canonical Ohm source-risk rejects",
        ),
    )
    # Current has already materialized into f.J, so reuse curl for the admitted
    # Faraday rate, including at init/endpoint/restore preparation.
    _sheared_faraday!(w.curl, st)
    return _sheared_divergence_receipt(st, beta)
end

function _sheared_faraday!(out, st)
    _shear_covariant!(st.work.Ecov, st.work.fields.E, st.mesh)
    curl!(out, st.work.Ecov, st.g)
    _faraday_negate!(out)
    all(a -> all(isfinite, a), out) ||
        throw(ArgumentError("sheared magnetic rate is not representable"))
    return nothing
end
