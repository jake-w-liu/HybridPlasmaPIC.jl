# Spectra.jl — power spectrum (from diagnostics.jl)

"""
    power_spectrum(field, g) -> (k, P)

1-D spatial power spectrum |f̂(k)|² of a periodic field along axis 1, summed over
the remaining axes; `k` are the non-negative angular wavenumbers.
"""
function power_spectrum(field::AbstractArray{T,D}, g::FourierGrid{D,T}) where {T,D}
    _require_grid_array(:field, field, g)
    fh = fft(field)
    n1 = g.n[1]
    nk = n1 ÷ 2 + 1
    P = zeros(T, nk)
    @inbounds for I in CartesianIndices(fh)
        m = Tuple(I)[1] - 1                # 0-based index along axis 1
        km = m <= n1 ÷ 2 ? m : n1 - m       # fold to non-negative
        km + 1 <= nk && (P[km+1] += abs2(fh[I]))
    end
    k = [T(2π) * (m) / g.L[1] for m = 0:nk-1]
    return k, P
end

function power_spectrum(field::AbstractArray, g::FourierGrid)
    _require_grid_array(:field, field, g)
    throw(ArgumentError("field eltype $(eltype(field)) does not match grid eltype"))
end

# ---------------------------------------------------------------- pressure–strain
