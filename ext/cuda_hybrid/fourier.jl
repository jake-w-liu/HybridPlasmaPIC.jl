import LinearAlgebra

function _ch_fourier(n)
    input = CUDA.CuArray{ComplexF64,1,CUDA.DeviceMemory}(undef, n)
    accumulator = CUDA.CuArray{ComplexF64,1,CUDA.DeviceMemory}(undef, n)
    forward = CUDA.cuFFT.plan_fft!(input, (1,))
    inverse = CUDA.cuFFT.plan_bfft!(accumulator, (1,))
    return _CUDAHybridFourier(input, accumulator, forward, inverse)
end

function _ch_fft_scale(st, fields)
    magnitude = 0.0
    for field in fields
        value = maximum(abs, field) # One bounded scalar control transfer per component.
        isfinite(value) || throw(ArgumentError("CUDA hybrid Fourier input must be finite"))
        magnitude = max(magnitude, value)
    end
    safe = _CH.SpectralOperators._fft_safe_magnitude(Float64, st.g.n[1], st.g.gain)
    return _CH.SpectralOperators._fft_input_scale(magnitude, safe)
end

function _ch_forward!(st, field, scale)
    _ch_launch!(st, "Fourier input", _ch_fft_copy!, st.g.n[1], st.fourier.input, field, scale)
    LinearAlgebra.mul!(st.fourier.input, st.fourier.forward, st.fourier.input)
    return nothing
end

function _ch_inverse_store!(st, out, scale)
    LinearAlgebra.mul!(st.fourier.accumulator, st.fourier.inverse, st.fourier.accumulator)
    _ch_launch!(
        st,
        "Fourier inverse/store",
        _ch_fft_store!,
        st.g.n[1],
        out,
        st.fourier.accumulator,
        1.0 / st.g.n[1],
        scale,
    )
    return nothing
end

function _ch_gradient!(st, out, field)
    scale = _ch_fft_scale(st, (field,))
    _ch_forward!(st, field, scale)
    _ch_launch!(
        st,
        "Fourier gradient mode",
        _ch_fft_mode!,
        st.g.n[1],
        st.fourier.accumulator,
        st.fourier.input,
        st.g.ik,
        1.0,
        false,
    )
    _ch_inverse_store!(st, out, scale)
    return nothing
end

function _ch_curl!(st, out, field)
    # Release0.1.8 intentionally includes inactive A1 in the common scale.
    scale = _ch_fft_scale(st, field)
    terms = _CH.SpectralOperators._d1_curl_terms()
    for c = 1:3
        if isempty(terms[c])
            fill!(out[c], 0.0)
        else
            fill!(st.fourier.accumulator, zero(ComplexF64))
            for term in terms[c]
                term.axis == 1 ||
                    throw(ArgumentError("CUDA hybrid requires canonical D1 curl geometry"))
                _ch_forward!(st, field[term.input], scale)
                _ch_launch!(
                    st,
                    "Fourier curl mode",
                    _ch_fft_mode!,
                    st.g.n[1],
                    st.fourier.accumulator,
                    st.fourier.input,
                    st.g.ik,
                    Float64(term.sign),
                    true,
                )
            end
            _ch_inverse_store!(st, out[c], scale)
        end
    end
    return nothing
end
