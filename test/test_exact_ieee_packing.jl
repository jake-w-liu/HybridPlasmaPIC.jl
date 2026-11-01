module ExactIEEEPackingTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
dyadic(n, e) = e >= 0 ? (BigInt(n) << e) // BigInt(1) : BigInt(n) // (BigInt(1) << -e)

@noinline function observed(operation::F, raw::Int32, flush::Bool) where {F}
    operation() # Warm this concrete operation before changing hardware controls.
    saved, zeros = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
    try
        Base.Rounding.setrounding_raw(Float64, raw)
        Base.set_zero_subnormals(flush)
        return operation(), Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
    finally
        Base.Rounding.setrounding_raw(Float64, saved)
        Base.set_zero_subnormals(zeros)
    end
end

function ieee_checks(::Type{T}, ::Type{U}, p, bias) where {T,U}
    f = p - 1
    signmask = one(U) << (8sizeof(T) - 1)
    normal = one(U) << f
    onecode, twocode = U(bias) << f, U(bias + 1) << f
    infcode = ((one(U) << (8sizeof(T) - p)) - one(U)) << f
    minq = 1 - bias - f
    # IEEE format boundaries, derived from precision/exponent fields.
    inputs = (
        zero(U), one(U), normal - one(U), normal, onecode, twocode, infcode - one(U),
    )
    overflowtie = dyadic((BigInt(1) << (p + 1)) - 1, bias - p)
    outputs = (
        (dyadic(1, minq), one(U)),
        (dyadic(1, minq - 1), zero(U)),
        (dyadic((BigInt(1) << p) - 1, minq - 1), normal),
        (dyadic((BigInt(1) << (p + 1)) - 1, -p), twocode),
        (dyadic((BigInt(1) << p) - 1, bias - f), infcode - one(U)),
        (overflowtie - dyadic(1, bias - p - 1), infcode - one(U)),
        (overflowtie, infcode),
    )
    saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
    for hardware in (RoundNearest, RoundDown, RoundUp, RoundToZero), flush in (false, true)
        raw = Base.Rounding.to_fenv(hardware)
        for negative in (false, true)
            sign = negative ? signmask : zero(U)
            for code in inputs
                field, fraction = Int(code >> f), code & (normal - one(U))
                n = BigInt(fraction) + (field == 0 ? 0 : BigInt(1) << f)
                magnitude = dyadic(n, (field == 0 ? 1 : field) - bias - f)
                exact = negative ? -magnitude : magnitude
                value = reinterpret(T, code | sign)
                @test observed(() -> H._exact_rational(value), raw, flush) ==
                      (exact, raw, flush)
            end
            @test observed(
                () -> reinterpret(U, H._pack_rounded_ieee(T, 0, minq, negative)),
                raw,
                flush,
            ) == (sign, raw, flush)
            for (magnitude, code) in outputs
                exact = negative ? -magnitude : magnitude
                scaled, square = exact / 8, exact * exact
                @test observed(
                    () -> reinterpret(U, H._from_exact(T, exact)), raw, flush,
                ) == (code | sign, raw, flush)
                @test observed(
                    () -> reinterpret(U, H._from_exact(T, scaled, Int128(3))), raw, flush,
                ) == (code | sign, raw, flush)
                @test observed(
                    () -> reinterpret(
                        U,
                        H._sqrt_from_exact(T, square, Int128(0), RoundNearest, negative),
                    ),
                    raw,
                    flush,
                ) == (code | sign, raw, flush)
            end
            midpoint = 1 + dyadic(1, -p)
            for mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)
                away = mode === RoundFromZero || (negative ? mode === RoundDown : mode === RoundUp)
                expected = (onecode + U(away)) | sign
                @test observed(
                    () -> reinterpret(
                        U,
                        H._sqrt_from_exact(T, midpoint^2, Int128(0), mode, negative),
                    ),
                    raw,
                    flush,
                ) == (expected, raw, flush)
            end
            if T !== Float16
                n = negative ? Int128(-1) : Int128(1)
                @test observed(
                    () -> reinterpret(U, H._from_exact(T, n, minq)), raw, flush,
                ) == (one(U) | sign, raw, flush)
                @test observed(
                    () -> reinterpret(U, H._from_exact(T, n, UInt64(1), minq)), raw, flush,
                ) == (one(U) | sign, raw, flush)
            end
        end
    end
    @test (Base.Rounding.rounding_raw(Float64), get_zero_subnormals()) == saved
    @test H._exact_rational(T(Inf)) == BigInt(1) // BigInt(0)
    @test H._exact_rational(-T(Inf)) == BigInt(-1) // BigInt(0)
    @test_throws ArgumentError H._exact_rational(T(NaN))
end

@noinline fixed_dyadic(::Type{T}, n, e) where {T} = H._from_exact(T, n, e)
@noinline fixed_ratio(::Type{T}, n, d, e) where {T} = H._from_exact(T, n, d, e)
function fixed_bytes(::Type{T}) where {T}
    n, d, e = Int128(1234567), UInt64(37), -19
    fixed_dyadic(T, n, e)
    fixed_ratio(T, n, d, e)
    return @allocated(fixed_dyadic(T, n, e)), @allocated(fixed_ratio(T, n, d, e))
end

@testset "Canonical IEEE exact conversion boundaries" begin
    # Binary16/32/64 precision and maximum unbiased exponent are format constants.
    for (T, U, p, bias) in (
        (Float16, UInt16, 11, 15),
        (Float32, UInt32, 24, 127),
        (Float64, UInt64, 53, 1023),
    )
        ieee_checks(T, U, p, bias)
    end
    @test H._exact_rational(Int128(-17)) == BigInt(-17) // BigInt(1)
    @test H._exact_rational(7 // 11) == BigInt(7) // BigInt(11)
    stored = setprecision(BigFloat, 512) do
        BigFloat(1) + ldexp(BigFloat(1), -400)
    end
    for bits in (24, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                @test H._exact_rational(stored) == 1 + dyadic(1, -400)
                @test precision(BigFloat) == bits && rounding(BigFloat) === mode
            end
        end
    end
    for T in (Float32, Float64)
        fixed_bytes(T)
        @test fixed_bytes(T) == (0, 0)
    end
end
end
