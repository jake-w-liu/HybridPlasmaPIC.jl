module ArrayChecks
using HybridPlasmaPIC, Test, LinearAlgebra
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# A real CPU array with unknown storage ownership exercises the generic backend
# contract. It does not stand in for a GPU implementation.
struct ExternalStorage{T,N} <: AbstractArray{T,N}
    data::Array{T,N}
end
Base.size(a::ExternalStorage) = size(a.data)
Base.IndexStyle(::Type{<:ExternalStorage}) = IndexLinear()
Base.getindex(a::ExternalStorage, i::Int) = a.data[i]
Base.setindex!(a::ExternalStorage, v, i::Int) = (a.data[i] = v)

function exact_square(::Type{S}, n, c) where {S}
    return H._from_exact(S, Q(c.pe0) * (Q(max(n, zero(n))) / Q(c.n0))^2)
end
function accurate_square(pe, n, c)
    expected = exact_square(typeof(pe), n, c)
    isfinite(expected) && !iszero(expected) || return isequal(pe, expected)
    radius = 4max(Q(nextfloat(expected)) - Q(expected), Q(expected) - Q(prevfloat(expected)))
    return isfinite(pe) && abs(Q(pe) - Q(expected)) <= radius
end

@testset "Public pressure preserves array storage and broadcast contracts" begin
    for S in (Float32, Float64), T in (Float32, Float64), C in (Float32, Float64)
        c = PolytropicElectrons(C(0.7), C(1.25), C(2))
        n = reshape(T[0.5, 0.75, 1, 1.5, 2, -1], 2, 3)
        before = copy(n)
        for (source, destination) in (
            (n, zeros(S, 2, 3)),
            (view(n, :, :), view(fill(S(-9), 4, 3), 1:2:4, :)),
            (PermutedDimsArray(n, (2, 1)), PermutedDimsArray(zeros(S, 2, 3), (2, 1))),
            (transpose(n), transpose(zeros(S, 2, 3))),
            (adjoint(n), adjoint(zeros(S, 2, 3))),
            (reshape(view(n, :), 3, 2), reshape(view(zeros(S, 6), :), 3, 2)),
            (
                reinterpret(T, copy(reinterpret((T === Float32 ? UInt32 : UInt64), n))),
                reinterpret(
                    S,
                    copy(reinterpret((S === Float32 ? UInt32 : UInt64), zeros(S, 2, 3))),
                ),
            ),
        )
            @test H._polytropic_host_storage(source)
            @test H._polytropic_host_storage(destination)
            @test electron_pressure!(destination, source, c) === destination
            @test all(accurate_square(p, x, c) for (p, x) in zip(destination, source))
        end
        @test n == before
        destparent = fill(S(-9), 4, 3)
        destination = view(destparent, 1:2:4, :)
        electron_pressure!(destination, n, c)
        @test all(==(S(-9)), view(destparent, 2:2:4, :))
        column = T[0.5, 1.5]
        wide = zeros(S, 2, 3)
        @test electron_pressure!(wide, column, c) === wide
        @test all(accurate_square(wide[i, j], column[i], c) for i = 1:2, j = 1:3)
        row = reshape(T[0.5, 1, 1.5], 1, 3)
        @test electron_pressure!(wide, row, c) === wide
        @test all(accurate_square(wide[i, j], row[j], c) for i = 1:2, j = 1:3)
        empty = S[]
        @test electron_pressure!(empty, T[], c) === empty
        scalar = fill(zero(S))
        @test electron_pressure!(scalar, fill(T(1.5)), c) === scalar
        @test accurate_square(scalar[], T(1.5), c)
        invalid = fill(S(-9), 2, 2)
        @test_throws DimensionMismatch electron_pressure!(invalid, n, c)
        @test all(==(S(-9)), invalid)
    end
end

@testset "Public pressure handles aliases and extreme values through wrappers" begin
    for T in (Float32, Float64)
        c = PolytropicElectrons(T(0.5), one(T), T(2))
        n = T[0.5, 1, 1.5, 2]
        expected = exact_square.(T, n, Ref(c))
        @test electron_pressure!(n, n, c) === n
        @test n == expected
        parent = T[0.5, 1, 1.5, 2, 3]
        source = view(parent, 1:4)
        destination = view(parent, 2:5)
        expected = exact_square.(T, copy(source), Ref(c))
        @test electron_pressure!(destination, source, c) === destination
        @test destination == expected
        k = T === Float32 ? 100 : 600
        for exponent in (-k, k)
            density = ldexp(one(T), exponent)
            extreme = PolytropicElectrons(inv(density), one(T), T(2))
            source = view(fill(density, 2, 3), :, :)
            destination = transpose(zeros(T, 3, 2))
            @test electron_pressure!(destination, source, extreme) === destination
            @test all(==(density), destination)
        end
    end
end

@testset "Other element types and unknown storage retain their arithmetic" begin
    for T in (Float16, BigFloat)
        c = PolytropicElectrons(T(0.7), T(1.25), T(1.5))
        n = T[-1, -0.0, 0.5, 1.5, 2]
        expected = c.pe0 .* (max.(n, zero(T)) ./ c.n0) .^ c.γ
        actual = similar(n)
        @test electron_pressure!(actual, n, c) === actual
        @test isequal(actual, expected)
        raw = T[-1, -0.0, 0.5, 1.5]
        pe = similar(raw)
        @test H._polytropic_unfloored_pressure!(pe, raw, T(0.7), T(2)) === pe
        @test isequal(pe, T(0.7) .* raw .^ T(2))
        for density in raw
            @test isequal(
                H._polytropic_unreferenced_pressure(T, density, T(0.7), T(2)),
                T(0.7) * density^T(2),
            )
        end
    end
    for T in (Float32, Float64)
        k = T === Float32 ? 100 : 600
        c = PolytropicElectrons(ldexp(one(T), -k), one(T), T(2))
        n = ExternalStorage(fill(ldexp(one(T), k), 3))
        pe = ExternalStorage(zeros(T, 3))
        @test !H._polytropic_host_storage(n)
        @test !H._polytropic_host_storage(view(n, :))
        @test electron_pressure!(pe, n, c) === pe
        expected = similar(n.data)
        H._legacy_polytropic_pressure!(expected, n.data, c)
        @test isequal(pe.data, expected)
        @test all(isinf, pe)
    end
end
end
