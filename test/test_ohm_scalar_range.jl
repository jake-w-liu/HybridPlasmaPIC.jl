using HybridPlasmaPIC, Test, LinearAlgebra

# Independent physical vector expression over exact stored operands. J, gradp
# and inverse density are represented inputs; this is not an exact FFT oracle.
function _ohm_range_vector_oracle(u, B, J, ninv, η, ηH, lap, gradp)
    Q = Rational{BigInt}
    qu, qb, qj, qlap, qgp = map(v -> Q.(collect(v)), (u, B, J, lap, gradp))
    return -cross(qu, qb) + cross(qj, qb) * Q(ninv) + Q(η) * qj - Q(ηH) * qlap - Q(ninv) * qgp
end

function _ohm_range_components(u, B, J, ninv, η, ηH, lap, gradp)
    return ntuple(3) do c
        a, b = ((2, 3), (3, 1), (1, 2))[c]
        HybridPlasmaPIC._ohm_scalar_field_value(
            u[a], u[b], B[a], B[b], J[a], J[b], ninv, η, J[c], ηH, lap[c], gradp[c],
        )
    end
end

@testset "scalar Ohm represented-factor range completion" begin
    for T in (Float32, Float64)
        z, o = zero(T), one(T)
        zeros3 = (z, z, z)
        large = ldexp(o, exponent(floatmax(T)) ÷ 2 + 1)
        small = ldexp(o, T === Float32 ? -100 : -600)
        invsmall = ldexp(o, T === Float32 ? 100 : 600)
        cases = (
            # Individually overflowing motional products cancel exactly.
            ((z, large, large), (z, large, large), zeros3, o, z, z, zeros3, zeros3),
            # The original Hall product is below the subnormal range, but its
            # complete inverse-density-scaled field is normal and nonzero.
            (zeros3, (z, z, small), (z, small, z), invsmall, z, z, zeros3, zeros3),
            (zeros3, (z, z, -small), (z, small, z), invsmall, z, z, zeros3, zeros3),
            # Motional/Hall and resistive terms must combine before rounding.
            ((z, large, z), (z, z, large), (z, large, z), o, z, z, zeros3, zeros3),
            (zeros3, zeros3, (large, large, large), o, large, large, (large, large, large), zeros3),
            (zeros3, zeros3, (large, z, z), large, large, z, zeros3, (large, z, z)),
            # A complete overflowing result retains native signed infinity.
            ((z, large, z), (z, z, large), zeros3, o, z, z, zeros3, zeros3),
        )
        for args in cases
            exact = _ohm_range_vector_oracle(args...)
            actual = _ohm_range_components(args...)
            @test all(isequal(actual[c], T(exact[c])) for c = 1:3)
        end
        @test iszero(small * small)
        @test _ohm_range_components(cases[2]...)[1] == small
        @test isinf(large * large)

        # The ordinary path preserves the released scalar operation order,
        # including the absence of a pressure subtraction off spatial axes.
        args = T.((0.3, -0.2, 0.7, -0.4, 0.9, 0.6, 1.2, 0.03, -0.8, 0.02, 0.5))
        v1, v2, b1, b2, j1, j2, ni, η, jc, ηH, lap = args
        for pressure in (nothing, T(0.17))
            native = -(v1 * b2 - v2 * b1) + (j1 * b2 - j2 * b1) * ni + η * jc - ηH * lap
            pressure === nothing || (native -= pressure * ni)
            @test isequal(HybridPlasmaPIC._ohm_scalar_field_value(args..., pressure), native)
        end
        @test isnan(HybridPlasmaPIC._ohm_scalar_field_value(T(NaN), z, z, o, z, z, o, z, z, z, z, nothing))

        # The public scalar and CGL Ohm owners reach the shared repair for
        # uniform data in every supported spatial dimension.
        for D = 1:3, closure in (IsothermalElectrons(z), CGLElectrons(z, z, o, o))
            n = ntuple(_ -> 4, D)
            g = FourierGrid(n, ntuple(_ -> T(2π), D))
            f = HybridFields{D,T}(n; anisotropic = closure isa CGLElectrons)
            fill!(f.n, o)
            for c in (2, 3)
                fill!(f.ui[c], large)
                fill!(f.B[c], large)
            end
            @test ohms_law!(f, HybridModel(closure), g) === f.E
            @test all(all(iszero, component) for component in f.E)
        end
    end
end
