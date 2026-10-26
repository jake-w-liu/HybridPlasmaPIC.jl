# Independent interpolation from component coordinate lattices. This does not
# call the diagnostic's index/average/product helpers. `absolute=true` gives
# the interpolant of |samples|, from which the normwise scale is formed.
function poynting_reference_sample(a, g, electric, component, target; absolute = false)
    D = length(g.n)
    stencils = ntuple(D) do d
        nodal = electric ? component != d : component == d
        offset = nodal ? 0 // 1 : 1 // 2
        coordinate = target[d] - offset
        if !g.periodic[d] && !nodal && (coordinate < 0 || coordinate > g.n[d] - 1)
            # Even extension to the conductor face.
            return ((coordinate < 0 ? 1 : g.n[d], BigInt(1) // BigInt(1)),)
        end
        lower = floor(Int, coordinate)
        fraction = Rational{BigInt}(coordinate - lower)
        index(k) = g.periodic[d] ? mod(k, g.n[d]) + 1 : k + 1
        iszero(fraction) && return ((index(lower), one(fraction)),)
        return ((index(lower), 1 - fraction), (index(lower + 1), fraction))
    end
    total = zero(Rational{BigInt})
    for entries in Iterators.product(stencils...)
        index = CartesianIndex(ntuple(d -> entries[d][1], D))
        sample = Rational{BigInt}(a[index])
        total += (absolute ? abs(sample) : sample) * prod(entry -> entry[2], entries)
    end
    return total
end

function poynting_reference_parts(m, target; absolute = false)
    E = ntuple(c -> poynting_reference_sample(m.fields.E[c], m.g, true, c, target; absolute), 3)
    B = ntuple(c -> poynting_reference_sample(m.fields.B[c], m.g, false, c, target; absolute), 3)
    return E, B
end

function poynting_reference_at(m, target)
    E, B = poynting_reference_parts(m, target)
    return Rational{BigInt}(m.c)^2 .*
           (E[2] * B[3] - E[3] * B[2], E[3] * B[1] - E[1] * B[3], E[1] * B[2] - E[2] * B[1])
end

# Normwise scale c²(|E_j||B_k| + |E_k||B_j|) of each component, with absolute interpolants.
function poynting_scale_at(m, target)
    E, B = poynting_reference_parts(m, target; absolute = true)
    return Rational{BigInt}(m.c)^2 .*
           (E[2] * B[3] + E[3] * B[2], E[3] * B[1] + E[1] * B[3], E[1] * B[2] + E[2] * B[1])
end

# ((exact, scale) per side) per normal; area times the face-cell sums.
function poynting_reference_faces(m)
    D = length(m.g.n)
    return ntuple(D) do normal
        area = prod(d -> d == normal ? one(Rational{BigInt}) : Rational{BigInt}(m.g.dx[d]), 1:D)
        sizes = ntuple(d -> d == normal ? 1 : m.g.n[d], D)
        ntuple(2) do side
            total = zero(Rational{BigInt})
            scale = zero(Rational{BigInt})
            for index in CartesianIndices(sizes)
                target = ntuple(
                    d -> d == normal ? (side == 1 ? 0 // 1 : m.g.n[d] // 1) : index[d] - 1 // 2,
                    D,
                )
                total += poynting_reference_at(m, target)[normal]
                scale += poynting_scale_at(m, target)[normal]
            end
            ((side == 1 ? -1 : 1) * total * area, scale * area)
        end
    end
end

# Normwise diagnostic contract: 8eps(T)*scale (plus N eps for a face sum over
# N cells) and half a subnormal quantum for the final rounding. A nonfinite
# result must be the rounded exact value (outside the floating-point range).
poynting_allowance(::Type{T}, scale, cells = 0) where {T} =
    (8 + cells) * Rational{BigInt}(eps(T)) * scale + Rational{BigInt}(nextfloat(zero(T))) / 2
poynting_close(observed::T, exact, scale, cells = 0) where {T} =
    isfinite(observed) ?
    abs(Rational{BigInt}(observed) - exact) <= poynting_allowance(T, scale, cells) :
    isequal(observed, T(exact))

function poynting_check_all(m)
    D = length(m.g.n)
    observed = yee_poynting(m)
    for index in CartesianIndices(m.g.n)
        target = ntuple(d -> index[d] - 1 // 2, D)
        expected, scale = poynting_reference_at(m, target), poynting_scale_at(m, target)
        for c = 1:3
            @test poynting_close(observed[c][index], expected[c], scale[c])
        end
    end
    faces, expected_faces = yee_poynting_flux(m), poynting_reference_faces(m)
    for d = 1:D, side = 1:2
        cells = prod(e -> e == d ? 1 : m.g.n[e], 1:D)
        @test poynting_close(faces[d][side], expected_faces[d][side]..., cells)
    end
    return observed, faces
end

@testset "Complete interpolated Poynting coordinate oracle" begin
    for T in (Float32, Float64)
        # Subnormal interpolants: the allowance is below one subnormal quantum.
        tiny = nextfloat(zero(T))
        m = MaxwellVacuum(YeeGrid((2,), (one(T),); boundaries = (:periodic,)); c = T(2))
        m.fields.E[2][1] = tiny
        fill!(m.fields.B[3], one(T))
        observed, _ = poynting_check_all(m)
        @test observed[1] == fill(2tiny, 2)
        m2 = MaxwellVacuum(
            YeeGrid((2, 2), (one(T), one(T)); boundaries = (:periodic, :periodic));
            c = T(2),
        )
        m2.fields.E[3][:, 1] .= tiny
        fill!(m2.fields.B[2], one(T))
        _, faces = poynting_check_all(m2)
        @test faces[1][2] == -2tiny
    end
    for T in (Float32, Float64),
        D = 1:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...)

        for profile in
            (:tiny_electric, :tiny_magnetic, :large_speed, :small_speed, :large_area, :ordinary)
            speed =
                profile === :large_speed ? ldexp(one(T), T === Float32 ? 60 : 500) :
                profile === :small_speed ? ldexp(one(T), T === Float32 ? -80 : -550) : T(2)
            amplitude_e =
                profile in (:tiny_electric, :large_speed, :large_area) ? nextfloat(zero(T)) :
                profile === :small_speed ? ldexp(one(T), T === Float32 ? 90 : 600) : one(T)
            amplitude_b =
                profile === :tiny_magnetic ? nextfloat(zero(T)) :
                profile === :small_speed ? amplitude_e : one(T)
            lengths = ntuple(
                d ->
                    profile === :large_area && d > 1 ? ldexp(one(T), T === Float32 ? 70 : 550) :
                    one(T),
                D,
            )
            g = YeeGrid(
                ntuple(_ -> 2, D),
                lengths;
                boundaries = map(p -> p ? :periodic : :pec, periodic),
            )
            m = MaxwellVacuum(g; c = speed)
            for c = 1:3, index in CartesianIndices(m.fields.E[c])
                m.fields.E[c][index] = amplitude_e * T(mod(sum(d -> d * index[d], 1:D) + c, 3) - 1)
            end
            for c = 1:3, index in CartesianIndices(m.fields.B[c])
                m.fields.B[c][index] =
                    amplitude_b * T(mod(sum(d -> (D - d + 1) * index[d], 1:D) + 2c, 3) - 1)
            end
            apply_yee_pec!(m.fields, g)
            arrays = (m.fields.E..., m.fields.B..., m._candidate.E..., m._candidate.B...)
            before = map(copy, arrays)
            poynting_check_all(m)
            @test all(isequal.(arrays, before))
            @test iszero(m.time[]) && iszero(m.step[])
        end
    end
end

@testset "Scaled cross-product cancellation across components and boundaries" begin
    for T in (Float32, Float64),
        D = 1:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        component = 1:3,
        sign in (-1, 1)

        p = precision(T) - 1
        for (electric_scale, magnetic_scale, speed_scale) in
            ((0, 0, p), (p, -p, p), (p, p, 0), (-p, -p, 2p))
            g = YeeGrid(
                ntuple(_ -> 2, D),
                ntuple(_ -> one(T), D);
                boundaries = map(b -> b ? :periodic : :pec, periodic),
            )
            m = MaxwellVacuum(g; c = ldexp(one(T), speed_scale))
            ea, ba = ldexp(one(T), electric_scale), ldexp(one(T), magnetic_scale)
            j, k = mod1(component + 1, 3), mod1(component + 2, 3)
            fill!(m.fields.E[j], sign * ea * (one(T) + eps(T)))
            fill!(m.fields.E[k], sign * ea)
            fill!(m.fields.B[j], ba)
            fill!(m.fields.B[k], ba * (one(T) - eps(T)))
            apply_yee_pec!(m.fields, g)
            arrays = (m.fields.E..., m.fields.B..., m._candidate.E..., m._candidate.B...)
            before = map(copy, arrays)
            poynting_check_all(m)
            if all(periodic)
                # The exact component is -sign after an eps² cancellation.
                expected = poynting_reference_at(m, ntuple(_ -> 1 // 2, D))
                @test expected[component] == -sign
            end
            @test all(isequal.(arrays, before))
            @test iszero(m.time[]) && iszero(m.step[])
        end
    end
end

@testset "Whole-face cancellation retains the normwise face bound" begin
    for T in (Float32, Float64), D = 2:3, normal = 1:D, sign in (-1, 1)
        p = precision(T) - 1
        for (electric_power, magnetic_power, speed_power) in ((0, 0, 0), (p, -p, 0), (0, 0, p))
            tangent = normal < D ? normal + 1 : 1
            magnetic = 6 - normal - tangent
            speed = ldexp(one(T), speed_power)
            lengths = ntuple(d -> d == tangent ? T(2) / (eps(T)^2 * speed^2) : one(T), D)
            m = MaxwellVacuum(
                YeeGrid(ntuple(_ -> 2, D), lengths; boundaries = ntuple(_ -> :periodic, D));
                c = speed,
            )
            ea, ba = ldexp(one(T), electric_power), ldexp(one(T), magnetic_power)
            for index in CartesianIndices(m.fields.E[tangent])
                m.fields.E[tangent][index] =
                    sign * ea * (index[tangent] == 1 ? one(T) + eps(T) : -one(T))
            end
            for index in CartesianIndices(m.fields.B[magnetic])
                m.fields.B[magnetic][index] = ba * (index[tangent] == 1 ? one(T) - eps(T) : one(T))
            end
            arrays = (m.fields.E..., m.fields.B..., m._candidate.E..., m._candidate.B...)
            before = map(copy, arrays)
            expected = poynting_reference_faces(m)
            orientation = tangent == mod1(normal + 1, 3) ? sign : -sign
            @test first.(expected[normal]) == (orientation, -orientation)
            poynting_check_all(m)
            @test all(isequal.(arrays, before))
            @test iszero(m.time[]) && iszero(m.step[])
        end
    end
end

@testset "Subnormal interpolation and complete face integral" begin
    # Subnormal interpolants amplified by c² into the normal range.
    for T in (Float32, Float64)
        g = YeeGrid((4,), (one(T),); boundaries = (:periodic,))
        speed = T === Float32 ? T(1e18) : T(1e154)
        m = MaxwellVacuum(g; c = speed)
        tiny = nextfloat(zero(T))
        fill!(m.fields.E[2], tiny)
        fill!(m.fields.B[3], one(T))
        expected = Rational{BigInt}(speed)^2 * Rational{BigInt}(tiny)
        @test all(x -> poynting_close(x, expected, expected), yee_poynting(m)[1])
        poynting_check_all(m)
    end
    # Overflowing and underflowing face products before area scaling.
    for lengths in ((1.0, 1e-200), (1.0, 1e-200, 1e-200), (1.0, 1e200, 1e100))
        D = length(lengths)
        g = YeeGrid(ntuple(_ -> 4, D), lengths; boundaries = ntuple(_ -> :periodic, D))
        m = MaxwellVacuum(g)
        amplitude = lengths[2] > 1 ? 1e-200 : 1e200
        fill!(m.fields.E[2], amplitude)
        fill!(m.fields.B[3], amplitude)
        expected = Rational{BigInt}(amplitude)^2 * prod(Rational{BigInt}.(lengths[2:end]))
        @test poynting_close(yee_poynting_flux(m)[1][2], expected, expected, 4^(D - 1))
        poynting_check_all(m)
    end
end
