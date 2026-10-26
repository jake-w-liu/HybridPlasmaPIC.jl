# Pointwise stage accuracy is separate from a full-step global-norm comparison.
# Each reference below starts from that stage's actual stored inputs.
function yee_check_stage_accuracy(m, dt, electric)
    source = electric ? m.fields.E : m.fields.B
    original = electric ? m.fields.B : m.fields.E
    before = map(copy, (m.fields.E..., m.fields.B...))
    destination = map(copy, original)
    HybridPlasmaPIC._yee_advance!(destination, original, source, m.g, electric, dt, m.c)
    @test normwise_stage_ratio(destination, original, source, m.g, electric, dt, m.c) <= 1
    @test all(isequal.((m.fields.E..., m.fields.B...), before))
    @test iszero(m.time[]) && iszero(m.step[])
    return destination
end

@testset "Normal nonzero curl cancellation has normwise stage accuracy" begin
    for T in (Float32, Float64),
        D = 2:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        component in (D == 2 ? (3,) : (1, 2, 3)),
        electric in (false, true),
        offset in (-4, -2, -1, 1, 2, 4),
        sign in (-1, 1)

        j, k = mod1(component + 1, 3), mod1(component + 2, 3)
        dx = T(1.3)
        q = one(T) / dx
        neighbor = offset < 0 ? prevfloat(q, -offset) : nextfloat(q, offset)
        g = YeeGrid(
            ntuple(_ -> 4, D),
            ntuple(d -> T(4) * (d == j ? dx : one(T)), D);
            boundaries = map(p -> p ? :periodic : :pec, periodic),
        )
        m = MaxwellVacuum(g; c = electric ? eps(T) : inv(eps(T)))
        source = electric ? m.fields.E : m.fields.B
        index = CartesianIndex(ntuple(_ -> 2, D))
        for (axis, field, amplitude) in ((j, k, T(sign)), (k, j, T(sign) * neighbor))
            positive = electric ? CartesianIndex(ntuple(d -> d == axis ? 3 : 2, D)) : index
            source[field][positive] = amplitude
        end
        apply_yee_pec!(m.fields, g)
        dt = inv(m.c) / T(4)
        @test !iszero(reference_curl(source, m.g, electric)[component][index])
        yee_check_stage_accuracy(m, dt, electric)
    end
end

@testset "Normal stage-sum cancellation across components and walls" begin
    for T in (Float32, Float64),
        D = 1:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        direction = 1:D,
        offset in (-4, -1, 1, 4),
        sign in (-1, 1)

        g = YeeGrid(
            ntuple(_ -> 4, D),
            ntuple(d -> d == direction ? T(5.2) : T(4), D);
            boundaries = map(p -> p ? :periodic : :pec, periodic),
        )
        m = MaxwellVacuum(g)
        dt = T(0.25)
        ec, bc = mod1(direction + 1, 3), mod1(direction + 2, 3)
        q = dt / g.dx[direction]
        old = offset < 0 ? prevfloat(q, -offset) : nextfloat(q, offset)
        fill!(m.fields.E[ec], T(sign) * old)
        for index in CartesianIndices(m.fields.B[bc])
            m.fields.B[bc][index] = index[direction] == 2 ? T(sign) : zero(T)
        end
        apply_yee_pec!(m.fields, g)
        yee_check_stage_accuracy(m, dt, false)
    end
end

function yee_check_diagnostic_accuracy(m)
    arrays = (m.fields.E..., m.fields.B..., m._candidate.E..., m._candidate.B...)
    before = map(copy, arrays)
    poynting_check_all(m)
    @test all(isequal.(arrays, before))
    @test iszero(m.time[]) && iszero(m.step[])
end

@testset "Normal product cancellation retains the normwise diagnostic bound" begin
    for T in (Float32, Float64),
        D = 1:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        component = 1:3,
        offset in (-4, -2, -1, 1, 2, 4),
        sign in (-1, 1)

        a, b = T(1.3), T(1.7)
        product = a * b
        neighbor = offset < 0 ? prevfloat(product, -offset) : nextfloat(product, offset)
        m = MaxwellVacuum(
            YeeGrid(
                ntuple(_ -> 2, D),
                ntuple(_ -> T(2), D);
                boundaries = map(p -> p ? :periodic : :pec, periodic),
            );
            c = inv(eps(T)),
        )
        j, k = mod1(component + 1, 3), mod1(component + 2, 3)
        fill!(m.fields.E[j], sign * a)
        fill!(m.fields.B[k], b)
        fill!(m.fields.E[k], sign * neighbor)
        fill!(m.fields.B[j], one(T))
        apply_yee_pec!(m.fields, m.g)
        yee_check_diagnostic_accuracy(m)
    end
end

@testset "Normal whole-face cancellation across normals and scales" begin
    for T in (Float32, Float64),
        D = 2:3,
        normal = 1:D,
        offset in (-4, -2, -1, 1, 2, 4),
        sign in (-1, 1),
        power in (-3, 0, 3)

        transverse = mod1(normal + 1, D)
        a, b = ldexp(T(1.3), power), ldexp(T(1.7), -power)
        product = a * b
        neighbor = offset < 0 ? prevfloat(product, -offset) : nextfloat(product, offset)
        m = MaxwellVacuum(
            YeeGrid(
                ntuple(_ -> 2, D),
                ntuple(d -> d == transverse ? T(2) / eps(T) : T(2), D);
                boundaries = ntuple(_ -> :periodic, D),
            ),
        )
        # The selected E component is centered along the transverse axis;
        # B's component is also centered there, so the alternating products
        # reach distinct quadrature cells before the complete face sum.
        ec = transverse
        bc = 6 - normal - ec
        for index in CartesianIndices(m.fields.E[ec])
            m.fields.E[ec][index] = sign * (index[transverse] == 1 ? a : -neighbor)
        end
        for index in CartesianIndices(m.fields.B[bc])
            m.fields.B[bc][index] = index[transverse] == 1 ? b : one(T)
        end
        yee_check_diagnostic_accuracy(m)
    end
end
