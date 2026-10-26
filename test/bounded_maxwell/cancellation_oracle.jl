# These fixtures specify source samples and an independent Rational stage oracle;
# they do not call production difference, completion or exact-zero helpers.
function cancellation_step_fixture(T, D, component, electric, profile, periodic, sign)
    epsilon = eps(T)
    j, k = mod1(component + 1, 3), mod1(component + 2, 3)
    power = profile === :division ? 2 : 1
    speed = electric ? epsilon^power : inv(epsilon)^power
    spacing = ntuple(d -> profile === :division && d == j ? one(T) + epsilon : one(T), D)
    g = YeeGrid(
        ntuple(_ -> 4, D),
        map(dx -> T(4) * dx, spacing);
        boundaries = map(b -> b ? :periodic : :pec, periodic),
    )
    m = MaxwellVacuum(g; c = speed)
    source = electric ? m.fields.E : m.fields.B
    index = CartesianIndex(ntuple(_ -> 2, D))
    function pair(axis)
        neighbor =
            CartesianIndex(ntuple(d -> d == axis ? index[d] + (electric ? 1 : -1) : index[d], D))
        electric ? (neighbor, index) : (index, neighbor)
    end
    positive, negative = pair(j)
    source[k][positive] = sign * (profile === :numerator ? one(T) + epsilon : one(T))
    source[k][negative] = sign * (profile === :numerator ? epsilon / 2 : zero(T))
    positive, negative = pair(k)
    source[j][positive] = sign * (profile === :division ? one(T) - epsilon : one(T))
    source[j][negative] = zero(T)
    m, inv(speed) / T(4), index
end

function zero_curl_gradient_fixture(T, D, periodic, anisotropic)
    cells = ntuple(d -> 2d + 2, D)
    spacing = ntuple(d -> anisotropic ? T(2)^(d - 1) : one(T), D)
    g = YeeGrid(
        cells,
        ntuple(d -> T(cells[d]) * spacing[d], D);
        boundaries = map(b -> b ? :periodic : :pec, periodic),
    )
    m = MaxwellVacuum(g)
    function potential(index)
        prod(1:D) do d
            coordinate = periodic[d] ? mod(index[d] - 1, cells[d]) : index[d] - 1
            T(coordinate * (cells[d] - coordinate))
        end
    end
    for component = 1:D, index in CartesianIndices(m.fields.E[component])
        next = CartesianIndex(ntuple(d -> d == component ? index[d] + 1 : index[d], D))
        m.fields.E[component][index] = (potential(next) - potential(index)) / spacing[component]
    end
    m, T(0.25) * maxwell_cfl(m)
end

@noinline cancellation_allocated_step(m, dt) = @allocated step_maxwell!(m, dt)

@testset "Curl cancellation in both Maxwell stages meets the normwise bound" begin
    for T in (Float32, Float64),
        D = 2:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        component in (D == 2 ? (3,) : (1, 2, 3)),
        electric in (false, true),
        profile in (:numerator, :division),
        sign in (-1, 1)

        m, dt, index = cancellation_step_fixture(T, D, component, electric, profile, periodic, sign)
        source = electric ? m.fields.E : m.fields.B
        exact_curl = reference_curl(source, m.g, electric)
        @test !iszero(exact_curl[component][index])
        expectedE, expectedB = reference_step(m, dt)
        @test all(a -> all(isfinite, a), (expectedE..., expectedB...))
        normwise_step_check(m, dt)
        @test m.time[] == dt && m.step[] == 1
    end
end

@testset "Maxwell stage-sum cancellation meets the normwise bound" begin
    for T in (Float32, Float64), sign in (-1, 1), direction in (1, 2, 3)
        epsilon = eps(T)
        g = YeeGrid(
            (2, 2, 2),
            ntuple(d -> T(2) * (d == direction ? one(T) + epsilon : one(T)), 3);
            boundaries = (:periodic, :periodic, :periodic),
        )
        m = MaxwellVacuum(g)
        electric_component = mod1(direction + 1, 3)
        magnetic_component = mod1(direction + 2, 3)
        dt = T(0.25)
        fill!(m.fields.E[electric_component], sign * dt / g.dx[direction])
        for index in CartesianIndices(m.fields.B[magnetic_component])
            m.fields.B[magnetic_component][index] = index[direction] == 2 ? T(sign) : zero(T)
        end
        index = CartesianIndex(ntuple(_ -> 2, 3))
        expected = T(
            Rational{BigInt}(m.fields.E[electric_component][index]) -
            sign * Rational{BigInt}(dt) / Rational{BigInt}(g.dx[direction]),
        )
        @test !iszero(expected) && isfinite(expected)
        normwise_step_check(m, dt)
    end
end

@testset "Exact-zero curl keeps fields exactly and allocates nothing" begin
    for T in (Float32, Float64),
        D = 2:3,
        periodic in Iterators.product(ntuple(_ -> (false, true), D)...),
        anisotropic in (false, true)

        m, dt = zero_curl_gradient_fixture(T, D, periodic, anisotropic)
        expected = map(copy, (m.fields.E..., m.fields.B...))
        @test all(a -> all(iszero, a), reference_curl(m.fields.E, m.g, true))
        for _ = 1:8
            m.time[] = zero(T)
            m.step[] = 0
            cancellation_allocated_step(m, dt)
        end
        for _ = 1:3
            m.time[] = zero(T)
            m.step[] = 0
            @test cancellation_allocated_step(m, dt) == 0
            @test all(isequal.((m.fields.E..., m.fields.B...), expected))
            @test m.time[] == dt && m.step[] == 1
        end
    end
end
