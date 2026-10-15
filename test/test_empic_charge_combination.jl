module EMPICChargeCombinationTests

using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function exact_density(ps, g, shape)
    cells = fill(zero(Q), g.n)
    for p in eachindex(ps.weight)
        stencils =
            ntuple(d -> H._stencil1d(shape, H._particle_cell_position(ps, g, d, p)), length(g.n))
        for offsets in Iterators.product((1:H.width(shape) for _ in g.n)...)
            value = Q(ps.weight[p])
            destination = ntuple(length(g.n)) do d
                value *= Q(stencils[d][2][offsets[d]])
                mod(stencils[d][1] + offsets[d] - 1, g.n[d]) + 1
            end
            cells[destination...] += value
        end
    end
    return cells ./ prod(Q(dx) for dx in g.dx)
end

# Independent nearest-even rounding cells, including subnormals and overflow.
function in_rounding_cell(value::T, exact::Q) where {T<:Union{Float32,Float64}}
    isnan(value) && return false
    maximum = Q(floatmax(T))
    overflow = maximum + (maximum - Q(prevfloat(floatmax(T)))) / 2
    isinf(value) && return signbit(value) ? exact <= -overflow : exact >= overflow
    before, after = prevfloat(value), nextfloat(value)
    lower = isfinite(before) ? (Q(before) + Q(value)) / 2 : -overflow
    upper = isfinite(after) ? (Q(after) + Q(value)) / 2 : overflow
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, value))
    return (exact > lower || (even && exact == lower)) &&
           (exact < upper || (even && exact == upper))
end

function previous_charge!(rho, ne, es, e, ions)
    if es.mobile && ions !== nothing
        fill!(rho, zero(eltype(rho)))
        for ps in (ions, e)
            density!(ne, ps, es.g, es.shape)
            for I in eachindex(rho)
                rho[I] += ps.q * ne[I]
            end
        end
    else
        density!(ne, e, es.g, es.shape)
        for I in eachindex(rho)
            rho[I] = es.n0 - ne[I]
        end
    end
    return rho
end

function particle_view(ps::ParticleSet{D,T}) where {D,T}
    ParticleSet{D,T}(
        map(a -> view(a, :), ps.x),
        map(a -> view(a, :), ps.v),
        view(ps.weight, :),
        view(ps.id, :),
        view(ps.tag, :),
        ps.q,
        ps.m,
    )
end

function check_complete(es, e, ions)
    T = eltype(e.weight)
    particles = ions === nothing ? (e,) : (e, ions)
    storage = Tuple(a for ps in particles for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
    before = map(copy, storage)
    exact = exact_density(e, es.g, es.shape)
    if es.mobile && ions !== nothing
        exact .= Q(e.q) .* exact .+ Q(ions.q) .* exact_density(ions, es.g, es.shape)
    else
        exact .= Q(es.n0) .- exact
    end
    rho, ne = zeros(T, es.g.n), zeros(T, es.g.n)
    @test H._charge_density!(rho, ne, es, e, ions) === rho
    @test all(in_rounding_cell(rho[I], exact[I]) for I in eachindex(rho))
    scratch = similar(ne)
    @test density!(scratch, e, es.g, es.shape) === scratch
    @test isequal(ne, scratch)
    @test all(isequal(a, b) for (a, b) in zip(storage, before))
    return rho
end

@testset "Complete charge across species, background and cell range" begin
    for T in (Float32, Float64),
        D = 1:3,
        legacy in (false, true),
        shape in (NGP(), CIC(), TSC()),
        views in (false, true)

        legacy && D != 1 && continue
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> T(2), D))
        e = ParticleSet{D,T}(16; q = -one(T))
        ions = ParticleSet{D,T}(17; q = one(T))
        for ps in (e, ions)
            foreach(x -> fill!(x, zero(T)), ps.x)
            fill!(ps.weight, floatmax(T) / 4)
        end
        ions.weight[end] = one(T)
        # The surviving ion charge occupies a different stencil from the
        # cancelling large populations, so replay must use its own positions.
        for d = 1:D
            ions.x[d][end] = T(1.25) * g.dx[d]
        end
        views && (e = particle_view(e); ions = particle_view(ions))
        es = legacy ? EMPIC1D(g, 16; mobile = true, shape) : EMPIC(g, 16; mobile = true, shape)
        check_complete(es, e, ions)
        # Sixteen half-maximum weights overflow even the largest 3D TSC
        # coefficient; the three-particle control below stays in range.
        bg =
            legacy ? EMPIC1D(g, 16; n0 = floatmax(T), shape) : EMPIC(g, 16; n0 = floatmax(T), shape)
        e = ParticleSet{D,T}(16; q = -one(T))
        foreach(x -> fill!(x, zero(T)), e.x)
        fill!(e.weight, floatmax(T) / 2)
        check_complete(bg, views ? particle_view(e) : e, nothing)
        # Equal finite populations cancel before an unrepresentable volume.
        scale = T === Float32 ? T(1e-25) : T(1e-200)
        tiny = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> scale, D))
        e = ParticleSet{D,T}(2; q = -one(T))
        ions = ParticleSet{D,T}(2; q = one(T))
        for ps in (e, ions)
            foreach(x -> fill!(x, zero(T)), ps.x)
            fill!(ps.weight, floatmax(T) / 2)
        end
        es = legacy ? EMPIC1D(tiny, 2; mobile = true, shape) : EMPIC(tiny, 2; mobile = true, shape)
        @test all(iszero, check_complete(es, e, ions))
    end
end

@testset "Complete charge subnormal midpoint and empty populations" begin
    for T in (Float32, Float64), legacy in (false, true)
        u = nextfloat(zero(T))
        g = FourierGrid((2,), (T(16),))
        e = ParticleSet{1,T}(1; q = -one(T))
        ions = ParticleSet{1,T}(1; q = one(T))
        fill!(e.x[1], zero(T))
        fill!(ions.x[1], zero(T))
        e.weight[1] = -3u
        ions.weight[1] = 3u
        es =
            legacy ? EMPIC1D(g, 1; mobile = true, shape = NGP()) :
            EMPIC(g, 1; mobile = true, shape = NGP())
        @test check_complete(es, e, ions) == T[u, 0]
        e.weight[1] = 4u
        bg = legacy ? EMPIC1D(g, 1; n0 = u, shape = NGP()) : EMPIC(g, 1; n0 = u, shape = NGP())
        @test check_complete(bg, e, nothing) == T[0, u]
        empty_e = ParticleSet{1,T}(0; q = -one(T))
        empty_i = ParticleSet{1,T}(0; q = one(T))
        @test all(iszero, check_complete(es, empty_e, empty_i))
        @test all(==(u), check_complete(bg, empty_e, nothing))
        check_complete(es, empty_e, ions)
        check_complete(es, e, empty_i)
    end
end

@testset "Ordinary charge and nonfinite propagation retain baseline order" begin
    for T in (Float32, Float64),
        D = 1:3,
        legacy in (false, true),
        mobile in (false, true),
        shape in (NGP(), CIC(), TSC()),
        mode in (:ordinary, :nan_e, :inf_i)

        legacy && D != 1 && continue
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2.3), D))
        rng = MersenneTwister(2319 + D)
        e = ParticleSet{D,T}(13; q = -one(T))
        ions = ParticleSet{D,T}(9; q = one(T))
        for ps in (e, ions)
            for d = 1:D
                ps.x[d] .= rand(rng, T, length(ps.weight)) .* g.L[d]
            end
            ps.weight .= randn(rng, T, length(ps.weight))
        end
        mode === :nan_e && (e.weight[2] = T(NaN))
        mode === :inf_i && (ions.weight[2] = T(Inf))
        es =
            legacy ? EMPIC1D(g, 13; mobile, shape, n0 = T(0.7)) :
            EMPIC(g, 13; mobile, shape, n0 = T(0.7))
        rho, ne = zeros(T, g.n), zeros(T, g.n)
        expected, expected_ne = similar(rho), similar(ne)
        selected = mobile ? ions : nothing
        previous_charge!(expected, expected_ne, es, e, selected)
        H._charge_density!(rho, ne, es, e, selected)
        @test isequal(rho, expected)
        @test isequal(ne, expected_ne)
    end
end

@testset "Neutral extreme populations initialize and step both EMPIC families" begin
    for T in (Float32, Float64), D = 1:3, legacy in (false, true), n_sub in (1, 2)
        legacy && D != 1 && continue
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        e = ParticleSet{D,T}(16; q = -one(T))
        ions = ParticleSet{D,T}(16; q = one(T))
        for ps in (e, ions)
            foreach(x -> fill!(x, zero(T)), ps.x)
            fill!(ps.weight, floatmax(T) / 4)
        end
        es =
            legacy ? EMPIC1D(g, 16; mobile = true, mi = ions.m, shape = NGP(), n_sub) :
            EMPIC(g, 16; mobile = true, mi = ions.m, shape = NGP(), n_sub)
        init_empic!(es, e, ions)
        @test all(iszero, es.rho_n)
        @test all(iszero, es.rho_np1)
        # Equal stationary species have zero total charge and current, even
        # when their separate densities overflow. Check the public substeps.
        for dt in (T(1) / 16, T(1) / 8)
            @test step_empic!(es, e, ions, dt) === es
            @test all(iszero, es.rho_n)
            @test all(iszero, es.rho_np1)
            fields = legacy ? (es.Ex, es.Ey, es.Bz, es.Jx, es.Jy) : (es.E..., es.B..., es.J...)
            @test all(a -> all(iszero, a), fields)
            @test all(ps -> all(a -> all(iszero, a), (ps.x..., ps.v...)), (e, ions))
            @test charge_conservation_residual(es, dt / n_sub) == zero(T)
        end
        @test es.time[] == T(3) / 16
        @test es.step[] == 2
    end
end

@testset "In-range extreme background retains ordinary rounding" begin
    for T in (Float32, Float64), D in (2, 3)
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> T(2), D))
        e = ParticleSet{D,T}(3; q = -one(T))
        foreach(x -> fill!(x, zero(T)), e.x)
        fill!(e.weight, floatmax(T) / 2)
        es = EMPIC(g, 3; n0 = floatmax(T), shape = TSC())
        rho, ne = zeros(T, g.n), zeros(T, g.n)
        state = H._density_with_range!(ne, e, g, TSC())
        @test !state.lost
        expected = es.n0 .- ne
        H._charge_density!(rho, ne, es, e, nothing)
        @test isequal(rho, expected)
    end
end

struct ChargeCoefficientShape{T} <: H.ShapeFunction
    coefficients::Tuple{T,T}
end
H.width(::ChargeCoefficientShape) = 2
H._stencil1d(shape::ChargeCoefficientShape, s) = (floor(Int, s), shape.coefficients)

struct ChargeSymbolicShape <: H.ShapeFunction end
H.width(::ChargeSymbolicShape) = 2
H._stencil1d(::ChargeSymbolicShape, s) = (floor(Int, s), (π, 1 - π))

struct ChargeStencilFailure <: Exception end
struct ChargeThrowingShape <: H.ShapeFunction end
H.width(::ChargeThrowingShape) = 1
function H._stencil1d(::ChargeThrowingShape, s)
    s > one(s) && throw(ChargeStencilFailure())
    return (floor(Int, s), (one(s),))
end

@testset "Charge coefficients retain represented precision and eligibility" begin
    for T in (Float32, Float64), legacy in (false, true)
        g = FourierGrid((2,), (T(2),))
        u = nextfloat(zero(T))
        coefficient = setprecision(BigFloat, 256) do
            BigFloat(u) * (BigFloat(3) / BigFloat(10))
        end
        for shape in (
            ChargeCoefficientShape((coefficient, zero(coefficient))),
            ChargeCoefficientShape((1 // 3, 0 // 1)),
        )
            e = ParticleSet{1,T}(1; q = -one(T))
            ions = ParticleSet{1,T}(1; q = one(T))
            fill!(e.x[1], zero(T))
            fill!(ions.x[1], zero(T))
            weight = shape.coefficients[1] isa BigFloat ? one(T) : u
            e.weight[1] = -weight
            ions.weight[1] = weight
            es = legacy ? EMPIC1D(g, 1; mobile = true, shape) : EMPIC(g, 1; mobile = true, shape)
            for bits in (53, 127, 1024)
                setprecision(BigFloat, bits) do
                    @test check_complete(es, e, ions) == T[u, 0]
                    @test precision(BigFloat) == bits
                end
            end
        end
        for shape in (
            ChargeSymbolicShape(),
            ChargeCoefficientShape((T(Inf), zero(T))),
            ChargeCoefficientShape((T(NaN), one(T))),
        )
            e = ParticleSet{1,T}(2; q = -one(T))
            ions = ParticleSet{1,T}(2; q = one(T))
            for ps in (e, ions)
                fill!(ps.x[1], zero(T))
                fill!(ps.weight, floatmax(T))
            end
            for mobile in (false, true)
                es = legacy ? EMPIC1D(g, 2; mobile, shape) : EMPIC(g, 2; mobile, shape)
                rho, ne = zeros(T, g.n), zeros(T, g.n)
                expected, expected_ne = similar(rho), similar(ne)
                selected = mobile ? ions : nothing
                previous_charge!(expected, expected_ne, es, e, selected)
                @test H._charge_density!(rho, ne, es, e, selected) === rho
                @test isequal(rho, expected)
                @test isequal(ne, expected_ne)
            end
        end
    end
end

@testset "Charge stencil errors propagate without changing particles" begin
    for T in (Float32, Float64), legacy in (false, true), failing in (:electron, :ion)
        g = FourierGrid((2,), (T(2),))
        e = ParticleSet{1,T}(1; q = -one(T))
        ions = ParticleSet{1,T}(1; q = one(T))
        fill!(e.x[1], failing === :electron ? T(1.5) : T(0.5))
        fill!(ions.x[1], failing === :ion ? T(1.5) : T(0.5))
        before = (copy(e.x[1]), copy(ions.x[1]))
        shape = ChargeThrowingShape()
        es = legacy ? EMPIC1D(g, 1; mobile = true, shape) : EMPIC(g, 1; mobile = true, shape)
        rho, ne = zeros(T, g.n), zeros(T, g.n)
        @test_throws ChargeStencilFailure H._charge_density!(rho, ne, es, e, ions)
        @test e.x[1] == before[1]
        @test ions.x[1] == before[2]
    end
end

end
