module EMPICCurrentCombinationTests

using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

struct CurrentFixedShape{T} <: H.ShapeFunction
    coefficients::Tuple{T,T}
end
H.width(::CurrentFixedShape) = 2
H._stencil1d(shape::CurrentFixedShape, s) = (floor(Int, s), shape.coefficients)
struct CurrentSymbolicShape <: H.ShapeFunction end
H.width(::CurrentSymbolicShape) = 2
H._stencil1d(::CurrentSymbolicShape, s) = (floor(Int, s), (π, 1 - π))
struct CurrentFailureShape <: H.ShapeFunction end
H.width(::CurrentFailureShape) = 1
function H._stencil1d(::CurrentFailureShape, s)
    s > typeof(s)(0.5) && error("injected midpoint stencil failure")
    return H._stencil1d(NGP(), s)
end

function species(D, T, N, q, mass, views)
    views || return ParticleSet{D,T}(N; q, m = mass)
    storage(::Type{S}) where {S} = view(zeros(S, N + 2), 2:(N+1))
    return ParticleSet{D,T}(
        ntuple(_ -> storage(T), D),
        ntuple(_ -> storage(T), 3),
        storage(T),
        storage(UInt64),
        storage(UInt32),
        q,
        mass,
    )
end
function setup(
    D,
    T,
    legacy;
    counts = (32, 33),
    shape = NGP(),
    views = false,
    width = one(T),
    mobile = true,
    cells = 3,
)
    g = FourierGrid(ntuple(_ -> cells, D), ntuple(_ -> T(cells) * width, D))
    e = species(D, T, counts[1], -one(T), one(T), views)
    ions = species(D, T, counts[2], one(T), T(1836), views)
    es = legacy ? EMPIC1D(g, counts[1]; mobile, shape) : EMPIC(g, counts[1]; mobile, shape)
    legacy ? H._ensure_ion_buffers!(es, counts[2]) : H._ensure_empic_ion_buffers!(es, counts[2])
    emid = legacy ? (es.work,) : es.mide
    imid = legacy ? (es.worki,) : es.midi
    for (ps, mid) in ((e, emid), (ions, imid))
        for d = 1:D
            fill!(ps.x[d], T(1.75) * width)
            fill!(mid[d], zero(T))
        end
        for c = 1:3
            fill!(ps.v[c], (one(T), T(-2), T(0.5))[c])
        end
        fill!(ps.weight, one(T))
    end
    return (; es, e, ions, emid, imid)
end

function deposit(s)
    ions = s.es.mobile ? s.ions : nothing
    return s.es isa EMPIC1D ? (H._deposit_Jy_at_midpoints!(s.es, s.e, ions),) :
           H._deposit_current_at_midpoints!(s.es, s.e, ions)
end

# Independent exact grid expression: enumerate the stored midpoint stencils,
# then combine both species before normalization or conversion. No production
# range flags or exact-moment accumulation helper is used.
function exact_current(s, component)
    g, shape = s.es.g, s.es.shape
    acc = fill(zero(Q), g.n)
    for (ps, mid) in ((s.e, s.emid), (s.ions, s.imid))
        ps === s.ions && !s.es.mobile && continue
        for p = 1:nparticles(ps)
            stencils = ntuple(length(g.n)) do d
                H._stencil1d(shape, mod(mid[d][p], g.L[d]) / g.dx[d])
            end
            for offsets in Iterators.product((1:H.width(shape) for _ in g.n)...)
                term = Q(ps.q) * Q(ps.weight[p]) * Q(ps.v[component][p])
                index = ntuple(length(g.n)) do d
                    base, coefficients = stencils[d]
                    term *= Q(coefficients[offsets[d]])
                    mod(base + offsets[d] - 1, g.n[d]) + 1
                end
                acc[index...] += term
            end
        end
    end
    volume = prod(Q.(g.dx))
    return acc ./ volume
end
function rounding_cell(x::T, exact::Q) where {T<:Union{Float32,Float64}}
    isfinite(x) || return false
    lo, hi = prevfloat(x), nextfloat(x)
    isfinite(lo) && isfinite(hi) || return false
    lower, upper = (Q(lo) + Q(x)) / 2, (Q(x) + Q(hi)) / 2
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, x))
    return (exact > lower || (even && exact == lower)) &&
           (exact < upper || (even && exact == upper))
end
function check_range(s)
    before = map(
        copy,
        (
            s.e.x...,
            s.ions.x...,
            s.e.v...,
            s.ions.v...,
            s.e.weight,
            s.ions.weight,
            s.emid...,
            s.imid...,
        ),
    )
    components = s.es isa EMPIC1D ? (2,) : (1, 2, 3)
    exact = map(c -> exact_current(s, c), components)
    output = deposit(s)
    for (actual, expected) in zip(output, exact)
        @test isequal(actual, map(x -> H._from_exact(eltype(actual), x), expected))
        @test all(rounding_cell(x, y) for (x, y) in zip(actual, expected))
    end
    after = (
        s.e.x...,
        s.ions.x...,
        s.e.v...,
        s.ions.v...,
        s.e.weight,
        s.ions.weight,
        s.emid...,
        s.imid...,
    )
    @test all(isequal(a, b) for (a, b) in zip(before, after))
    if s.es isa EMPIC1D
        @test isequal(s.es.Eyp, s.e.weight .* s.e.v[2])
        @test isequal(s.es.Eyi, s.ions.weight .* s.ions.v[2])
        @test isequal(s.es.Bzp, s.e.x[1])
        @test isequal(s.es.Bzi, s.ions.x[1])
        exactdc =
            sum(
                (
                    Q(ps.q) * Q(ps.weight[p]) * Q(ps.v[1][p]) for ps in (s.e, s.ions) for
                    p = 1:nparticles(ps)
                );
                init = zero(Q),
            ) / Q(s.es.g.L[1])
        value = H._combined_Jx0(s.es, s.e, s.ions)
        @test rounding_cell(value, exactdc)
        H._esirkepov_Jx!(s.es, s.e, s.ions, one(eltype(s.e)))
        @test all(==(value), s.es.Jx)
    else
        @test isequal(s.es.worke, s.e.weight .* s.e.v[3])
        @test isequal(s.es.worki, s.ions.weight .* s.ions.v[3])
    end
end

@testset "Complete species currents use their own midpoints" begin
    for T in (Float32, Float64), D = 1:3, legacy in (false, true)
        legacy && D != 1 && continue
        for shape in (
                NGP(),
                CIC(),
                TSC(),
                CurrentFixedShape((T(0.25), T(0.75))),
                CurrentFixedShape((1 // 3, 2 // 3)),
            ),
            views in (false, true)

            s = setup(D, T, legacy; shape, views)
            fill!(s.e.weight, floatmax(T) / 4)
            fill!(s.ions.weight, floatmax(T) / 4)
            s.ions.weight[end] = T(1.3)
            for d = 1:D
                s.imid[d][end] = T(0.375) * s.es.g.dx[d]
            end
            check_range(s)
        end
    end
end

@testset "Combined current subnormal and tie-to-even boundaries" begin
    for T in (Float32, Float64), legacy in (false, true)
        u = nextfloat(zero(T))
        s = setup(1, T, legacy; counts = (1, 1), width = T(8))
        foreach(v -> fill!(v, -one(T)), s.e.v)
        foreach(v -> fill!(v, one(T)), s.ions.v)
        s.e.weight[1] = 3u
        s.ions.weight[1] = 3u
        output = deposit(s)
        @test all(first(a) == u for a in output)
        @test exact_current(s, 2)[1] == 3Q(u) / 4
        @test rounding_cell(first(first(output)), 3Q(u) / 4)
        s.e.weight[1] = 6u
        s.ions.weight[1] = 2u
        foreach(v -> fill!(v, -one(T)), s.ions.v)
        output = deposit(s)
        @test all(iszero(first(a)) for a in output)
        @test exact_current(s, 2)[1] == Q(u) / 2
        @test rounding_cell(first(first(output)), Q(u) / 2)
        if legacy
            s.e.weight[1] = 9u
            s.ions.weight[1] = 9u
            foreach(v -> fill!(v, one(T)), s.ions.v)
            @test H._combined_Jx0(s.es, s.e, s.ions) == u
            s.e.weight[1] = 18u
            s.ions.weight[1] = 6u
            foreach(v -> fill!(v, -one(T)), s.ions.v)
            @test iszero(H._combined_Jx0(s.es, s.e, s.ions))
        end
    end
end

@testset "Range status survives cell volume and retained coefficient precision" begin
    for T in (Float32, Float64)
        for width in (T === Float32 ? T(1e25) : T(1e200), T === Float32 ? T(1e-25) : T(1e-200))
            s = setup(2, T, false; counts = (1, 1), width)
            w =
                width > one(T) ? (T === Float32 ? T(1e30) : T(1e300)) :
                (T === Float32 ? T(1e-30) : T(1e-300))
            s.e.weight[1] = w
            s.ions.weight[1] = nextfloat(w)
            check_range(s)
        end
        shape = setprecision(BigFloat, 256) do
            CurrentFixedShape((BigFloat(1) / 3, BigFloat(2) / 3))
        end
        s = setup(1, T, false; shape)
        fill!(s.e.weight, floatmax(T) / 4)
        fill!(s.ions.weight, floatmax(T) / 4)
        s.ions.weight[end] = T(1.3)
        s.imid[1][end] = T(0.375)
        # At two-bit ambient precision the original third-component sums incur
        # no range event. The archived implementation also returns zero here;
        # ordinary rounded cancellation remains part of the contract.
        setprecision(BigFloat, 2) do
            @test all(iszero, deposit(s)[3])
        end
        # Force a native weighted-velocity overflow before BigFloat stencil
        # arithmetic, so every component must retain the stored exact operands.
        for ps in (s.e, s.ions), c = 1:3
            fill!(ps.v[c], (T(8), T(-16), T(32))[c])
        end
        @test all(isinf(s.e.weight[1] * s.e.v[c][1]) for c = 1:3)
        expected = setprecision(BigFloat, 256) do
            map(c -> exact_current(s, c), (1, 2, 3))
        end
        for bits in (2, 4, 8, 64, 256)
            setprecision(BigFloat, bits) do
                output = deposit(s)
                @test all(
                    all(rounding_cell(x, y) for (x, y) in zip(a, b)) for
                    (a, b) in zip(output, expected)
                )
                @test precision(shape.coefficients[1]) == 256
                @test precision(BigFloat) == bits
            end
        end
    end
end

@testset "Current combination keeps empty and nonfinite behavior" begin
    for T in (Float32, Float64), legacy in (false, true), mobile in (false, true)
        s = setup(1, T, legacy; counts = (0, 0), mobile)
        @test all(all(iszero, a) for a in deposit(s))
        s = setup(1, T, legacy; counts = (1, 1), mobile)
        s.e.v[2][1] = T(Inf)
        s.ions.v[2][1] = T(Inf)
        output = deposit(s)
        @test mobile ? isnan(legacy ? output[1][1] : output[2][1]) :
              (legacy ? output[1][1] : output[2][1]) == -T(Inf)
        @test all(==(T(1.75)), s.e.x[1])
        @test all(==(T(1.75)), s.ions.x[1])
    end
    for T in (Float32, Float64), legacy in (false, true), extreme in (false, true)
        s = setup(1, T, legacy; counts = (32, 32), shape = CurrentSymbolicShape())
        if extreme
            fill!(s.e.weight, floatmax(T) / 4)
            fill!(s.ions.weight, floatmax(T) / 4)
        end
        output = deposit(s)
        @test extreme ? any(any(isnan, a) for a in output) : all(all(iszero, a) for a in output)
        @test all(==(T(1.75)), s.e.x[1])
        @test all(==(T(1.75)), s.ions.x[1])
    end
end

@testset "Midpoint borrowing restores state and propagates failures" begin
    for legacy in (false, true), failing_electron in (false, true)
        s = setup(1, Float64, legacy; counts = (1, 1), shape = CurrentFailureShape())
        s.emid[1][1] = failing_electron ? 0.75 : 0.25
        s.imid[1][1] = failing_electron ? 0.25 : 0.75
        before = map(copy, (s.e.x..., s.ions.x..., s.emid..., s.imid...))
        @test_throws ErrorException deposit(s)
        @test all(
            isequal(a, b) for (a, b) in zip(before, (s.e.x..., s.ions.x..., s.emid..., s.imid...))
        )
        @test_throws ErrorException H._exact_midpoint_current(
            failing_electron ? s.e : s.ions,
            failing_electron ? s.emid : s.imid,
            s.es.g,
            s.es.shape,
            2,
        )
        @test all(
            isequal(a, b) for (a, b) in zip(before, (s.e.x..., s.ions.x..., s.emid..., s.imid...))
        )
    end
    s = setup(2, Float64, false; counts = (1, 1))
    resize!(s.emid[2], 0)
    before = map(copy, (s.e.x..., s.emid...))
    @test_throws DimensionMismatch H._swap_positions_with_midpoints!(s.e, s.emid)
    @test all(isequal(a, b) for (a, b) in zip(before, (s.e.x..., s.emid...)))
end

end
