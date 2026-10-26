module MomentNormalizationTests

using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
struct FixedMomentShape{T} <: H.ShapeFunction
    coefficients::Tuple{T,T}
end
H.width(::FixedMomentShape) = 2
H._stencil1d(shape::FixedMomentShape, s) = (floor(Int, s), shape.coefficients)

# The oracle accumulates complete rational cell expressions. Production flags,
# the production replay, and its intermediate products are never called here.
function exact_moment(ps, g, shape, component = 0, charge = one(eltype(ps.weight)))
    acc = fill(zero(Q), g.n)
    for p in eachindex(ps.weight)
        axes_stencils = map(1:length(g.n)) do d
            H._stencil1d(shape, H._particle_cell_position(ps, g, d, p))
        end
        for offsets in Iterators.product((1:H.width(shape) for _ in g.n)...)
            contribution = Q(ps.weight[p]) * Q(charge)
            component == 0 || (contribution *= Q(ps.v[component][p]))
            destination = ntuple(length(g.n)) do d
                base, coeffs = axes_stencils[d]
                contribution *= Q(coeffs[offsets[d]])
                mod(base + offsets[d] - 1, g.n[d]) + 1
            end
            acc[destination...] += contribution
        end
    end
    volume = prod(Q(dx) for dx in g.dx)
    return map(x -> H._from_exact(eltype(ps.weight), x / volume), acc)
end

function original_deposit!(out, ps, vals, g, shape)
    fill!(out, zero(eltype(out)))
    D = length(g.n)
    for (p, vi) in enumerate(eachindex(vals))
        st = ntuple(d -> H._stencil1d(shape, H._particle_cell_position(ps, g, d, p)), D)
        for c in CartesianIndices(ntuple(_ -> H.width(shape), D))
            o = Tuple(c)
            w = vals[vi]
            for d = 1:D
                w *= st[d][2][o[d]]
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, g.n[d]) + 1, D)
            out[idx...] += w
        end
    end
    out
end
function original_moment(ps, g, shape, component = 0, charge = one(eltype(ps.weight)))
    vals = component == 0 ? ps.weight : ps.weight .* ps.v[component]
    out = zeros(eltype(ps.weight), g.n)
    original_deposit!(out, ps, vals, g, shape)
    out ./= prod(g.dx)
    out .*= charge
    out
end

# Normwise contract: every finite exact cell agrees with the actual value to a
# few ulps of the field scale. A tail cell far below that scale may round
# differently, but real range loss (an underflowed or overflowed field) still
# fails this comparison. Cells whose exact value is non-finite (the true result
# exceeds the representable range) must match exactly.
function normwise_match(actual, expected)
    T = eltype(actual)
    scale = maximum((BigFloat(abs(e)) for e in expected if isfinite(e)); init = zero(BigFloat))
    iszero(scale) && return all(iszero, expected) ? all(iszero, actual) :
           all(i -> isequal(actual[i], expected[i]), eachindex(actual, expected))
    tolerance = 8 * BigFloat(eps(T)) * scale
    return all(
        i ->
            isfinite(expected[i]) ? abs(BigFloat(actual[i]) - BigFloat(expected[i])) <= tolerance :
            isequal(actual[i], expected[i]),
        eachindex(actual, expected),
    )
end

function check_exact(ps, g, shape; density = true, momentum = true, current = true)
    before = map(copy, (ps.x..., ps.v..., ps.weight))
    T = eltype(ps.weight)
    if density
        expected = exact_moment(ps, g, shape)
        n = zeros(T, g.n)
        density!(n, ps, g, shape)
        @test normwise_match(n, expected)
        nt = similar(n)
        density_threaded!(nt, ps, g, shape)
        @test normwise_match(nt, expected)
    end
    for scratch in (false, true)
        work = scratch ? fill(T(7), nparticles(ps)) : nothing
        if momentum
            m = ntuple(_ -> zeros(T, g.n), 3)
            momentum!(m, ps, g, shape; work)
            @test all(normwise_match(m[c], exact_moment(ps, g, shape, c)) for c = 1:3)
        end
        if current
            J = ntuple(_ -> zeros(T, g.n), 3)
            current!(J, ps, g, shape; work)
            @test all(normwise_match(J[c], exact_moment(ps, g, shape, c, ps.q)) for c = 1:3)
        end
        if scratch && (momentum || current)
            m2 = ntuple(_ -> zeros(T, g.n), 3)
            momentum!(m2, ps, g, shape)
            work .= T(NaN)
            m3 = ntuple(_ -> zeros(T, g.n), 3)
            momentum!(m3, ps, g, shape; work)
            @test all(isequal(m2[c], m3[c]) for c = 1:3)
        end
    end
    @test all(isequal(a, b) for (a, b) in zip((ps.x..., ps.v..., ps.weight), before))
end

@testset "MOMENT-001 complete finite expressions" begin
    for T in (Float32, Float64)
        large = T === Float32 ? T(1e25) : T(1e200)
        small = T === Float32 ? T(1e-25) : T(1e-200)
        largew = T === Float32 ? T(1e30) : T(1e300)
        smallw = T === Float32 ? T(1e-30) : T(1e-300)
        for D in (2, 3),
            L in (large, small),
            shape in (NGP(), CIC(), TSC(), FixedMomentShape((T(0.25), T(0.75))))

            g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> L, D))
            ps = ParticleSet{D,T}(3; q = T(0.5))
            for d = 1:D
                ps.x[d] .= g.dx[d] .* T[0.0, 0.25, 0.75]
            end
            ps.weight .= (L == large ? largew : smallw) .* T[1, -1, 0.5]
            for c = 1:3
                ps.v[c] .= T[c, -c, 2c]
            end
            check_exact(ps, g, shape)
        end
        # Intermediate volume range loss followed by recovery to normal size.
        for widths in ((large, large, inv(large)), (small, small, inv(small)))
            g = FourierGrid((2, 2, 2), map(x -> T(2) * x, widths))
            ps = ParticleSet{3,T}(1; q = T(2))
            ps.weight[1] = one(T)
            foreach(x -> fill!(x, zero(T)), ps.x)
            foreach(v -> fill!(v, one(T)), ps.v)
            check_exact(ps, g, NGP())
        end
        # Overflowing partial sums, with and without exact signed cancellation.
        for signed in (false, true), shape in (NGP(), CIC(), TSC())
            g = FourierGrid((2,), (T(2e8),))
            ps = ParticleSet{1,T}(17; q = T(0.25))
            fill!(ps.x[1], zero(T))
            fill!(ps.weight, floatmax(T) / 4)
            signed && (ps.weight[9:16] .*= -1)
            ps.weight[17] = one(T)
            foreach(v -> fill!(v, one(T)), ps.v)
            check_exact(ps, g, shape)
        end
        # Weighted velocity overflow before stencil splitting/normalization.
        for shape in (
            NGP(),
            CIC(),
            TSC(),
            FixedMomentShape((1 // 4, 3 // 4)),
            FixedMomentShape((0, 1)),
            FixedMomentShape((BigFloat(0.25), BigFloat(0.75))),
        )
            g = FourierGrid((2, 2), (T(16), T(16)))
            ps = ParticleSet{2,T}(1; q = T(0.125))
            foreach(x -> fill!(x, T(4)), ps.x)
            ps.weight[1] = floatmax(T) / 2
            foreach(v -> fill!(v, T(4)), ps.v)
            check_exact(ps, g, shape; density = false)
        end
        # Underflow inside a tensor-product stencil, with ordinary finite volume.
        dx, fraction, weight =
            T === Float32 ? (T(1e-15), T(1e-20), T(1e-10)) : (T(1e-75), T(1e-150), T(1e-100))
        g = FourierGrid((2, 2), (2dx, 2dx))
        ps = ParticleSet{2,T}(1; q = T(2))
        for d = 1:2
            ps.x[d][1] = g.dx[d] * fraction
        end
        ps.weight[1] = weight
        foreach(v -> fill!(v, T(2)), ps.v)
        check_exact(ps, g, CIC())
        # An underflowing raw product can be rescued by a large finite charge.
        g = FourierGrid((2,), (T(2),))
        ps = ParticleSet{1,T}(1; q = T === Float32 ? T(1e30) : T(1e300))
        ps.x[1][1] = zero(T)
        ps.weight[1] = nextfloat(zero(T))
        foreach(v -> fill!(v, T(0.25)), ps.v)
        check_exact(ps, g, NGP())
        # Momentum itself can overflow while the current remains finite or zero.
        for q in (T === Float32 ? T(1e-20) : T(1e-200), zero(T))
            ps = ParticleSet{1,T}(16; q)
            fill!(ps.x[1], zero(T))
            fill!(ps.weight, floatmax(T) / 4)
            foreach(v -> fill!(v, one(T)), ps.v)
            check_exact(ps, g, NGP())
        end
        # Empty populations remain zero despite unrepresentable cell volumes.
        g = FourierGrid((2, 2), (small, small))
        ps = ParticleSet{2,T}(0; q = zero(T))
        check_exact(ps, g, TSC())
        # Predividing the weight is invalid: preserve the finite split result.
        g = FourierGrid((2, 2), (one(T), one(T)))
        ps = ParticleSet{2,T}(1)
        foreach(x -> fill!(x, T(0.25)), ps.x)
        ps.weight[1] = floatmax(T) / 2
        foreach(v -> fill!(v, one(T)), ps.v)
        check_exact(ps, g, CIC())
        # Custom stencil coefficients are used as stored, including products
        # narrowed from Float64 into a Float32 output accumulator.
        custom =
            T === Float32 ? FixedMomentShape((1e-30, 1.0 - 1e-30)) :
            FixedMomentShape((1e-200, 1.0 - 1e-200))
        g = FourierGrid((2, 2), (small, small))
        ps = ParticleSet{2,T}(1)
        foreach(x -> fill!(x, zero(T)), ps.x)
        ps.weight[1] = smallw
        foreach(v -> fill!(v, one(T)), ps.v)
        check_exact(ps, g, custom)
    end
end

@testset "MOMENT ordinary arithmetic and nonfinite compatibility" begin
    for T in (Float32, Float64),
        D = 1:3,
        shape in (NGP(), CIC(), TSC(), FixedMomentShape((BigFloat(0.25), BigFloat(0.75)))),
        mode in (:finite, :nanweight, :infweight, :infvelocity, :nancharge)

        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2.3), D))
        ps = ParticleSet{D,T}(9; q = T(0.7))
        rng = MersenneTwister(43 + D)
        for d = 1:D
            ps.x[d] .= rand(rng, T, 9) .* g.L[d]
        end
        ps.weight .= randn(rng, T, 9)
        for c = 1:3
            ps.v[c] .= randn(rng, T, 9)
        end
        mode === :nanweight && (ps.weight[2] = T(NaN))
        mode === :infweight && (ps.weight[2] = T(Inf))
        mode === :infvelocity && (ps.v[2][3] = T(Inf))
        mode === :nancharge && (ps.q = T(NaN))
        n = zeros(T, g.n)
        density!(n, ps, g, shape)
        @test isequal(n, original_moment(ps, g, shape))
        for scratch in (false, true)
            work = scratch ? zeros(T, 9) : nothing
            m = ntuple(_ -> zeros(T, g.n), 3)
            J = map(similar, m)
            momentum!(m, ps, g, shape; work)
            current!(J, ps, g, shape; work)
            @test all(isequal(m[c], original_moment(ps, g, shape, c)) for c = 1:3)
            @test all(isequal(J[c], original_moment(ps, g, shape, c, ps.q)) for c = 1:3)
            if scratch
                work .= T(NaN)
                m4 = ntuple(_ -> zeros(T, g.n), 3)
                momentum!(m4, ps, g, shape; work)
                @test all(isequal(m4[c], m[c]) for c = 1:3)
            end
        end
    end
end

@testset "MOMENT ambient precision does not alter hardware fallback" begin
    for T in (Float32, Float64)
        g = FourierGrid((2,), (T(2e8),))
        ps = ParticleSet{1,T}(17; q = T(0.1))
        fill!(ps.x[1], zero(T))
        fill!(ps.weight, floatmax(T) / 4)
        ps.weight[9:16] .*= -1
        ps.weight[17] = T(1.3)
        foreach(v -> fill!(v, T(0.7)), ps.v)
        for bits in (2, 4, 8, 32, 256)
            setprecision(BigFloat, bits) do
                check_exact(ps, g, NGP())
                @test precision(BigFloat) == bits
            end
        end
    end
end

@testset "Per-species EM current integration and scratch compatibility" begin
    for T in (Float32, Float64), D = 1:3, shape in (NGP(), CIC(), TSC()), extreme in (false, true)
        lengths = ntuple(_ -> extreme ? T(2e8) : T(2.3), D)
        g = FourierGrid(ntuple(_ -> 4, D), lengths)
        es = EMPIC(g, 16; shape)
        ps = ParticleSet{D,T}(16; q = -one(T))
        for d = 1:D
            ps.x[d] .= extreme ? zero(T) : T(0.1d)
        end
        fill!(ps.weight, extreme ? floatmax(T) / 4 : T(0.3))
        for c = 1:3
            fill!(ps.v[c], T(c))
        end
        scratch = zeros(T, g.n)
        work = zeros(T, 16)
        J = ntuple(_ -> zeros(T, g.n), 3)
        for accumulate in (false, true)
            foreach(a -> fill!(a, T(0.5)), J)
            H._add_species_current!(J, es, ps, work, scratch, accumulate)
            for c = 1:3
                oldraw = original_deposit!(zeros(T, g.n), ps, ps.weight .* ps.v[c], g, shape)
                expected =
                    extreme ? exact_moment(ps, g, shape, c, ps.q) : (ps.q / prod(g.dx)) .* oldraw
                accumulate && (expected .+= T(0.5))
                @test isequal(J[c], expected)
            end
            @test isequal(work, ps.weight .* ps.v[3])
            @test isequal(
                scratch,
                original_deposit!(zeros(T, g.n), ps, ps.weight .* ps.v[3], g, shape),
            )
        end
        if D == 1
            legacy = EMPIC1D(g, 16; shape)
            for accumulate in (false, true)
                fill!(legacy.Jy, T(0.5))
                H._add_transverse_Jy!(legacy, ps, work, scratch, accumulate)
                oldraw = original_deposit!(zeros(T, g.n), ps, ps.weight .* ps.v[2], g, shape)
                expected =
                    extreme ? exact_moment(ps, g, shape, 2, ps.q) : (ps.q / prod(g.dx)) .* oldraw
                accumulate && (expected .+= T(0.5))
                @test isequal(legacy.Jy, expected)
                @test isequal(scratch, oldraw)
                @test isequal(work, ps.weight .* ps.v[2])
            end
            oldacc = zero(T)
            for p = 1:16
                oldacc += ps.weight[p] * ps.v[1][p]
            end
            olddc = ps.q * oldacc / g.L[1]
            exactdc = H._from_exact(
                T,
                Q(ps.q) * sum(Q(ps.weight[p]) * Q(ps.v[1][p]) for p = 1:16) / Q(g.L[1]),
            )
            @test H._species_Jx0(ps, g.L[1]) == (extreme ? exactdc : olddc)
        end
    end
end
@testset "EM q/volume coefficient range is checked before recovery" begin
    for T in (Float32, Float64)
        g = FourierGrid((2,), (floatmax(T),))
        ps = ParticleSet{1,T}(1; q = -one(T))
        ps.x[1][1] = zero(T)
        ps.weight[1] = floatmax(T) / 4
        foreach(v -> fill!(v, one(T)), ps.v)
        @test issubnormal(ps.q / g.dx[1])
        es = EMPIC(g, 1; shape = NGP())
        J = ntuple(_ -> zeros(T, g.n), 3)
        scratch = zeros(T, g.n)
        work = zeros(T, 1)
        H._add_species_current!(J, es, ps, work, scratch, false)
        @test all(isequal(J[c], exact_moment(ps, g, NGP(), c, ps.q)) for c = 1:3)
        es1 = EMPIC1D(g, 1; shape = NGP())
        H._add_transverse_Jy!(es1, ps, work, scratch, false)
        @test isequal(es1.Jy, exact_moment(ps, g, NGP(), 2, ps.q))
    end
end

# This check never invokes the production rounding converter. The exact answer
# must lie in the selected float's nearest-even rounding cell.
function rounding_cell_contains(x::T, exact::Q) where {T<:Union{Float32,Float64}}
    isfinite(x) || return false
    lo, hi = prevfloat(x), nextfloat(x)
    (isfinite(lo) && isfinite(hi)) || return false
    lower = (Q(lo) + Q(x)) / 2
    upper = (Q(x) + Q(hi)) / 2
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, x))
    return (exact > lower || (even && exact == lower)) &&
           (exact < upper || (even && exact == upper))
end
@testset "Independent rounding cells and preserved boundaries" begin
    for T in (Float32, Float64)
        # Raw sum overflows, but signed remainder and normalization are finite.
        g = FourierGrid((2,), (T(2e8),))
        ps = ParticleSet{1,T}(17; q = T(0.1))
        fill!(ps.x[1], zero(T))
        fill!(ps.weight, floatmax(T) / 4)
        ps.weight[9:16] .*= -1
        ps.weight[17] = T(1.3)
        foreach(v -> fill!(v, T(0.7)), ps.v)
        n = zeros(T, g.n)
        density!(n, ps, g, NGP())
        nt = similar(n)
        density_threaded!(nt, ps, g, NGP())
        m = ntuple(_ -> similar(n), 3)
        J = map(similar, m)
        momentum!(m, ps, g, NGP())
        current!(J, ps, g, NGP())
        exact = sum(Q(w) for w in ps.weight) / Q(g.dx[1])
        @test rounding_cell_contains(n[1], exact)
        @test rounding_cell_contains(nt[1], exact)
        @test rounding_cell_contains(m[1][1], exact * Q(ps.v[1][1]))
        @test rounding_cell_contains(J[1][1], exact * Q(ps.v[1][1]) * Q(ps.q))
        # Tensor product vanishes before a tiny volume rescues the density.
        dx, fraction, w =
            T === Float32 ? (T(1e-15), T(1e-20), T(1e-10)) : (T(1e-75), T(1e-150), T(1e-100))
        g = FourierGrid((2, 2), (2dx, 2dx))
        ps = ParticleSet{2,T}(1)
        foreach(x -> fill!(x, dx * fraction), ps.x)
        ps.weight[1] = w
        n = zeros(T, g.n)
        density!(n, ps, g, CIC())
        coeffs = ntuple(d -> H._stencil1d(CIC(), H._particle_cell_position(ps, g, d, 1))[2][2], 2)
        exact = Q(w) * Q(coeffs[1]) * Q(coeffs[2]) / prod(Q.(g.dx))
        @test !iszero(exact)
        @test rounding_cell_contains(n[2, 2], exact)
        # The ordinary arithmetic contract deliberately retains in-range cancellation.
        g = FourierGrid((2,), (T(2),))
        ps = ParticleSet{1,T}(3)
        fill!(ps.x[1], zero(T))
        ps.weight .= [T(2)^precision(T), one(T), -T(2)^precision(T)]
        n = zeros(T, g.n)
        density!(n, ps, g, NGP())
        @test n[1] == zero(T)
        @test sum(Q.(ps.weight)) == one(Q)
        # Nonfinite shape operands keep the original propagation behavior.
        for shape in (FixedMomentShape((T(Inf), zero(T))), FixedMomentShape((T(NaN), one(T))))
            density!(n, ps, g, shape)
            @test isequal(n, original_moment(ps, g, shape))
        end
    end
end
@testset "Recovered subnormal volume and narrowed storage witnesses" begin
    for T in (Float32, Float64)
        dx = T === Float32 ? (T(1e-22), T(1e-22), T(1e30)) : (T(1e-160), T(1e-160), T(1e200))
        @test issubnormal(dx[1] * dx[2])
        g = FourierGrid((2, 2, 2), map(x -> 2x, dx))
        ps = ParticleSet{3,T}(1)
        foreach(x -> fill!(x, zero(T)), ps.x)
        ps.weight[1] = one(T)
        n = zeros(T, g.n)
        density!(n, ps, g, NGP())
        @test rounding_cell_contains(n[1, 1, 1], inv(prod(Q.(g.dx))))
    end
    g = FourierGrid((2,), (2.0f-15,))
    ps = ParticleSet{1,Float32}(1)
    ps.x[1][1] = 0
    ps.weight[1] = 1.0f-20
    shape = FixedMomentShape((1e-30, 1.0))
    n = zeros(Float32, g.n)
    @test iszero(Float32(Float64(ps.weight[1]) * shape.coefficients[1]))
    density!(n, ps, g, shape)
    @test rounding_cell_contains(n[1], Q(ps.weight[1]) * Q(shape.coefficients[1]) / Q(g.dx[1]))
end


@testset "Stored BigFloat coefficients survive narrowing and ambient precision" begin
    for T in (Float32, Float64)
        # The coefficient is stored before lowering ambient precision. Its raw
        # native-grid contribution rounds to zero, but cell normalization rescues
        # a finite nonzero result. The oracle retains its exact stored rational.
        shape, coefficient = setprecision(BigFloat, 256) do
            c = BigFloat(nextfloat(zero(T))) * (BigFloat(3) / BigFloat(10))
            FixedMomentShape((c, one(c))), Q(c)
        end
        @test iszero(T(shape.coefficients[1]))
        dx = T === Float32 ? T(1e-20) : T(1e-200)
        g = FourierGrid((2,), (2dx,))
        ps = ParticleSet{1,T}(1; q = T(0.5))
        ps.x[1][1] = zero(T)
        ps.weight[1] = one(T)
        foreach(v -> fill!(v, T(2)), ps.v)
        expected = coefficient / Q(g.dx[1])
        for bits in (2, 4, 8, 32, 256)
            setprecision(BigFloat, bits) do
                n = zeros(T, g.n)
                density!(n, ps, g, shape)
                @test !iszero(n[1])
                @test rounding_cell_contains(n[1], expected)
                density_threaded!(n, ps, g, shape)
                @test rounding_cell_contains(n[1], expected)
                J = ntuple(_ -> zeros(T, g.n), 3)
                current!(J, ps, g, shape)
                @test all(rounding_cell_contains(J[c][1], expected) for c = 1:3)
                @test precision(shape.coefficients[1]) == 256
                @test precision(BigFloat) == bits
            end
        end
    end
end


# The shape interface permits symbolic Real coefficients. Exact replay must not
# introduce a Rational-conversion error where ordinary deposition propagates its
# rounded range result. No symbolic coefficient is silently approximated.
struct SymbolicMomentShape <: H.ShapeFunction end
H.width(::SymbolicMomentShape) = 2
H._stencil1d(::SymbolicMomentShape, s) = (floor(Int, s), (π, 1 - π))
@testset "Symbolic Real shape coefficients preserve legacy arithmetic" begin
    @test sum((π, 1 - π)) == 1
    for T in (Float32, Float64), extreme in (false, true)
        g = FourierGrid((2,), (T(2e8),))
        ps = ParticleSet{1,T}(2; q = -one(T))
        fill!(ps.x[1], zero(T))
        fill!(ps.weight, extreme ? floatmax(T) / 4 : T(0.3))
        for c = 1:3
            fill!(ps.v[c], T(c))
        end
        shape = SymbolicMomentShape()
        n = zeros(T, g.n)
        density!(n, ps, g, shape)
        @test isequal(n, original_moment(ps, g, shape))
        density_threaded!(n, ps, g, shape)
        @test isequal(n, original_moment(ps, g, shape))
        for scratch in (false, true)
            work = scratch ? zeros(T, 2) : nothing
            m = ntuple(_ -> zeros(T, g.n), 3)
            J = map(similar, m)
            momentum!(m, ps, g, shape; work)
            current!(J, ps, g, shape; work)
            @test all(isequal(m[c], original_moment(ps, g, shape, c)) for c = 1:3)
            @test all(isequal(J[c], original_moment(ps, g, shape, c, ps.q)) for c = 1:3)
            if scratch
                work .= T(NaN)
                m4 = ntuple(_ -> zeros(T, g.n), 3)
                momentum!(m4, ps, g, shape; work)
                @test all(isequal(m4[c], m[c]) for c = 1:3)
            end
        end
        es = EMPIC(g, 2; shape)
        J = ntuple(_ -> zeros(T, g.n), 3)
        work = zeros(T, 2)
        raw = zeros(T, g.n)
        H._add_species_current!(J, es, ps, work, raw, false)
        @test all(
            isequal(
                J[c],
                (ps.q / prod(g.dx)) .*
                original_deposit!(zeros(T, g.n), ps, ps.weight .* ps.v[c], g, shape),
            ) for c = 1:3
        )
        es1 = EMPIC1D(g, 2; shape)
        H._add_transverse_Jy!(es1, ps, work, raw, false)
        @test isequal(
            es1.Jy,
            (ps.q / prod(g.dx)) .*
            original_deposit!(zeros(T, g.n), ps, ps.weight .* ps.v[2], g, shape),
        )
    end
end

end
