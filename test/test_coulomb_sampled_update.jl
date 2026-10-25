module TestCoulombSampledUpdate
using Test, Random, LinearAlgebra, HybridPlasmaPIC

mutable struct ReplayAngle{T} <: AbstractRNG
    draw::T
    uniforms::Vector{T}
    index::Int
    events::Vector{Symbol}
end
function Random.shuffle!(rng::ReplayAngle, values::AbstractVector)
    push!(rng.events, :shuffle)
    return values
end
function normal(rng::ReplayAngle{T}) where {T}
    push!(rng.events, :normal)
    return rng.draw
end
function uniform(rng::ReplayAngle{T}) where {T}
    push!(rng.events, :uniform)
    rng.index += 1
    return rng.uniforms[rng.index]
end
for T in (Float16, Float32, Float64, BigFloat)
    @eval Random.randn(rng::ReplayAngle{$T}, ::Type{$T}) = normal(rng)
    @eval Random.rand(rng::ReplayAngle{$T}, ::Type{$T}) = uniform(rng)
end

# Independent Rodrigues rotation with a spherical orthogonal basis, followed
# by the physical half kick. No production arithmetic helper is called.
# The random angle uses the original finite rate factors; geometry uses the
# exact input-vector direction. Native azimuth samples are retained.
function reference(a, b, rate, dt, draw, phase, isotropic, polar)
    T = eltype(a)
    native_g = a .- b
    speed = hypot(native_g...)
    sp, cp = sincos((2T(pi)) * phase)
    setprecision(BigFloat, 1024) do
        ga, gb = BigFloat.(a), BigFloat.(b)
        g = ga - gb
        magnitude = norm(g)
        direction = g / magnitude
        radius = hypot(g[1], g[2])
        tangent_direction = if iszero(radius)
            BigFloat[cp, sp, 0]
        else
            longitude = atan(g[2], g[1])
            sl, cl = sincos(longitude)
            polar_cosine, polar_sine = g[3] / magnitude, radius / magnitude
            e_theta = [polar_cosine * cl, polar_cosine * sl, -polar_sine]
            e_phi = [-sl, cl, zero(BigFloat)]
            BigFloat(cp) * e_theta + BigFloat(sp) * e_phi
        end
        axis = cross(direction, tangent_direction)
        if isotropic
            cosine = BigFloat(2polar - one(T))
            sine = BigFloat(sqrt(max(zero(T), one(T) - T(cosine)^2)))
        else
            variance = BigFloat(rate) * BigFloat(dt) / BigFloat(speed)^3
            angle = 2atan(sqrt(variance) * BigFloat(draw))
            sine, cosine = sincos(angle)
        end
        rotation = g * cosine + cross(axis, g) * sine + axis * dot(axis, g) * (1 - cosine)
        kick = (rotation - g) / 2
        return ga + kick, gb - kick
    end
end


# Independent matrix action for the arithmetic contract: all sampled scalar
# operands are captured at particle precision before the high-precision block.
function sampled_reference(a, b, rate, draw, phase, isotropic, polar)
    T = eltype(a)
    g = a .- b
    speed = hypot(g...)
    radius = hypot(g[1], g[2])
    sp, cp = sincos((2T(pi)) * phase)
    tangent = sqrt(rate / ((speed * speed) * speed)) * draw
    cosine = 2polar - one(T)
    sine = sqrt(max(zero(T), one(T) - cosine * cosine))
    setprecision(BigFloat, 1024) do
        ga, gb = BigFloat.(a), BigFloat.(b)
        difference = ga .- gb
        r, m = BigFloat(radius), BigFloat(speed)
        matrix =
            iszero(r) ? hcat(difference, BigFloat[m, 0, 0], BigFloat[0, m, 0]) :
            hcat(
                difference,
                BigFloat[difference[1]*difference[3]/r, difference[2]*difference[3]/r, -r],
                BigFloat[-difference[2]*m/r, difference[1]*m/r, 0],
            )
        sine_big, cosine_big =
            isotropic ? (BigFloat(sine), BigFloat(cosine)) : sincos(2atan(BigFloat(tangent)))
        change = matrix * BigFloat[cosine_big-1, sine_big*BigFloat(cp), sine_big*BigFloat(sp)]
        ga .+ change ./ 2, gb .- change ./ 2
    end
end

@testset "Coulomb sampled-expression accuracy and independent physical rotation" begin
    # Retain the full predeclared seed set, including every cancellation case.
    for T in (Float16, Float32, Float64, BigFloat), seed = 1:256, isotropic in (false, true)
        generator = MersenneTwister(seed)
        a, b = T.(randn(generator, 3)), T.(randn(generator, 3))
        draw, phase, polar = T(randn(generator)), T(rand(generator)), T(rand(generator))
        speed = hypot((a .- b)...)
        rate = speed^3 * (isotropic ? T(2) : T(0.125))
        expected = sampled_reference(a, b, rate, draw, phase, isotropic, polar)
        physical = reference(a, b, rate, one(T), draw, phase, isotropic, polar)
        ps = ParticleSet{1,T}(2)
        fill!(ps.x[1], T(0.25))
        for c = 1:3
            ps.v[c] .= (a[c], b[c])
        end
        rng = ReplayAngle(draw, isotropic ? T[polar, phase] : T[phase], 0, Symbol[])
        collide_coulomb!(ps, rate, one(T); rng, index_work = zeros(Int, 2))
        @test rng.events ==
              (isotropic ? [:shuffle, :uniform, :uniform] : [:shuffle, :normal, :uniform])
        physical_scale = max(norm(BigFloat.(a)), norm(BigFloat.(b)))
        for p = 1:2, c = 1:3
            actual, target = BigFloat(ps.v[c][p]), expected[p][c]
            bound = max(BigFloat(32eps(T)) * abs(target), BigFloat(eps(T(target))) / 2)
            @test abs(actual - target) <= bound
            # Native angle/geometry rounding can amplify relative error in a
            # nearly zero final component; the physical rotation uses a
            # normwise bound, separately from sampled-expression arithmetic.
            @test abs(actual - physical[p][c]) <= BigFloat(32eps(T)) * physical_scale
        end
    end
end
end
