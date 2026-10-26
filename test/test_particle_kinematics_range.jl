module ParticleKinematicsRangeTests
# Normwise range contract of the particle kinematics (kinematics_range.jl):
# certified native inputs retain historical bits; repaired IEEE priming uses
# rounding-cell certification before wider replay. Range-flagged results agree
# with independent exact references to a few ulps of the operation scale, and
# public steppers reject unresolved timesteps.
using HybridPlasmaPIC, Test, Random, LinearAlgebra, Serialization
import MPI
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# ------------------------------------------------------------------ tolerance

function binade(s::Q)
    n, d = numerator(s), denominator(s)
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    (e >= 0 ? n < d << e : n << -e < d) && (e -= 1)
    return e
end

# Spacing of T at magnitude sqrt(s2): the ulp of an operation whose 2-norm
# scale is sqrt(s2). Subnormal scales use the subnormal quantum.
function ulp_sqrt(::Type{T}, s2::Q) where {T}
    emin = exponent(floatmin(T))
    e = iszero(s2) ? emin : max(fld(binade(s2), 2), emin)
    return Q(2)^(e - precision(T) + 1)
end
ulp_big(::Type{T}, s::BigFloat) where {T} =
    BigFloat(
        2,
    )^(
        max(iszero(s) ? exponent(floatmin(T)) : exponent(s), exponent(floatmin(T))) - precision(T) +
        1
    )

norm2(x) = sum(abs2, x)

# Componentwise check against an exact result with an absolute normwise
# tolerance. A component whose exact value rounds beyond T must be that ±Inf.
function normwise_ok(actual::NTuple{N,T}, exact, tol::Q) where {N,T}
    return all(1:N) do i
        rounded = H._from_exact(T, Q(exact[i]))
        isfinite(rounded) ? isfinite(actual[i]) && abs(Q(actual[i]) - exact[i]) <= tol :
        isequal(actual[i], rounded)
    end
end

# ------------------------------------------------------------------ independent references

# Exact Cayley form of the Boris rotation (algebraically the two-step Boris).
function boris_exact(v, E, B, qm, dt)
    v, E, B = map(x -> Q.(collect(x)), (v, E, B))
    h = Q(qm) * Q(dt) / 2
    t = h .* B
    w = v .+ h .* E
    t2 = dot(t, t)
    rotated = ((1 - t2) .* w .+ 2dot(t, w) .* t .+ 2cross(w, t)) ./ (1 + t2)
    return Tuple(rotated .+ h .* E), max(norm2(v), norm2(h .* E))
end

# Relativistic Boris in 4096-bit BigFloat from the exact Lorentz gap.
function boris_rel_reference(v, E, B, qm, dt, c)
    setprecision(BigFloat, 4096) do
        gap = (Q(c)^2 - sum(x -> Q(x)^2, v)) / Q(c)^2
        u = BigFloat.(collect(v)) ./ sqrt(BigFloat(gap))
        h = BigFloat(Q(qm) * Q(dt) / 2)
        hE = h .* BigFloat.(collect(E))
        c2 = BigFloat(c)^2
        w = u .+ hE
        t = h .* BigFloat.(collect(B)) ./ sqrt(1 + dot(w, w) / c2)
        t2 = dot(t, t)
        u2 = ((1 - t2) .* w .+ 2dot(t, w) .* t .+ 2cross(w, t)) ./ (1 + t2) .+ hE
        γf = sqrt(1 + dot(u2, u2) / c2)
        out = u2 ./ γf
        return Tuple(out), max(norm(out), max(norm(u), norm(hE)) / γf)
    end
end

# Exact forward-Euler half-kick u ∓ h*qm*(E + u×B) and its term scale.
function euler_exact(u, E, B, h, qm, sign)
    u, E, B = map(x -> Q.(collect(x)), (u, E, B))
    hq = Q(h) * Q(qm)
    out = u .+ sign * hq .* (E .+ cross(u, B))
    return Tuple(out), max(norm2(u), hq^2 * norm2(E), hq^2 * norm2(u) * norm2(B))
end

# Relativistic momentum back-up, 4096-bit reference.
function rel_prime_reference(v, E, B, qm, h, c)
    setprecision(BigFloat, 4096) do
        gap = (Q(c)^2 - sum(x -> Q(x)^2, v)) / Q(c)^2
        g = BigFloat.(collect(v)) ./ sqrt(BigFloat(gap))
        a =
            BigFloat(qm) .*
            (BigFloat.(collect(E)) .+ cross(BigFloat.(collect(v)), BigFloat.(collect(B))))
        um = g .- BigFloat(h) .* a
        gm = sqrt(1 + dot(um, um) / BigFloat(c)^2)
        out = um ./ gm
        return Tuple(out), max(norm(out), max(norm(g), abs(BigFloat(h)) * norm(a)) / gm)
    end
end

function rel_ok(actual::NTuple{3,T}, reference, scale::BigFloat) where {T}
    tol = 4ulp_big(T, scale)
    return all(i -> isfinite(actual[i]) && abs(BigFloat(actual[i]) - reference[i]) <= tol, 1:3)
end

# ------------------------------------------------------------------ historical formulas

# Copies of the bodies on main before the range repair (bitwise references).
function old_boris(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
    h = qm * dt / 2
    vmx = vx + h * Ex
    vmy = vy + h * Ey
    vmz = vz + h * Ez
    tx = h * Bx
    ty = h * By
    tz = h * Bz
    t2 = tx * tx + ty * ty + tz * tz
    f = 2 / (1 + t2)
    sx = f * tx
    sy = f * ty
    sz = f * tz
    vpx = vmx + (vmy * tz - vmz * ty)
    vpy = vmy + (vmz * tx - vmx * tz)
    vpz = vmz + (vmx * ty - vmy * tx)
    vnx = vmx + (vpy * sz - vpz * sy)
    vny = vmy + (vpz * sx - vpx * sz)
    vnz = vmz + (vpx * sy - vpy * sx)
    return (vnx + h * Ex, vny + h * Ey, vnz + h * Ez)
end

function old_boris_rel(vx::T, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt, c) where {T}
    c2 = c * c
    βx = vx / c
    βy = vy / c
    βz = vz / c
    β2 = βx * βx + βy * βy + βz * βz
    ginit = one(T) / sqrt(one(T) - β2)
    ux = ginit * vx
    uy = ginit * vy
    uz = ginit * vz
    h = qm * dt / 2
    umx = ux + h * Ex
    umy = uy + h * Ey
    umz = uz + h * Ez
    γm = sqrt(one(T) + (umx * umx + umy * umy + umz * umz) / c2)
    tx = h * Bx / γm
    ty = h * By / γm
    tz = h * Bz / γm
    t2 = tx * tx + ty * ty + tz * tz
    f = 2 / (1 + t2)
    sx = f * tx
    sy = f * ty
    sz = f * tz
    upx = umx + (umy * tz - umz * ty)
    upy = umy + (umz * tx - umx * tz)
    upz = umz + (umx * ty - umy * tx)
    unx = umx + (upy * sz - upz * sy)
    uny = umy + (upz * sx - upx * sz)
    unz = umz + (upx * sy - upy * sx)
    ux2 = unx + h * Ex
    uy2 = uny + h * Ey
    uz2 = unz + h * Ez
    γf = sqrt(one(T) + (ux2 * ux2 + uy2 * uy2 + uz2 * uz2) / c2)
    return (ux2 / γf, uy2 / γf, uz2 / γf)
end

function old_prime!(v, Ep, Bp, qm, h, np, rel, c)
    vx, vy, vz = v
    ex, ey, ez = Ep
    bx, by, bz = Bp
    if rel
        c2 = c * c
        o = one(c)
        for p = 1:np
            v0x, v0y, v0z = vx[p], vy[p], vz[p]
            β2 = (v0x / c) * (v0x / c) + (v0y / c) * (v0y / c) + (v0z / c) * (v0z / c)
            g0 = o / sqrt(o - β2)
            ax = qm * (ex[p] + (v0y * bz[p] - v0z * by[p]))
            ay = qm * (ey[p] + (v0z * bx[p] - v0x * bz[p]))
            az = qm * (ez[p] + (v0x * by[p] - v0y * bx[p]))
            umx = g0 * v0x - h * ax
            umy = g0 * v0y - h * ay
            umz = g0 * v0z - h * az
            gm = sqrt(o + (umx * umx + umy * umy + umz * umz) / c2)
            vx[p] = umx / gm
            vy[p] = umy / gm
            vz[p] = umz / gm
        end
    else
        for p = 1:np
            ux, uy, uz = vx[p], vy[p], vz[p]
            vx[p] = ux - h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
            vy[p] = uy - h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
            vz[p] = uz - h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
        end
    end
    return v
end

function old_predict!(vp, v, Ep, Bp, qm, h, np)
    vx, vy, vz = v
    ex, ey, ez = Ep
    bx, by, bz = Bp
    for p = 1:np
        ux, uy, uz = vx[p], vy[p], vz[p]
        vp[1][p] = ux + h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
        vp[2][p] = uy + h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
        vp[3][p] = uz + h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
    end
    return vp
end

function old_prime_empic1d!(v, Exp, Eyp, Bzp, qm, h, np, rel, c)
    vx, vy, vz = v
    if rel
        c2 = c * c
        o = one(c)
        for p = 1:np
            bz = Bzp[p]
            v0x, v0y, v0z = vx[p], vy[p], vz[p]
            β2 = (v0x / c) * (v0x / c) + (v0y / c) * (v0y / c) + (v0z / c) * (v0z / c)
            g0 = o / sqrt(o - β2)
            ax = qm * (Exp[p] + v0y * bz)
            ay = qm * (Eyp[p] - v0x * bz)
            umx = g0 * v0x - h * ax
            umy = g0 * v0y - h * ay
            umz = g0 * v0z
            gm = sqrt(o + (umx * umx + umy * umy + umz * umz) / c2)
            vx[p] = umx / gm
            vy[p] = umy / gm
            vz[p] = umz / gm
        end
    else
        for p = 1:np
            bz = Bzp[p]
            ux, uy = vx[p], vy[p]
            vx[p] = ux - h * qm * (Exp[p] + uy * bz)
            vy[p] = uy - h * qm * (Eyp[p] - ux * bz)
        end
    end
    return v
end

# ------------------------------------------------------------------ ordinary inputs

# Signed magnitudes in [2^-span, 2^span) with exact (signed) zeros, well inside
# the range of every type (the Float16 kick-scale threshold is 256*floatmin ≈ 0.016).
function ordinary(rng, ::Type{T}, span) where {T}
    r = rand(rng)
    r < 0.02 && return zero(T)
    r < 0.04 && return -zero(T)
    return T((rand(rng, Bool) ? -1 : 1) * (1 + rand(rng)) * exp2(rand(rng, (-span):(span-1))))
end
spans(::Type{Float16}) = 1
spans(::Type{Float32}) = 8
spans(::Type{Float64}) = 16
const NORDINARY = 100_000

# This classification uses only the physical reference and the documented
# output-conditioning threshold. It does not inspect candidate flags/differences.
function check_rel_prime_compatibility(actual, historical, start, Ep, Bp, qm, h, c)
    certified = repaired = unchanged_failures = rounding_failures = 0
    for p in eachindex(start[1])
        v, E, B = map(arrays -> ntuple(i -> arrays[i][p], 3), (start, Ep, Bp))
        reference, _ = rel_prime_reference(v, E, B, qm, h, c)
        setprecision(BigFloat, 4096) do
            beta2 = sum(abs2, reference) / BigFloat(c)^2
            # This fixture's classification must be far from the uncertainty
            # of its independent reference arithmetic. Boundary cases belong
            # to dedicated interval/rounding controls instead of guesswork.
            @assert abs(beta2 - BigFloat(0.5)) > 64eps(BigFloat)
            observed = ntuple(i -> actual[i][p], 3)
            if beta2 <= BigFloat(0.5)
                certified += 1
                unchanged_failures += !isequal(observed, ntuple(i -> historical[i][p], 3))
            else
                repaired += 1
                rounding_failures += !isequal(observed, eltype(actual[1]).(reference))
            end
        end
    end
    @test certified + repaired == length(start[1])
    @test unchanged_failures == 0
    @test rounding_failures == 0
    @info "Priming physical conditioning coverage" type = eltype(actual[1]) certified repaired
end

@testset "Kinematics preserve native bits or correctly round conditioned priming" begin
    for T in (Float16, Float32, Float64)
        rng = Xoshiro(0x6b696e)
        s = spans(T)
        boris = rel = 0
        flagged = 0
        for _ = 1:NORDINARY
            v = ntuple(_ -> ordinary(rng, T, s), 3)
            E = ntuple(_ -> ordinary(rng, T, s), 3)
            B = ntuple(_ -> ordinary(rng, T, s), 3)
            qm = ordinary(rng, T, s)
            dt = ordinary(rng, T, s)
            out, parts = H._boris_kick_parts(v..., E..., B..., qm, dt)
            flagged += !H._boris_kick_in_range(out, parts, v..., qm, dt)
            boris +=
                !isequal(
                    H.boris_kick(v..., E..., B..., qm, dt),
                    old_boris(v..., E..., B..., qm, dt),
                )
            # Subluminal relativistic inputs.
            c = T(1 + 3rand(rng))
            w = ntuple(_ -> T(c * (rand(rng) - 0.5) / 2), 3)
            rel +=
                !isequal(
                    H._boris_rel(w..., E..., B..., qm, dt, c),
                    old_boris_rel(w..., E..., B..., qm, dt, c),
                )
        end
        @test flagged == 0
        @test boris == 0
        @test rel == 0

        n = NORDINARY
        gen(f) = ntuple(_ -> [f() for _ = 1:n], 3)
        v = gen(() -> ordinary(rng, T, s))
        Ep = gen(() -> ordinary(rng, T, s))
        Bp = gen(() -> ordinary(rng, T, s))
        qm, h = T(-1.25), T(0.0625)
        for rel_branch in (false, true)
            c = T(4)
            start = rel_branch ? map(a -> clamp.(a, -T(1), T(1)), v) : v
            new, old = map(copy, start), map(copy, start)
            H._prime_leapfrog!(new, Ep, Bp, qm, h, n, rel_branch, c)
            old_prime!(old, Ep, Bp, qm, h, n, rel_branch, c)
            if rel_branch
                check_rel_prime_compatibility(new, old, start, Ep, Bp, qm, h, c)
            else
                @test isequal(new, old)
            end
            new, old = map(copy, start), map(copy, start)
            H._prime_empic1d!(new, Ep[1], Ep[2], Bp[3], qm, h, n, rel_branch, c)
            old_prime_empic1d!(old, Ep[1], Ep[2], Bp[3], qm, h, n, rel_branch, c)
            if rel_branch
                z = fill(zero(T), n)
                planar_E, planar_B = (Ep[1], Ep[2], z), (z, z, Bp[3])
                check_rel_prime_compatibility(new, old, start, planar_E, planar_B, qm, h, c)
            else
                @test isequal(new, old)
            end
        end
        new, old = map(similar, v), map(similar, v)
        H._predict_half_kick!(new, v, Ep, Bp, qm, h, n)
        old_predict!(old, v, Ep, Bp, qm, h, n)
        @test isequal(new, old)
        hq = h * qm
        @test isequal(H._es_prime_component.(v[1], Ep[1], h, qm, hq), (v[1] .- h * qm .* Ep[1]))
        x, dt = Ep[1], T(0.375)
        @test isequal(H._particle_coordinate.(x, v[1], dt), x .+ dt .* v[1])
        @test isequal(H._particle_coordinate.(x, v[1], dt / 2), x .+ (dt / 2) .* v[1])
    end
end

# ------------------------------------------------------------------ Boris witnesses

function boris_cases(::Type{T}) where {T}
    large, larger, small, strong =
        T === Float16 ? (12, 5, -12, 9) :
        T === Float32 ? (100, 40, -100, 80) : (800, 300, -1000, 600)
    p(k) = ldexp(one(T), k)
    z = (zero(T), zero(T), zero(T))
    fm = floatmax(T)
    gyro = T === Float16 ? T(16) : T(1024)
    return [
        (z, z, z, p(large), p(larger)),                                  # h overflows in zero field
        (z, (p(small), zero(T), zero(T)), z, p(large), p(larger)),       # h overflows, finite kick
        (z, (p(-small), zero(T), zero(T)), z, p(-large), p(-larger)),    # qm*dt underflows
        ((T(1), T(2), T(3)), z, (zero(T), zero(T), p(strong)), T(2), T(1)),  # |t|² overflows
        ((fm / 4, fm / 8, zero(T)), z, (zero(T), zero(T), gyro), T(2), T(1)),  # v⁻×t overflows
        (
            (fm * T(0.75), zero(T), zero(T)),
            (fm * T(0.75), zero(T), zero(T)),
            (zero(T), zero(T), T(4)),
            T(2),
            T(1),
        ),
        ((T(1), T(2), T(3)), z, (zero(T), zero(T), sqrt(fm) / T(1.25)), T(2), T(1)),  # f < floatmin
        (
            (p(exponent(floatmin(T)) - 3), zero(T), zero(T)),
            z,
            (zero(T), zero(T), one(T)),
            one(T),
            one(T),
        ),
    ]
end

@testset "Boris range witnesses agree with the exact rotation" begin
    for T in (Float16, Float32, Float64), (v, E, B, qm, dt) in boris_cases(T)
        out, parts = H._boris_kick_parts(v..., E..., B..., qm, dt)
        @test !H._boris_kick_in_range(out, parts, v..., qm, dt)
        for order in ((1, 2, 3), (2, 3, 1), (3, 1, 2), (2, 1, 3)), sign in (-one(T), one(T))
            vv, ee, bb = map(x -> ntuple(i -> x[order[i]], 3), (v, E, B))
            exact, scale2 = boris_exact(vv, ee, bb, sign * qm, dt)
            actual = H.boris_kick(vv..., ee..., bb..., sign * qm, dt)
            @test normwise_ok(actual, exact, 4ulp_sqrt(T, scale2))
            # Float64 takes the exact rational evaluation of the same expression.
            T === Float64 && @test actual == map(x -> H._from_exact(T, x), exact)
        end
    end
    for T in (Float16, Float32, Float64)
        z = (zero(T), zero(T), zero(T))
        # A true overflow stays ±Inf for the caller; non-finite input propagates.
        @test H.boris_kick(z..., one(T), zero(T), zero(T), z..., floatmax(T), T(4)) ==
              (T(Inf), zero(T), zero(T))
        @test all(isnan, H.boris_kick(z..., z..., z..., T(NaN), one(T)))
        @test all(isnan, H.boris_kick(T(NaN), zero(T), zero(T), z..., z..., floatmax(T), T(4)))
        # Extreme pure rotations keep the speed and reverse to the start.
        exponents =
            T === Float16 ? (-20, -12, 9, 14) :
            T === Float32 ? (-140, -80, 80, 120) : (-1050, -600, 600, 1000)
        for e in exponents
            v = (T(0.5), -T(0.25), T(0.125))
            B = (ldexp(one(T), e), -ldexp(one(T), e), zero(T))
            moved = H.boris_kick(v..., z..., B..., T(2), one(T))
            @test all(isfinite, moved)
            @test hypot(Float64.(moved)...) ≈ hypot(Float64.(v)...) rtol = 16eps(T)
            back = H.boris_kick(moved..., z..., B..., T(2), -one(T))
            @test all(i -> abs(back[i] - v[i]) <= 16eps(T), 1:3)
        end
    end
    # Other number types keep the plain formula (exact for rationals).
    v, E, B = (1 // 1, 2 // 1, 3 // 1), (2 // 1, -1 // 1, 0 // 1), (0 // 1, 0 // 1, 4 // 1)
    @test H.boris_kick(v..., E..., B..., 2 // 1, 1 // 1) == first(boris_exact(v, E, B, 2, 1))
    big = BigFloat.((1, 2, 3, 0, 0, 0, 0, 0, 1, 1, 1))
    @test H.boris_kick(big...) == old_boris(big...)
end

# ------------------------------------------------------------------ relativistic witnesses

function rel_cases(::Type{T}) where {T}
    p(k) = ldexp(one(T), k)
    z = (zero(T), zero(T), zero(T))
    large = T === Float16 ? 9 : T === Float32 ? 70 : 520
    tiny_c = p(-(exponent(floatmax(T)) ÷ 2 + 6))
    return [
        ((T(0.25), -T(0.125), T(0.5)), z, z, floatmax(T), T(4), one(T)),          # h overflows
        (z, (p(large), zero(T), zero(T)), (zero(T), zero(T), one(T)), T(2), one(T), one(T)),  # γ⁻ overflows
        ((T(0.5), zero(T), zero(T)), z, (zero(T), zero(T), p(large)), T(2), one(T), one(T)),  # |t|² overflows
        (
            (p(exponent(floatmin(T)) - 2), zero(T), zero(T)),
            z,
            (zero(T), zero(T), one(T)),
            one(T),
            one(T),
            one(T),
        ),
        ((tiny_c / 2, zero(T), zero(T)), z, (zero(T), zero(T), one(T)), one(T), one(T), tiny_c),  # c² underflows
    ]
end

@testset "Relativistic Boris range witnesses agree with a 4096-bit reference" begin
    for T in (Float16, Float32, Float64), (v, E, B, qm, dt, c) in rel_cases(T)
        for sign in (-one(T), one(T))
            reference, scale = boris_rel_reference(v, E, B, sign * qm, dt, c)
            actual = H._boris_rel(v..., E..., B..., sign * qm, dt, c)
            @test rel_ok(actual, reference, scale)
        end
    end
    for T in (Float32, Float64)
        # The zero-field push leaves the velocity unchanged.
        v = (T(0.25), -T(0.125), T(0.5))
        z = (zero(T), zero(T), zero(T))
        @test H._boris_rel(v..., z..., z..., floatmax(T), T(4), one(T)) == v
        # Input validation is unchanged.
        @test_throws ArgumentError H._boris_rel(
            one(T),
            zero(T),
            zero(T),
            z...,
            z...,
            one(T),
            one(T),
            one(T),
        )
    end
end

# ------------------------------------------------------------------ Euler, ES and drift witnesses

function euler_cases(::Type{T}) where {T}
    p(k) = ldexp(one(T), k)
    z = (zero(T), zero(T), zero(T))
    fm = floatmax(T)
    E_big = T === Float16 ? p(12) : T === Float32 ? p(100) : p(800)
    return [
        # (u, E, B, h, qm)
        ((T(1), -T(2), T(3)), z, z, T(2), fm),                                  # h*qm overflows
        (z, (E_big, zero(T), zero(T)), z, floatmin(T), p(-10)),                 # h*qm underflows
        ((fm / 4, zero(T), zero(T)), z, (zero(T), zero(T), T(16)), p(-10), one(T)),  # u×B overflows
        (
            (p(exponent(floatmin(T)) - 2), zero(T), zero(T)),
            z,
            (zero(T), zero(T), one(T)),
            one(T),
            one(T),
        ),
    ]
end

@testset "Euler half-kick and drift range witnesses" begin
    for T in (Float16, Float32, Float64), (u, E, B, h, qm) in euler_cases(T)
        vec(x) = ntuple(c -> [x[c]], 3)
        for (mode, sign) in ((:prime, -1), (:predict, 1))
            exact, scale2 = euler_exact(u, E, B, h, qm, sign)
            tol = 4ulp_sqrt(T, scale2)
            if mode === :prime
                v = vec(u)
                H._prime_leapfrog!(v, vec(E), vec(B), qm, h, 1)
                @test normwise_ok(map(only, v), exact, tol)
                T === Float64 && @test map(only, v) == map(x -> H._from_exact(T, x), exact)
                # The planar EMPIC1D back-up is the 3-D expression with B = (0, 0, Bz).
                planar_exact, planar_scale2 = euler_exact(
                    (u[1], u[2], zero(T)),
                    (E[1], E[2], zero(T)),
                    (zero(T), zero(T), B[3]),
                    h,
                    qm,
                    -1,
                )
                w = vec((u[1], u[2], zero(T)))
                H._prime_empic1d!(w, [E[1]], [E[2]], [B[3]], qm, h, 1, false, one(T))
                @test normwise_ok(map(only, w), planar_exact, 4ulp_sqrt(T, planar_scale2))
            else
                predicted = vec((zero(T), zero(T), zero(T)))
                H._predict_half_kick!(predicted, vec(u), vec(E), vec(B), qm, h, 1)
                @test normwise_ok(map(only, predicted), exact, tol)
                # An output that aliases the input takes the single guarded pass.
                aliased = vec(u)
                H._predict_half_kick!(aliased, aliased, vec(E), vec(B), qm, h, 1)
                @test isequal(map(only, aliased), map(only, predicted))
            end
        end
        # Electrostatic priming of each component, v - (h*qm)*E.
        for c = 1:3
            exact = Q(u[c]) - Q(h) * Q(qm) * Q(E[c])
            scale2 = max(Q(u[c])^2, (Q(h) * Q(qm) * Q(E[c]))^2)
            primed = H._es_prime_component(u[c], E[c], h, qm, h * qm)
            @test normwise_ok((primed,), (exact,), 4ulp_sqrt(T, scale2))
        end
    end
    for T in (Float16, Float32, Float64)
        fm = floatmax(T)
        # A backed-up velocity beyond T is an explicit error; a prediction stays ±Inf.
        v = ([-T(0.75) * fm], [zero(T)], [zero(T)])
        E = ([fm / 4], [zero(T)], [zero(T)])
        B = ntuple(_ -> [zero(T)], 3)
        @test_throws ArgumentError H._prime_leapfrog!(map(copy, v), E, B, one(T), T(2), 1)
        @test_throws ArgumentError H._prime_empic1d!(
            map(copy, v),
            E[1],
            E[2],
            B[3],
            one(T),
            T(2),
            1,
            false,
            one(T),
        )
        @test_throws ArgumentError H._es_prime_component(v[1][1], E[1][1], T(2), one(T), T(2))
        predicted = ntuple(_ -> [zero(T)], 3)
        H._predict_half_kick!(predicted, map(a -> -a, v), E, B, one(T), T(2), 1)
        @test only(predicted[1]) == T(Inf)

        # Drift x + dt*v: overflowing products of finite coordinates, true
        # overflow, and finite random bit patterns, each with the full and half step.
        cases = (
            (-fm, T(0.75) * fm, T(2)),
            (T(0.75) * fm, -T(0.75) * fm, T(2)),
            (fm, fm, T(2)),
            (zero(T), fm / 4, nextfloat(zero(T))),
            (floatmin(T), -floatmin(T), one(T)),
        )
        bits = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        rng = Xoshiro(0x647269)
        random = [ntuple(_ -> reinterpret(T, rand(rng, bits)), 3) for _ = 1:20_000]
        failures = 0
        for (x, w, dt) in Iterators.flatten((cases, filter(c -> all(isfinite, c), random)))
            for interval in (dt, dt / 2)
                exact = Q(x) + Q(interval) * Q(w)
                scale = abs(Q(x)) + abs(Q(interval) * Q(w))
                coordinate = H._particle_coordinate(x, w, interval)
                failures += !normwise_ok((coordinate,), (exact,), 4ulp_sqrt(T, scale^2))
            end
        end
        @test failures == 0
    end
end

@testset "Relativistic priming range witnesses" begin
    for T in (Float32, Float64)
        p(k) = ldexp(one(T), k)
        z = (zero(T), zero(T), zero(T))
        large = T === Float32 ? 70 : 520
        cases = (
            # (v, E, B, qm, h, c)
            (
                (T(0.25), zero(T), zero(T)),
                (T(4), zero(T), zero(T)),
                z,
                floatmax(T) / 2,
                floatmin(T),
                one(T),
            ),
            (z, (p(large), zero(T), zero(T)), (zero(T), zero(T), one(T)), one(T), one(T), one(T)),
            (
                (p(exponent(floatmin(T)) - 2), zero(T), zero(T)),
                z,
                (zero(T), zero(T), one(T)),
                one(T),
                one(T),
                one(T),
            ),
        )
        for (v, E, B, qm, h, c) in cases, sign in (-one(T), one(T))
            reference, scale = rel_prime_reference(v, E, B, sign * qm, h, c)
            w = ntuple(i -> [v[i]], 3)
            H._prime_leapfrog!(
                w,
                ntuple(i -> [E[i]], 3),
                ntuple(i -> [B[i]], 3),
                sign * qm,
                h,
                1,
                true,
                c,
            )
            @test rel_ok(map(only, w), reference, scale)
            planar_reference, planar_scale = rel_prime_reference(
                v,
                (E[1], E[2], zero(T)),
                (zero(T), zero(T), B[3]),
                sign * qm,
                h,
                c,
            )
            w = ntuple(i -> [v[i]], 3)
            H._prime_empic1d!(w, [E[1]], [E[2]], [B[3]], sign * qm, h, 1, true, c)
            @test rel_ok(map(only, w), planar_reference, planar_scale)
        end
    end
end

# ------------------------------------------------------------------ dt validation

state_bytes(x...) = (io = IOBuffer(); serialize(io, x); take!(io))
solver_state(st) =
    (; (name => getfield(st, name) for name in fieldnames(typeof(st)) if name != :g)...)
unresolved(::Type{T}) where {T} = (nextfloat(zero(T)), floatmin(T), prevfloat(2floatmin(T)))

function particles(::Type{T}, D, n; q = one(T), m = one(T), L = T(2π)) where {T}
    ps = ParticleSet{D,T}(n; q, m)
    for d = 1:D, p = 1:n
        ps.x[d][p] = L * T(p - 1 // 2) / n
    end
    for c = 1:3, p = 1:n
        ps.v[c][p] = T(0.125) * T(mod(3p + c, 5) - 2)
    end
    return ps
end

function hybrid_case(::Type{T}, kind) where {T}
    g = FourierGrid((4, 4), (T(2π), T(2π)))
    ps = particles(T, 2, 16)
    set_density_weight!(ps, one(T), g)
    if kind === :camcl
        st = CAMCLStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), 16)
        init_camcl!(st, ps)
        return st, ps, dt -> step_camcl!(st, ps, dt; NB = 2)
    end
    population = kind === :collection ? [ps] : ps
    st = HybridStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), 16)
    fill!(st.fields.B[3], one(T))
    init!(st, population)
    return st, population, dt -> step!(st, population, dt; NB = 2)
end

function em_case(::Type{T}, State, n_sub) where {T}
    g = FourierGrid((8,), (T(2π),))
    e = particles(T, 1, 8; q = -one(T))
    set_density_weight!(e, one(T), g)
    st = State(g, 8; n0 = one(T), n_sub)
    init_empic!(st, e)
    return st, e, dt -> step_empic!(st, e, dt)
end

function es_case(::Type{T}, kind) where {T}
    if kind === :legacy
        g = FourierGrid((8,), (T(2π),))
        e = particles(T, 1, 8; q = -one(T))
        set_density_weight!(e, one(T), g)
        st = Electrostatic1D(g, 8; n0 = one(T))
        init_espic!(st, e)
        return st, e, dt -> step_espic!(st, e, dt)
    end
    g = FourierGrid((4, 4), (T(2π), T(2π)))
    e = particles(T, 2, 16; q = -one(T))
    set_density_weight!(e, one(T), g)
    if kind === :scalar
        st = ElectrostaticPIC(g, 16; n0 = one(T))
        init_espic!(st, e)
        return st, e, dt -> step_espic!(st, e, dt)
    end
    species = [e]
    st = ElectrostaticPIC(g, species; neutrality = :neutralize)
    init_espic!(st, species)
    return st, species, dt -> step_espic!(st, species, dt)
end

function shock_case(::Type{T}, D, kind = :rk4) where {T}
    L = T(8)
    sh =
        D == 1 ? PerpShock(8, L; Te = T(0.1)) :
        D == 2 ? PerpShock2D(8, 3, L, one(T); Te = T(0.1)) :
        PerpShock3D(8, 3, 3, L, one(T), one(T); Te = T(0.1))
    ps = ParticleSet{D,T}(6)
    for p = 1:6
        ps.x[1][p] = L * T(p) / 7
        for d = 2:D
            ps.x[d][p] = T(p) / 7
        end
        ps.v[1][p] = T(0.25) * T(p - 3)
        ps.v[2][p] = T(0.125)
    end
    fill!(ps.weight, T(0.5))
    (D == 1 ? init_shock! : D == 2 ? init_shock2d! : init_shock3d!)(sh, ps)
    advance = if kind === :leroy
        bc = LeroyBoundary(
            Xoshiro(7);
            V1 = T(0.5),
            vthi = T(0.1),
            V2 = T(0.2),
            vth2 = T(0.1),
            B_down = T(2),
            p_up = T(0.5),
        )
        dt -> step_leroy_shock!(sh, ps, dt; NB = 1, bc)
    elseif kind === :cn
        dt -> step_shock3d!(sh, ps, dt; NB = 1, field_method = :cn)
    else
        dt -> (D == 1 ? step_shock! : D == 2 ? step_shock2d! : step_shock3d!)(sh, ps, dt; NB = 1)
    end
    return sh, ps, advance
end

function hall_case(::Type{T}) where {T}
    g = FourierGrid((4, 4), (T(2π), T(2π)))
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(T(0.1))))
    st.fields.n .= one(T)
    st.fields.ui[1] .= T(0.2)
    return st, nothing, dt -> step_hall_mhd!(st, dt)
end

function transport_case(::Type{T}) where {T}
    g = FourierGrid((4, 4), (T(2π), T(2π)))
    ps = particles(T, 2, 16)
    # One particle per cell: the transport Ohm law needs n >= nfloor everywhere.
    for p = 1:16
        i, j = divrem(p - 1, 4)
        ps.x[1][p] = (i + T(0.5)) * g.dx[1]
        ps.x[2][p] = (j + T(0.5)) * g.dx[2]
    end
    set_density_weight!(ps, one(T), g)
    st = HybridStepper(g, HybridModel(TransportElectrons(T(5 / 3)); η = T(0.02)), CIC(), 16)
    fill!(st.fields.B[3], one(T))
    fill!(st.fields.pe, T(0.4))
    init!(st, ps)
    return st, ps, dt -> step!(st, ps, dt; NB = 2)
end

# A collisional stepper owns a fixed `dt`, which its constructor accepts even
# when unresolved; `step!` must then reject it before mutation.
function collisional_case(::Type{T}, dt) where {T}
    g = FourierGrid((1,), (one(T),))
    ps = ParticleSet{1,T}(2; q = one(T), m = one(T))
    ps.x[1] .= (T(0.25), T(0.75))
    ps.v[1] .= (T(0.5), -T(0.5))
    species = [ps]
    inner = HybridStepper(g, HybridModel(IsothermalElectrons(zero(T))), NGP(), 2)
    st = CollisionalHybridStepper(
        inner,
        species;
        dt,
        rng = Xoshiro(619),
        coulomb_log = one(T),
        epsilon0 = one(T),
        weight_to_number = one(T),
    )
    init!(st, species)
    return st, species, dt -> step!(st, species, dt)
end

function check_dt_contract(name, build, ::Type{T}; smallest = 2floatmin(T)) where {T}
    @testset "$name $T" begin
        st, ps, advance = build()
        before = state_bytes(solver_state(st), ps)
        for dt in unresolved(T)
            @test_throws ArgumentError advance(dt)
            @test state_bytes(solver_state(st), ps) == before
        end
        @test advance(zero(T)) === st
        @test state_bytes(solver_state(st), ps) == before
        @test advance(smallest) === st
        @test st.time[] == smallest
        @test st.step[] == 1
    end
end

@testset "Unresolved timesteps are rejected before mutation" begin
    for T in (Float32, Float64)
        for kind in (:scalar, :collection, :camcl)
            check_dt_contract("hybrid $kind", () -> hybrid_case(T, kind), T)
        end
        for State in (EMPIC1D, EMPIC)
            check_dt_contract("$State", () -> em_case(T, State, 1), T)
            # dt/n_sub is the electron step and must be resolved as well.
            st, e, advance = em_case(T, State, 2)
            before = state_bytes(solver_state(st), e)
            @test_throws ArgumentError advance(2floatmin(T))
            @test state_bytes(solver_state(st), e) == before
            @test advance(4floatmin(T)) === st
            @test st.step[] == 1
        end
        check_dt_contract("Hall-MHD", () -> hall_case(T), T)
        check_dt_contract("hybrid TransportElectrons", () -> transport_case(T), T)
        @testset "collisional hybrid $T" begin
            snapshot(st, species) = state_bytes(
                solver_state(st.inner),
                species,
                st.rng,
                st.source_time,
                st.source_count,
            )
            for dt in unresolved(T)
                st, species, advance = collisional_case(T, dt)
                before = snapshot(st, species)
                # The fixed-dt check would accept this dt; the resolution rule must not.
                err = try
                    advance(dt)
                    nothing
                catch e
                    e
                end
                @test err isa ArgumentError && occursin("2*floatmin", err.msg)
                @test snapshot(st, species) == before
            end
            st, species, advance = collisional_case(T, 2floatmin(T))
            before = snapshot(st, species)
            @test advance(zero(T)) === st
            @test snapshot(st, species) == before
            @test advance(2floatmin(T)) === st
            @test st.inner.time[] == 2floatmin(T)
            @test st.inner.step[] == 1
        end

        for kind in (:legacy, :scalar, :collection)
            @testset "electrostatic $kind $T" begin
                st, ps, advance = es_case(T, kind)
                before = state_bytes(solver_state(st), ps)
                for dt in unresolved(T)
                    @test_throws ArgumentError advance(dt)
                    @test state_bytes(solver_state(st), ps) == before
                end
                @test advance(zero(T)) === st
                @test state_bytes(solver_state(st), ps) == before
                @test advance(2floatmin(T)) === st
                @test st.primed[]
                @test st.last_dt[] == 2floatmin(T)
            end
        end

        for (D, kind) in ((1, :rk4), (1, :leroy), (2, :rk4), (3, :rk4), (3, :cn))
            @testset "shock D$D $kind $T" begin
                sh, ps, advance = shock_case(T, D, kind)
                before = state_bytes(sh, ps)
                for dt in unresolved(T)
                    @test_throws ArgumentError advance(dt)
                    @test state_bytes(sh, ps) == before
                end
                @test advance(2floatmin(T)) === sh
                @test all(a -> all(isfinite, a), (ps.x..., ps.v...))
            end
        end

        @testset "uniform and gathered pushers $T" begin
            for D = 1:3
                ps = particles(T, D, 5)
                E = (T(0.25), -T(0.5), T(0.125))
                B = (T(0.5), T(0.25), -T(1))
                gathered = (ntuple(c -> fill(E[c], 5), 3), ntuple(c -> fill(B[c], 5), 3))
                midpoint = ntuple(_ -> zeros(T, 5), D)
                before = state_bytes(ps, midpoint)
                for dt in (unresolved(T)..., map(-, unresolved(T))...)
                    @test_throws ArgumentError push_uniform!(ps, E, B, dt)
                    @test_throws ArgumentError push_gathered!(ps, gathered..., dt; xmid = midpoint)
                    @test state_bytes(ps, midpoint) == before
                end
                reference = deepcopy(ps)
                @test push_uniform!(ps, E, B, zero(T)) === ps
                @test ps.x == reference.x && ps.v == reference.v
                for dt in (2floatmin(T), -2floatmin(T))
                    @test push_uniform!(ps, E, B, dt) === ps
                    @test push_gathered!(ps, gathered..., dt; xmid = midpoint) === ps
                    @test all(a -> all(isfinite, a), (ps.x..., ps.v..., midpoint...))
                end
            end
        end

        @testset "mpi_step! $T" begin
            ensure_mpi_initialized!()
            layout = LogicalRankLayout((1, 1); periodic = (true, true))
            ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
            try
                st, ps, _ = hybrid_case(T, :scalar)
                before = state_bytes(solver_state(st), ps)
                for dt in unresolved(T)
                    @test_throws ArgumentError mpi_step!(st, ps, ctx, dt; NB = 1)
                    @test state_bytes(solver_state(st), ps) == before
                end
                @test mpi_step!(st, ps, ctx, zero(T); NB = 1) === st
                @test state_bytes(solver_state(st), ps) == before
                @test mpi_step!(st, ps, ctx, 2floatmin(T); NB = 1) === st
                @test st.step[] == 1
            finally
                free_mpi_communicator!(ctx)
            end
        end
    end
end

# ------------------------------------------------------------------ public end-to-end witnesses

@testset "Range repair reaches the public uniform and gathered pushers" begin
    for T in (Float16, Float32, Float64)
        large, larger, small =
            T === Float16 ? (12, 5, -12) : T === Float32 ? (100, 40, -100) : (800, 300, -1000)
        q, dt = ldexp(one(T), large), ldexp(one(T), larger)
        E = (ldexp(one(T), small), zero(T), zero(T))
        z = (zero(T), zero(T), zero(T))
        a = ParticleSet{1,T}(2; q)
        fill!(a.x[1], zero(T))
        b = deepcopy(a)
        midpoint = (zeros(T, 2),)
        push_uniform!(a, E, z, dt)
        push_gathered!(
            b,
            ntuple(c -> fill(E[c], 2), 3),
            ntuple(_ -> zeros(T, 2), 3),
            dt;
            xmid = midpoint,
        )
        kick = Q(q) * Q(dt) * Q(E[1])
        expected = H._from_exact(T, kick)
        @test a.v == b.v && a.x == b.x
        @test all(==(expected), a.v[1]) && all(iszero, a.v[2]) && all(iszero, a.v[3])
        @test all(==(H._from_exact(T, Q(dt) * Q(expected))), a.x[1])
        @test all(==(H._from_exact(T, Q(dt) * Q(expected) / 2)), midpoint[1])

        # A drift whose rounded product overflows lands on the exact coordinate.
        ps = ParticleSet{1,T}(1)
        ps.x[1][1] = -floatmax(T)
        ps.v[1][1] = T(0.75) * floatmax(T)
        push_uniform!(ps, z, z, T(2))
        @test only(ps.x[1]) == H._from_exact(T, Q(-floatmax(T)) + 2Q(T(0.75) * floatmax(T)))
    end
end

physical(es, species) = deepcopy((
    species = [(x = ps.x, v = ps.v, weight = ps.weight, q = ps.q, m = ps.m) for ps in species],
    E = es.E,
    ne = es.ne,
    primed = es.primed[],
    last_dt = es.last_dt[],
))

@testset "Zero-field extreme-charge equilibria stay finite and unchanged" begin
    for T in (Float32, Float64)
        q = ldexp(one(T), T === Float32 ? 100 : 800)
        dt = ldexp(one(T), T === Float32 ? 40 : 300)
        @testset "electrostatic collection $T" begin
            g = FourierGrid((4,), (one(T),))
            species = [ParticleSet{1,T}(1; q = sign * q) for sign in (-1, 1)]
            for ps in species
                ps.x[1][1] = T(0.25)
                ps.weight[1] = inv(q)
            end
            es = ElectrostaticPIC(g, species)
            init_espic!(es, species)
            initial = physical(es, species).species
            for interval in (dt, dt / 2, zero(T), dt)
                step_espic!(es, species, interval)
                @test physical(es, species).species == initial
                @test all(a -> all(iszero, a), es.E)
            end
            # A backed-up velocity beyond the type is an explicit error before publication.
            init_espic!(es, species)
            fill!(es.E[1], floatmax(T))
            for ps in species
                ps.v[1][1] = -sign(ps.q) * floatmax(T) / 2
            end
            before = physical(es, species)
            @test_throws ArgumentError step_espic!(es, species, T(4))
            @test physical(es, species) == before
        end
        for (kind, State, init, advance) in (
            (:hybrid, HybridStepper, init!, step!),
            (:camcl, CAMCLStepper, init_camcl!, step_camcl!),
        )
            @testset "$kind $T" begin
                g = FourierGrid((4, 4), (T(4), T(4)))
                ps = ParticleSet{2,T}(16; q, m = one(T))
                load_lattice!(ps, (zero(T), zero(T)), g.L, g.n)
                set_density_weight!(ps, one(T), g)
                st = State(g, HybridModel(IsothermalElectrons(zero(T))), CIC(), 16)
                init(st, ps)
                x0, v0 = deepcopy(ps.x), deepcopy(ps.v)
                for _ = 1:3
                    @test advance(st, ps, dt; NB = 2) === st
                    @test ps.v == v0 && ps.x == x0
                    @test all(a -> all(isfinite, a), (st.fields.E..., st.fields.B...))
                end
            end
        end
        for (State, D) in ((EMPIC1D, 1), (EMPIC, 2))
            for relativistic in (false, true)
                @testset "$State relativistic=$relativistic $T" begin
                    # Tiny ion mass gives an overflowing ion q/m*dt/2 in a zero field.
                    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
                    e = ParticleSet{D,T}(1; q = -one(T), m = one(T))
                    ions = ParticleSet{D,T}(1; q = one(T), m = floatmin(T))
                    for ps in (e, ions)
                        foreach(a -> fill!(a, T(1.5)), ps.x)
                        fill!(ps.weight, one(T))
                    end
                    st = State(
                        g,
                        1;
                        n0 = zero(T),
                        c = one(T),
                        mobile = true,
                        mi = floatmin(T),
                        relativistic,
                    )
                    init_empic!(st, e, ions)
                    snapshot = deepcopy((e.x, e.v, ions.x, ions.v))
                    for _ = 1:2
                        @test step_empic!(st, e, ions, T(8)) === st
                        @test (e.x, e.v, ions.x, ions.v) == snapshot
                    end
                end
            end
        end
    end
end

@testset "Shock steppers reflect a finite drift whose rounded product overflows" begin
    for T in (Float32, Float64), D = 1:3
        L = T(0.875) * floatmax(T)
        x = T(0.75) * floatmax(T)
        v = -x
        dt = T(2)
        sh =
            D == 1 ? PerpShock(8, L; Te = 0, η = 0, τ = 0, B0 = 0) :
            D == 2 ? PerpShock2D(8, 3, L, one(T); Te = 0, η = 0, τ = 0, B0 = 0) :
            PerpShock3D(8, 3, 3, L, one(T), one(T); Te = 0, η = 0, τ = 0, B0 = 0)
        ps = ParticleSet{D,T}(1)
        foreach(a -> fill!(a, zero(T)), ps.x)
        ps.x[1][1] = x
        ps.v[1][1] = v
        (D == 1 ? init_shock! : D == 2 ? init_shock2d! : init_shock3d!)(sh, ps)
        @test !isfinite(dt * v)
        (D == 1 ? step_shock! : D == 2 ? step_shock2d! : step_shock3d!)(sh, ps, dt; NB = 1)
        @test nparticles(ps) == 1
        @test only(ps.x[1]) == abs(H._from_exact(T, Q(x) + Q(dt) * Q(v)))
        @test only(ps.v[1]) == -v
    end
end

# ------------------------------------------------------------------ batch certificates

# One extreme particle in an ordinary batch forces the guarded loop: the batch
# without it is certified (native loop), its ordinary results are bitwise those
# of the full batch, and the extreme result equals the guarded singleton batch.
pick(a::Tuple, idx) = map(x -> x[idx], a)
others(n, k) = [p for p = 1:n if p != k]
function batch_ok(full, rest, single, n, k)
    o = others(n, k)
    return all(c -> isequal(full[c][o], rest[c]), eachindex(full)) &&
           all(c -> isequal(full[c][k], only(single[c])), eachindex(full)) &&
           all(c -> isfinite(full[c][k]), eachindex(full))
end
boris_certified(v, E, B, qm, dt, n) = first(
    H._boris_batch_certificate(
        qm * dt / 2,
        qm,
        dt,
        H._boris_batch_scales(v, E, B, qm * dt / 2, 1:n)...,
    ),
)
euler_certified(v, E, B, qm, h, n) = H._euler_batch_certificate(
    h,
    qm,
    h * qm,
    H._euler_batch_scales(v, E, B, h * qm, 1:n, typeof(h))...,
)

@testset "Batch owners repair one extreme particle and keep the others native" begin
    for T in (Float32, Float64)
        rng = Xoshiro(0x626174)
        n, k = 129, 57
        fm = floatmax(T)
        o = others(n, k)
        ord(scale = 1) = [T(scale * randn(rng) / 4) for _ = 1:n]
        big = T === Float32 ? ldexp(one(T), 70) : ldexp(one(T), 520)

        @testset "push_uniform! and push_gathered! $T" begin
            # A rotation whose v⁻×t overflows (kick) and a drift whose rounded
            # product overflows (drift) at particle k.
            for extreme in (:kick, :drift)
                ps = ParticleSet{3,T}(n; q = T(2))
                foreach(a -> a .= ord(), ps.v)
                foreach(a -> a .= ord(4), ps.x)
                z = (zero(T), zero(T), zero(T))
                E, B, dt = extreme === :kick ? (z, (zero(T), zero(T), T(4)), one(T)) : (z, z, T(2))
                ps.v[1][k] = T(0.75) * fm
                extreme === :drift && (ps.x[1][k] = -fm)
                @test !boris_certified(ps.v, E, B, ps.q / ps.m, dt, n) ||
                      !H._coordinate_batch_certified(H._coordinate_batch_max(ps.x, T), dt, fm)
                @test boris_certified(pick(ps.v, o), E, B, ps.q / ps.m, dt, n - 1)
                subsets = (
                    ps,
                    ParticleSet{3,T}(
                        pick(ps.x, o),
                        pick(ps.v, o),
                        ps.weight[o],
                        ps.id[o],
                        ps.tag[o],
                        ps.q,
                        ps.m,
                    ),
                    ParticleSet{3,T}(
                        pick(ps.x, [k]),
                        pick(ps.v, [k]),
                        ps.weight[[k]],
                        ps.id[[k]],
                        ps.tag[[k]],
                        ps.q,
                        ps.m,
                    ),
                )
                uniform = map(deepcopy, subsets)
                gathered = map(deepcopy, subsets)
                mids = map(q -> ntuple(_ -> zeros(T, nparticles(q)), 3), gathered)
                for (i, q) in enumerate(uniform)
                    push_uniform!(q, E, B, dt)
                    m = nparticles(gathered[i])
                    push_gathered!(
                        gathered[i],
                        ntuple(c -> fill(E[c], m), 3),
                        ntuple(c -> fill(B[c], m), 3),
                        dt;
                        xmid = mids[i],
                    )
                end
                for set in (uniform, gathered)
                    @test batch_ok(
                        (set[1].v..., set[1].x...),
                        (set[2].v..., set[2].x...),
                        (set[3].v..., set[3].x...),
                        n,
                        k,
                    )
                end
                @test batch_ok(mids[1], mids[2], mids[3], n, k)
                @test map(only, (uniform[3].v...,)) ==
                      H.boris_kick(map(a -> a[k], ps.v)..., E..., B..., ps.q / ps.m, dt)
            end
        end

        @testset "hybrid population push $T" begin
            g = FourierGrid((8,), (T(2π),))
            ps = ParticleSet{1,T}(n)
            ps.x[1] .= T(2π) .* rand(rng, T, n)
            foreach(a -> a .= ord(), ps.v)
            set_density_weight!(ps, one(T), g)
            st = HybridStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), n)
            fill!(st.fields.B[3], T(8))
            init!(st, ps)
            st.step[] = 1
            ps.v[1][k] = T(0.75) * fm    # t = 4: v⁻×t overflows at particle k
            sets = (
                ps,
                ParticleSet{1,T}(
                    pick(ps.x, o),
                    pick(ps.v, o),
                    ps.weight[o],
                    ps.id[o],
                    ps.tag[o],
                    ps.q,
                    ps.m,
                ),
                ParticleSet{1,T}(
                    pick(ps.x, [k]),
                    pick(ps.v, [k]),
                    ps.weight[[k]],
                    ps.id[[k]],
                    ps.tag[[k]],
                    ps.q,
                    ps.m,
                ),
            )
            for q in sets
                ws = H._hybrid_particle_workspace(st, 1:nparticles(q))
                H._push_hybrid_population!(st, ws, q, one(T), one(T))
            end
            @test batch_ok(sets[1].v, sets[2].v, sets[3].v, n, k)
        end

        @testset "EMPIC velocity pushes and drifts $T" begin
            for rel in (false, true), D in (1, 2)
                g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> T(2π), D))
                es =
                    D == 1 ? EMPIC1D(g, n; c = one(T), relativistic = rel) :
                    EMPIC(g, n; c = one(T), relativistic = rel)
                vel = ntuple(_ -> ord(rel ? 0.4 : 1), 3)
                Ep = ntuple(_ -> ord(), 3)
                Bp = ntuple(_ -> ord(), 3)
                if rel
                    Ep[1][k] = big            # |u⁻|²/c² overflows
                else
                    vel[1][k], Bp[3][k] = T(0.75) * fm, T(8)   # v⁻×t overflows
                end
                run(idx) = begin
                    s = ParticleSet{D,T}(length(idx); q = -one(T))
                    foreach(c -> s.v[c] .= vel[c][idx], 1:3)
                    if D == 1
                        H._push_velocities!(es, s, Ep[1][idx], Ep[2][idx], Bp[3][idx], T(0.5))
                    else
                        H._push_velocities!(es, s, pick(Ep, idx), pick(Bp, idx), T(0.5))
                    end
                    s.v
                end
                @test batch_ok(run(1:n), run(o), run([k]), n, k)
            end
            for D in (1, 2)
                xs = ntuple(_ -> ord(4), D)
                vs = ntuple(_ -> ord(), 3)
                xs[1][k], vs[1][k] = -fm, T(0.75) * fm   # x + dt v overflows, exact x is finite
                L = ntuple(_ -> T(2π), D)
                drift(idx) = begin
                    s = ParticleSet{D,T}(length(idx))
                    foreach(d -> s.x[d] .= xs[d][idx], 1:D)
                    foreach(c -> s.v[c] .= vs[c][idx], 1:3)
                    mid = ntuple(_ -> zeros(T, length(idx)), D)
                    D == 1 ? H._drift_record_mid!(s, mid[1], T(2), L[1]) :
                    H._drift_record_mid!(s, mid, T(2), L)
                    (s.x..., mid...)
                end
                @test batch_ok(drift(1:n), drift(o), drift([k]), n, k)
            end
        end

        @testset "electrostatic priming, kick and drift $T" begin
            g = FourierGrid((8,), (T(2π),))
            xs, vs, Es = T(2π) .* rand(rng, T, n), ord(), ord(0.1)
            # A subnormal velocity and field at particle k: S below the threshold.
            vs[k], Es[k] =
                ldexp(one(T), exponent(floatmin(T)) - 3), ldexp(one(T), exponent(floatmin(T)) - 2)
            push(idx) = begin
                s = ParticleSet{1,T}(length(idx); q = -one(T))
                s.x[1] .= xs[idx]
                s.v[1] .= vs[idx]
                H._espic_push!(s, (Es[idx],), g, T(0.25), false, zero(T))
                s.v
            end
            @test batch_ok(push(1:n), push(o), push([k]), n, k)
        end

        @testset "priming and predictor batches $T" begin
            for (kind, rel) in (
                (:leapfrog, false),
                (:leapfrog, true),
                (:empic1d, false),
                (:empic1d, true),
                (:predict, false),
            )
                u = ntuple(_ -> ord(rel ? 0.4 : 1), 3)
                Ep = ntuple(_ -> ord(), 3)
                Bp = ntuple(_ -> ord(), 3)
                if rel
                    Ep[1][k] = big            # |u⁻|²/c² overflows
                else
                    u[1][k] = fm / 4          # u×B overflows
                    Bp[3][k] = T(16)
                    Bp[2][k] = zero(T)
                end
                h, qm = T(2)^-10, one(T)
                rel || @test !euler_certified(u, Ep, Bp, qm, h, n)
                prime(idx) = begin
                    w = map(a -> a[idx], u)
                    m = length(idx)
                    if kind === :leapfrog
                        H._prime_leapfrog!(
                            w,
                            pick(Ep, idx),
                            pick(Bp, idx),
                            qm,
                            rel ? one(T) : h,
                            m,
                            rel,
                            one(T),
                        )
                    elseif kind === :empic1d
                        H._prime_empic1d!(
                            w,
                            Ep[1][idx],
                            Ep[2][idx],
                            Bp[3][idx],
                            qm,
                            rel ? one(T) : h,
                            m,
                            rel,
                            one(T),
                        )
                    else
                        out = map(similar, w)
                        H._predict_half_kick!(out, w, pick(Ep, idx), pick(Bp, idx), qm, h, m)
                        w = out
                    end
                    w
                end
                @test batch_ok(prime(1:n), prime(o), prime([k]), n, k)
            end
        end

        @testset "shock push batches $T" begin
            # A particle at rest has the unknown scale |h E_p|, so it forces the
            # guarded loop; the ordinary particles match the certified batch.
            for D = 1:3
                L = T(8)
                sh =
                    D == 1 ? PerpShock(8, L; Te = T(0.1)) :
                    D == 2 ? PerpShock2D(8, 3, L, one(T); Te = T(0.1)) :
                    PerpShock3D(8, 3, 3, L, one(T), one(T); Te = T(0.1))
                m = 12
                ps = ParticleSet{D,T}(m)
                for p = 1:m
                    ps.x[1][p] = L * T(p) / (m + 1)
                    for d = 2:D
                        ps.x[d][p] = T(p) / (m + 1)
                    end
                    foreach(c -> ps.v[c][p] = T(randn(rng) / 8), 1:3)
                end
                fill!(ps.weight, T(0.5))
                (D == 1 ? init_shock! : D == 2 ? init_shock2d! : init_shock3d!)(sh, ps)
                foreach(c -> ps.v[c][5] = zero(T), 1:3)
                rest = [p for p = 1:m if p != 5]
                Egrid, Bgrid = D == 3 ? (sh.E, sh.B) : ((sh.Ex, sh.Ey), (sh.Bz,))
                @test !H._shock_batch_certified(ps, Egrid, Bgrid, T(0.05))
                ps2 = ParticleSet{D,T}(
                    pick(ps.x, rest),
                    pick(ps.v, rest),
                    ps.weight[rest],
                    ps.id[rest],
                    ps.tag[rest],
                    ps.q,
                    ps.m,
                )
                @test H._shock_batch_certified(ps2, Egrid, Bgrid, T(0.05))
                sh2 = deepcopy(sh)
                step = D == 1 ? step_shock! : D == 2 ? step_shock2d! : step_shock3d!
                step(sh, ps, T(0.05); NB = 1)
                step(sh2, ps2, T(0.05); NB = 1)
                @test nparticles(ps) == m && nparticles(ps2) == m - 1
                @test all(c -> isequal(ps.v[c][rest], ps2.v[c]), 1:3)
                @test all(d -> isequal(ps.x[d][rest], ps2.x[d]), 1:D)
            end
        end
    end
end

include("particle_kinematics/lorentz_gap.jl")
include("particle_kinematics/priming_publication.jl")
end
