using Test, HybridPlasmaPIC, Serialization

const _C1H = HybridPlasmaPIC
const _C1Q = Rational{BigInt}

function _c1_physical_bytes(es, species...)
    # The undo workspace is scratch, not physical/checkpoint state. Include
    # every other solver field, the complete grid (all four FFT buffers), and
    # complete particles rather than reusing the implementation's inventory.
    names = filter(!=(:transaction), fieldnames(typeof(es)))
    state = NamedTuple{names}(map(name -> getfield(es, name), names))
    io = IOBuffer()
    serialize(io, (state, species))
    return take!(io)
end

function _c1_set_transverse!(es, value)
    es isa EMPIC1D ? fill!(es.Ey, value) : fill!(es.E[2], value)
    return nothing
end

function _c1_uniform(State, ::Type{T}; mobile = false, mi = one(T), n_sub = 1) where {T}
    g = FourierGrid((4,), (T(4),))
    es = State(g, 4; relativistic = true, c = one(T), mobile, mi, n_sub)
    electrons = ParticleSet{1,T}(4; q = -one(T))
    electrons.x[1] .= T.(0:3)
    ions = mobile ? ParticleSet{1,T}(4; m = mi) : nothing
    ions === nothing || (ions.x[1] .= T.(0:3))
    init_empic!(es, electrons, ions)
    fields = es isa EMPIC1D ? (es.Ex, es.Ey, es.Bz) : (es.E..., es.B...)
    @test all(a -> all(iszero, a), fields)
    return es, electrons, ions
end

@testset "Relativistic rounded outputs reject the unsupported velocity boundary" begin
    for T in (Float16, Float32, Float64)
        z = ntuple(_ -> zero(T), 3)
        c, dt = one(T), T(2)
        E = ldexp(one(T), cld(precision(T), 2))
        electric = (E, zero(T), zero(T))
        # Independent exact rounding-cell proof: the physical axial result is
        # strictly subluminal but its nearest target component is exactly c.
        u = 2 * _C1Q(E)
        midpoint = (_C1Q(prevfloat(c)) + _C1Q(c)) / 2
        @test midpoint^2 < u^2 / (1 + u^2) < 1
        @test_throws ArgumentError _C1H._boris_rel(z..., electric..., z..., one(T), dt, c)
        @test_throws ArgumentError _C1H._boris_rel_native(z..., electric..., z..., one(T), dt, c)
        # Raw cold arithmetic keeps its numerical role; publication owns the
        # rejection and does not move this correctly rounded result inward.
        @test _C1H._boris_rel_repair(z, z..., electric..., z..., one(T), dt, c) ==
              (c, zero(T), zero(T))
        ep, bp = ntuple(i -> [electric[i]], 3), ntuple(_ -> T[0], 3)
        for planar in (false, true)
            v = ntuple(_ -> T[0], 3)
            if planar
                @test_throws ArgumentError _C1H._prime_empic1d!(
                    v,
                    ep[1],
                    ep[2],
                    bp[3],
                    one(T),
                    one(T),
                    1,
                    true,
                    c,
                )
            else
                @test_throws ArgumentError _C1H._prime_leapfrog!(
                    v,
                    ep,
                    bp,
                    one(T),
                    one(T),
                    1,
                    true,
                    c,
                )
            end
        end
        if T !== Float16
            for State in (EMPIC1D, EMPIC)
                es = State(FourierGrid((4,), (T(4),)), 1; relativistic = true, c)
                ps = ParticleSet{1,T}(1)
                scales = _C1H._boris_batch_scales(ps.v, ep, bp, one(T), 1:1)
                # Branch coverage only: the expected rejection above comes
                # from the independent rounding-cell inequality.
                @test _C1H._boris_rel_batch_certificate(one(T), one(T), dt, c, scales...)
                if State === EMPIC1D
                    @test_throws ArgumentError _C1H._push_velocities!(
                        es,
                        ps,
                        ep[1],
                        ep[2],
                        bp[3],
                        dt,
                    )
                else
                    @test_throws ArgumentError _C1H._push_velocities!(es, ps, ep, bp, dt)
                end
            end
        end
    end
end

@testset "EMPIC rollback restores old gather lengths after population changes" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC), established in (false, true)
        g = FourierGrid((4,), (T(4),))
        es = State(g, 0; relativistic = true, c = one(T))
        empty = ParticleSet{1,T}(0; q = -one(T))
        init_empic!(es, empty)
        if established
            step_empic!(es, empty, T(0.125))
        end
        ps = ParticleSet{1,T}(4; q = -one(T))
        ps.x[1] .= T.(0:3)
        es.n_sub = 3
        strong = ldexp(one(T), cld(precision(T), 2) + 2)
        _c1_set_transverse!(es, strong)
        before = _c1_physical_bytes(es, ps)
        @test_throws ArgumentError step_empic!(es, ps, T(2))
        @test _c1_physical_bytes(es, ps) == before
        gathers =
            es isa EMPIC1D ? (es.Exp, es.Eyp, es.Bzp, es.work) :
            (es.Ep..., es.Bp..., es.mide..., es.worke)
        @test all(isempty, gathers)
        _c1_set_transverse!(es, T(0.001))
        @test step_empic!(es, ps, T(0.001)) === es
        @test all(a -> length(a) == 4, gathers)
    end
end

@testset "EMPIC rollback preserves borrowed particle arrays and ignored ions" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        backing = ParticleSet{1,T}(8; q = -one(T))
        fill!(backing.x[1], T(-13))
        index = 1:2:8
        ps = ParticleSet{1,T}(
            map(a -> view(a, index), backing.x),
            map(a -> view(a, index), backing.v),
            view(backing.weight, index),
            view(backing.id, index),
            view(backing.tag, index),
            -one(T),
            one(T),
        )
        ps.x[1] .= T.(0:3)
        es = State(FourierGrid((4,), (T(4),)), 4; relativistic = true, c = one(T))
        init_empic!(es, ps)
        ignored = ParticleSet{1,T}(3; q = T(NaN), m = T(Inf))
        fill!(ignored.x[1], T(-7))
        fill!(ignored.v[1], one(T))
        _c1_set_transverse!(es, ldexp(one(T), cld(precision(T), 2)))
        before = _c1_physical_bytes(es, ps, backing, ignored)
        @test_throws ArgumentError step_empic!(es, ps, ignored, T(2))
        @test _c1_physical_bytes(es, ps, backing, ignored) == before
        @test all(==(T(-13)), backing.x[1][2:2:8])
        _c1_set_transverse!(es, T(0.001))
        # Compare the ignored species by itself because the solver advances.
        io = IOBuffer()
        serialize(io, ignored)
        saved_ignored = take!(io)
        @test step_empic!(es, ps, ignored, T(0.001)) === es
        io = IOBuffer()
        serialize(io, ignored)
        @test take!(io) == saved_ignored
    end
end

# A deterministic throwing shape exercises the existing local borrowing
# boundary without mutable callback state. Initial/gather coordinates are
# outside the rejecting interval; the first drift midpoint is inside it.
struct _C1MidpointFailureShape <: _C1H.ShapeFunction end
_C1H.width(::_C1MidpointFailureShape) = 1
function _C1H._stencil1d(::_C1MidpointFailureShape, s)
    oftype(s, 0.5) < s < one(s) && error("C1 midpoint stencil rejection")
    return _C1H._stencil1d(NGP(), s)
end

@testset "EMPIC full undo includes local midpoint cleanup and diagnostic scratch" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        es = State(
            FourierGrid((4,), (T(4),)),
            1;
            relativistic = true,
            c = T(10),
            shape = _C1MidpointFailureShape(),
        )
        ps = ParticleSet{1,T}(1; q = -one(T))
        ps.x[1][1] = T(0.25)
        ps.v[1][1] = one(T)
        init_empic!(es, ps)
        for (i, name) in enumerate((:cbuf, :tbuf, :abuf, :sbuf))
            fill!(getproperty(es.g, name), Complex{T}(T(i), -T(i)))
        end
        before = _c1_physical_bytes(es, ps)
        @test_throws "C1 midpoint stencil rejection" step_empic!(es, ps, one(T))
        @test _c1_physical_bytes(es, ps) == before
    end
end

@testset "EMPIC restores complete public state after output-range rejection" begin
    for T in (Float32, Float64),
        State in (EMPIC1D, EMPIC),
        firststep in (true, false),
        mobile in (false, true),
        ns in (1, 3)
        # For mobile runs, a light ion fails after the electron has completed
        # its priming/kick. Immobile runs exercise the native electron bypass.
        exponent = cld(precision(T), 2)
        mi = mobile ? ldexp(one(T), -exponent - 2) : one(T)
        es, electrons, ions = _c1_uniform(State, T; mobile, mi, n_sub = ns)
        if !firststep
            step_empic!(es, electrons, ions, T(2))
            @test all(a -> all(iszero, a), electrons.v)
        end
        _c1_set_transverse!(es, mobile ? one(T) : ldexp(one(T), exponent))
        before = _c1_physical_bytes(es, electrons, ions)
        identities = map(
            objectid,
            (electrons.x..., electrons.v..., electrons.weight, electrons.id, electrons.tag),
        )
        @test_throws ArgumentError step_empic!(es, electrons, ions, T(2))
        @test _c1_physical_bytes(es, electrons, ions) == before
        @test map(
            objectid,
            (electrons.x..., electrons.v..., electrons.weight, electrons.id, electrons.tag),
        ) == identities
        @test all(a -> all(iszero, a), electrons.v)
        @test ions === nothing || all(a -> all(iszero, a), ions.v)
        # Rejection leaves a valid writable/restartable state, including the
        # exact pre-call half-level history. Retry with a supported field.
        mktempdir() do dir
            path = joinpath(dir, "rejected-state.bin")
            save_checkpoint(path, es, electrons, ions)
            @test load_checkpoint!(es, electrons, ions, path) === es
            fresh = State(
                FourierGrid((4,), (T(4),)),
                0;
                relativistic = true,
                c = one(T),
                mobile,
                mi,
                n_sub = ns,
            )
            fresh_e = ParticleSet{1,T}(0; q = -one(T))
            fresh_i = mobile ? ParticleSet{1,T}(0; m = mi) : nothing
            @test fresh.transaction === nothing
            @test load_checkpoint!(fresh, fresh_e, fresh_i, path) === fresh
            @test fresh.transaction === nothing
            loaded = _c1_physical_bytes(fresh, fresh_e, fresh_i)
            @test_throws ArgumentError step_empic!(fresh, fresh_e, fresh_i, T(2))
            @test _c1_physical_bytes(fresh, fresh_e, fresh_i) == loaded
            _c1_set_transverse!(fresh, T(0.01) * mi)
            @test step_empic!(fresh, fresh_e, fresh_i, T(0.001)) === fresh
            archive = joinpath(dir, "rejected-state.archive")
            archive_run(archive, es, electrons, ions; rng_seed = 1)
            @test load_archive(archive).state.step == es.step[]
        end
        _c1_set_transverse!(es, T(0.01) * mi)
        @test step_empic!(es, electrons, ions, T(0.001)) === es
        @test isfinite(kinetic_energy_relativistic(electrons, one(T)))
    end
end

@testset "Relativistic initialization rejects before state or buffer changes" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC), mobile in (false, true)
        es, electrons, ions = _c1_uniform(State, T; mobile)
        species = mobile ? ions : electrons
        species.v[1][end] = one(T)
        before = _c1_physical_bytes(es, electrons, ions)
        workspace = es.transaction
        @test_throws ArgumentError init_empic!(es, electrons, ions)
        @test _c1_physical_bytes(es, electrons, ions) == before
        @test es.transaction === workspace
        @test step_empic!(es, electrons, ions, zero(T)) === es
        @test _c1_physical_bytes(es, electrons, ions) == before
    end
end

@testset "General EMPIC output rejection restores 2D and 3D state" begin
    for T in (Float32, Float64), D in (2, 3), mobile in (false, true), firststep in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        n = prod(g.n)
        k = cld(precision(T), 2)
        mi = mobile ? ldexp(one(T), -k - 2) : one(T)
        es = EMPIC(g, n; relativistic = true, c = one(T), mobile, mi)
        e = ParticleSet{D,T}(n; q = -one(T))
        ions = mobile ? ParticleSet{D,T}(n; m = mi) : nothing
        for ps in (e, ions)
            ps === nothing && continue
            for (p, I) in enumerate(CartesianIndices(g.n)), d = 1:D
                ps.x[d][p] = T(I[d] - 1)
            end
        end
        init_empic!(es, e, ions)
        @test all(a -> all(iszero, a), (es.E..., es.B...))
        firststep || step_empic!(es, e, ions, T(2))
        _c1_set_transverse!(es, mobile ? one(T) : ldexp(one(T), k))
        before = _c1_physical_bytes(es, e, ions)
        @test_throws ArgumentError step_empic!(es, e, ions, T(2))
        @test _c1_physical_bytes(es, e, ions) == before
    end
end

# Frozen fd0a168 generic RR arithmetic, used only to check valid ambient bits.
# Domain correctness is classified separately from the original components.
function _c1_old_generic_rr_velocity(v, E, B, dt, K, c)
    T = typeof(c)
    c2 = c * c
    ux, uy, uz = v
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    v2 = ux * ux + uy * uy + uz * uz
    gamma = one(T) / sqrt(one(T) - v2 / c2)
    wx = Ex + (uy * Bz - uz * By)
    wy = Ey + (uz * Bx - ux * Bz)
    wz = Ez + (ux * By - uy * Bx)
    vE = ux * Ex + uy * Ey + uz * Ez
    chi2 = (wx * wx + wy * wy + wz * wz) - (vE * vE) / c2
    chi2 > 0 || return v
    s = gamma * exp(-K * gamma * chi2 * dt)
    px, py, pz = s * ux, s * uy, s * uz
    gn = sqrt(one(T) + (px * px + py * py + pz * pz) / c2)
    return (px / gn, py / gn, pz / gn)
end

@testset "Generic radiation preserves ambient bits and stages range rejection" begin
    stored = setprecision(BigFloat, 256) do
        (BigFloat(1) / 4, -BigFloat(1) / 8, BigFloat(1) / 16)
    end
    for bits in (11, 24, 53, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                c, dt, K = one(BigFloat), BigFloat(0.125), BigFloat(0.25)
                E, B = BigFloat.((0, 1, 0)), BigFloat.((0, 0, 2))
                reference = _c1_old_generic_rr_velocity(stored, E, B, dt, K, c)
                @test sum(abs2, _C1Q.(reference)) < _C1Q(c)^2
                for D = 1:3
                    ps = ParticleSet{D,BigFloat}(1)
                    foreach(a -> fill!(a, zero(BigFloat)), ps.x)
                    for d = 1:3
                        ps.v[d][1] = stored[d]
                    end
                    apply_radiation_reaction!(ps, E, B, dt; K, c)
                    @test isequal(map(only, ps.v), reference)
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                end
            end
        end
    end
    # At 11-bit upward rounding the native normalization returns an exact-
    # interior diagonal vector whose rounded beta² is one. It is outside the
    # unchanged public D_narrow even though its exact norm is below c. Freeze
    # the independently reproduced parent arithmetic rather than widening it.
    setprecision(BigFloat, 11) do
        setrounding(BigFloat, RoundUp) do
            c = one(BigFloat)
            z = zero(BigFloat)
            v = ntuple(_ -> BigFloat(1181) / 2048, 3)
            E, B = (z, c, z), (z, z, z)
            K, dt = ldexp(c, -44), c
            b = v ./ c
            @test b[1] * b[1] + b[2] * b[2] + b[3] * b[3] < 1
            @test sum(abs2, _C1Q.(v)) < _C1Q(c)^2
            native = _c1_old_generic_rr_velocity(v, E, B, dt, K, c)
            b = native ./ c
            @test b[1] * b[1] + b[2] * b[2] + b[3] * b[3] == 1
            @test sum(abs2, _C1Q.(native)) < _C1Q(c)^2
            ps = ParticleSet{1,BigFloat}(2)
            fill!(ps.x[1], z)
            ps.v[1][1] = c / 4
            for d = 1:3
                ps.v[d][2] = v[d]
            end
            io = IOBuffer()
            serialize(io, ps)
            before = take!(io)
            @test_throws ArgumentError apply_radiation_reaction!(ps, E, B, dt; K, c)
            io = IOBuffer()
            serialize(io, ps)
            @test take!(io) == before
            @test precision(BigFloat) == 11
            @test rounding(BigFloat) == RoundUp
        end
    end
end

# Independent D_narrow classifier for the generic public surface. In this
# near-boundary corpus the exact comparison is always required; it does not
# query the production predicate, flags or candidate output differences.
function _c1_independent_domain(v, c)
    all(isfinite, v) || return false
    b = v ./ c
    beta2 = b[1] * b[1] + b[2] * b[2] + b[3] * b[3]
    return isfinite(beta2) && beta2 < 1 && sum(abs2, _C1Q.(v)) < _C1Q(c)^2
end

@testset "Generic radiation scale, vector and stored-precision admission" begin
    valid = rejected_input = rejected_output = native_error = 0
    for bits in (11, 24, 53),
        mode in (RoundNearest, RoundDown, RoundUp, RoundToZero),
        stored_bits in (11, 256),
        scale in (3 // 4, 3 // 2, 11 // 10),
        direction in ((1, 0, 0), (1, 1, 1), (1, 2, 3)),
        margin in (2, 8)

        stored = setprecision(BigFloat, stored_bits) do
            c = BigFloat(scale)
            direction = BigFloat.(direction)
            factor =
                (1 - BigFloat(margin) * BigFloat(2)^(-bits + 1)) * c / sqrt(sum(abs2, direction))
            map(x -> x * factor, direction)
        end
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                # Public parameter conversion is mirrored before the frozen
                # native formula. Velocity storage deliberately retains its
                # independently selected precision.
                c = BigFloat(scale)
                z = zero(BigFloat)
                E, B = (z, one(BigFloat), z), (z, z, z)
                dt, K = one(BigFloat), ldexp(one(BigFloat), -4bits)
                ps = ParticleSet{1,BigFloat}(2)
                fill!(ps.x[1], z)
                ps.v[1][1] = c / 4
                for d = 1:3
                    ps.v[d][2] = stored[d]
                end
                io = IOBuffer()
                serialize(io, ps)
                before = take!(io)
                if !_c1_independent_domain(stored, c)
                    rejected_input += 1
                    @test_throws ArgumentError apply_radiation_reaction!(ps, E, B, dt; K, c)
                    io = IOBuffer()
                    serialize(io, ps)
                    @test take!(io) == before
                else
                    reference = try
                        _c1_old_generic_rr_velocity(stored, E, B, dt, K, c)
                    catch err
                        err isa DomainError || rethrow()
                        err
                    end
                    if reference isa DomainError
                        native_error += 1
                        @test_throws DomainError apply_radiation_reaction!(ps, E, B, dt; K, c)
                        io = IOBuffer()
                        serialize(io, ps)
                        @test take!(io) == before
                    elseif _c1_independent_domain(reference, c)
                        valid += 1
                        first = _c1_old_generic_rr_velocity((c / 4, z, z), E, B, dt, K, c)
                        apply_radiation_reaction!(ps, E, B, dt; K, c)
                        @test isequal(map(a -> a[2], ps.v), reference)
                        @test isequal(map(a -> a[1], ps.v), first)
                    else
                        rejected_output += 1
                        @test_throws ArgumentError apply_radiation_reaction!(ps, E, B, dt; K, c)
                        io = IOBuffer()
                        serialize(io, ps)
                        @test take!(io) == before
                    end
                end
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end
    @test valid + rejected_input + rejected_output + native_error == 432
    @info "C1 generic domain coverage" valid rejected_input rejected_output native_error
end
