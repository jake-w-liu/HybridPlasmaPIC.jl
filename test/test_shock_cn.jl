module ShockCNSolveTests
using HybridPlasmaPIC, LinearAlgebra, Test
const H = HybridPlasmaPIC

function periodic_matrix(n, length, ::Type{T}, ::Type{R} = Float64) where {T,R}
    # Explicit Fourier-series differentiation, independent of production FFTs.
    D = zeros(R, n, n)
    for i = 0:n-1, j = 0:n-1, mode = 0:n-1
        signed_mode = mode <= n ÷ 2 ? mode : mode - n
        k = iseven(n) && mode == n ÷ 2 ? 0.0 : R(T(2π) * signed_mode / length)
        D[i+1, j+1] += real(im * k * cis(R(2) * R(π) * mode * (i - j) / n)) / n
    end
    return D
end

function independent_operator(sh; R = Float64)
    nx, ny, nz = sh.nx, sh.ny, sh.nz
    N = nx * ny * nz
    T = eltype(sh.B[1])
    dx = R(sh.sbp.dx)
    dx1 = zeros(R, nx, nx)
    dx1[1, 1] = -1 / dx
    dx1[1, 2] = 1 / dx
    dx1[nx, nx-1] = -1 / dx
    dx1[nx, nx] = 1 / dx
    for i = 2:nx-1
        dx1[i, i-1] = -1 / (2dx)
        dx1[i, i+1] = 1 / (2dx)
    end
    eye(n) = Matrix{R}(I, n, n)
    Dx = kron(eye(nz), eye(ny), dx1)
    Dy = kron(eye(nz), periodic_matrix(ny, sh.Ly, T, R), eye(nx))
    Dz = kron(periodic_matrix(nz, sh.Lz, T, R), eye(ny), eye(nx))
    Z = zeros(R, N, N)
    C = [Z -Dz Dy; Dz Z -Dx; -Dy Dx Z]
    ux, uy, uz = map(u -> Diagonal(R.(vec(u))), sh.u)
    crossu = [Z -uz uy; uz Z -ux; -uy ux Z]
    A = C * crossu - R(sh.η) * C * C
    b = zeros(R, 3N)
    damping = R(sh.τ) / R(sh.sbp.H[nx])
    for k = 1:nz, j = 1:ny
        index = nx + (j - 1) * nx + (k - 1) * nx * ny
        for c = 0:2
            A[c*N+index, c*N+index] -= damping
        end
        b[2N+index] = damping * R(sh.B0)
    end
    return A, b
end

flatten(fields) = vcat(map(x -> Float64.(vec(x)), fields)...)

function setup(::Type{T}, kind) where {T}
    eta = kind == :diffusion ? T(1) : kind == :mixed ? T(0.1) : zero(T)
    tau = kind == :sat ? T(3) : kind == :mixed ? T(0.2) : zero(T)
    sh = PerpShock3D(3, 5, 3, T(2), T(2π), T(2π); η = eta, τ = tau, B0 = T(0.5), Te = zero(T))
    for k = 1:sh.nz, j = 1:sh.ny, i = 1:sh.nx
        sh.B[3][i, j, k] = T(0.5) + T(0.125) * cos(T(2π) * T(j - 1) / T(sh.ny))
        if kind == :advection
            sh.u[2][i, j, k] = one(T)
        elseif kind == :mixed
            sh.u[1][i, j, k] = T(0.1) * T(i)
            sh.u[2][i, j, k] = T(0.25) * sin(T(j))
            sh.u[3][i, j, k] = T(0.1) * cos(T(k))
            sh.B[1][i, j, k] = T(0.1) * sin(T(j + k))
            sh.B[2][i, j, k] = T(0.1) * cos(T(i + k))
        end
    end
    return sh
end

function compare_step(sh, dt; public = false, reference_type = Float64)
    T = eltype(sh.B[1])
    before = flatten(sh.B)
    A, b = independent_operator(sh; R = reference_type)
    K = map(similar, sh.B)
    H._b_rhs3d!(K, sh.B, sh)
    rhs_error = norm(flatten(K) - (A * before + b), Inf)
    rhs_scale = max(norm(A * before + b, Inf), 1.0)
    @test rhs_error <= 300eps(T) * rhs_scale
    h = Float64(T(dt))
    lhs = I - (h / 2) * A
    rhs = (I + (h / 2) * A) * before + h * b
    reference = lhs \ rhs
    reference_residual = norm(lhs * reference - rhs, Inf) / max(norm(rhs, Inf), 1.0)
    @test reference_residual < 2e-12
    if public
        ps = ParticleSet{3,T}(0)
        init_shock3d!(sh, ps)
        @test all(iszero, flatten(sh.u))
        step_shock3d!(sh, ps, T(dt); NB = 1, field_method = :cn)
    else
        H._cn_b3d!(sh, T(dt))
    end
    observed = flatten(sh.B)
    residual = norm(lhs * observed - rhs, Inf) / max(norm(rhs, Inf), 1.0)
    relative_error = norm(observed - reference, Inf) / max(norm(reference, Inf), 1.0)
    # The candidate must solve the declared trapezoidal equation, rather than
    # accepting the fixed-point iterate solely because the loop has ended.
    @test residual <= 500eps(T)
    @test relative_error <= 1000eps(T)
    return (;
        rhs_error,
        reference_residual,
        residual,
        relative_error,
        observed_max = maximum(abs, observed),
        reference_max = maximum(abs, reference),
    )
end


@testset "shock CN independent trapezoid" begin
    for T in (Float32, Float64),
        kind in (:advection, :diffusion, :sat, :mixed),
        dt in (0.25, 1.0, 4.0)

        compare_step(setup(T, kind), dt)
    end
    for T in (Float32, Float64), dt in (4.0, 10000.0)
        # The full SBP operator at h=10000 has condition number ~3.8e8.
        # Independently assemble at high precision to resolve the reference's
        # nullspace; Float64 dense assembly failed its analytic-mode control.
        compare_step(
            setup(T, :diffusion),
            dt;
            public = true,
            reference_type = dt == 10000 ? BigFloat : Float64,
        )
    end
    for T in (Float32, Float64), dims in ((3, 1, 1), (4, 4, 2), (3, 2, 5))
        sh = PerpShock3D(dims..., T(3), T(4), T(5); η = T(0.1), τ = T(0.2))
        for i in eachindex(sh.n)
            for c = 1:3
                sh.B[c][i] = T(sin(i + 0.3c))
                sh.u[c][i] = T(0.1cos(2i + c))
            end
        end
        compare_step(sh, T(0.7))
    end
end

@testset "shock CN analytic modes and accuracy" begin
    for T in (Float32, Float64), kind in (:advection, :diffusion), h in (0.125, 1.0, 4.0)
        sh = setup(T, kind)
        before = copy(sh.B[3])
        k = Float64(T(2π) / sh.Ly)
        lambda = kind == :advection ? -im * k : -k * k
        amp = (1 + h * lambda / 2) / (1 - h * lambda / 2)
        H._cn_b3d!(sh, T(h))
        for j = 1:sh.ny
            expected = 0.5 + 0.125real(amp * cis(2π * (j - 1) / sh.ny))
            @test maximum(abs, sh.B[3][:, j, :] .- expected) <= 1000eps(T)
        end
        energy_before = sum(abs2, before .- T(0.5))
        energy_after = sum(abs2, sh.B[3] .- T(0.5))
        @test isapprox(
            energy_after / energy_before,
            abs2(amp);
            atol = 2000eps(T),
            rtol = 2000eps(T),
        )
    end
    errors = Float64[]
    for nsteps in (4, 8, 16)
        sh = setup(Float64, :advection)
        ws = ShockCNWorkspace(sh)
        for _ = 1:nsteps
            H._cn_b3d!(sh, 1 / nsteps; workspace = ws)
        end
        exact = [0.5 + 0.125cos(2π * (j - 1) / sh.ny - 1) for j = 1:sh.ny]
        push!(errors, maximum(abs(sh.B[3][1, j, 1] - exact[j]) for j = 1:sh.ny))
    end
    @test 3.8 < errors[1] / errors[2] < 4.2
    @test 3.8 < errors[2] / errors[3] < 4.2
end

function public_fixture(::Type{T}) where {T}
    sh = setup(T, :mixed)
    ps = ParticleSet{3,T}(4)
    ps.x[1] .= T[0.01, 0.7, 1.4, 1.99]
    ps.x[2] .= T[0.2, 0.4, 1.2, 2.4]
    ps.x[3] .= T[0.3, 1.1, 2.1, 3.1]
    ps.v[1] .= T[-1, 0.2, -0.3, 1]
    ps.v[2] .= T[0.3, 0.4, 0.5, 0.6]
    ps.weight .= T[0.3, 0.7, 1.1, 0.9]
    ps.tag .= UInt32[1, 3, 7, 9]
    init_shock3d!(sh, ps)
    return sh, ps
end

same(sh, ps, a, b) = H._shock3d_restart_bitmatch(sh, ps, a, b)

@testset "shock CN public transaction, controls and restart" begin
    for T in (Float32, Float64)
        sh, ps = public_fixture(T)
        ws = ShockCNWorkspace(sh, ps)
        oldsh, oldps = deepcopy(sh), deepcopy(ps)
        @test_throws ShockCNConvergenceError step_shock3d!(
            sh,
            ps,
            T(0.2);
            NB = 2,
            field_method = :cn,
            cn_workspace = ws,
            cn_maxiter = 1,
        )
        @test same(sh, ps, oldsh, oldps)
        @test nparticles(ws.particles) < nparticles(ps)
        original_arrays = (sh.B..., sh.E..., ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
        before = flatten(sh.B)
        mktemp() do path, io
            close(io)
            checkpoint_shock3d(path, sh, ps)
            restored, particles = restore_shock3d(path)
            newws = ShockCNWorkspace(restored, particles)
            step_shock3d!(sh, ps, T(0.2); NB = 2, field_method = :cn, cn_workspace = ws)
            step_shock3d!(
                restored,
                particles,
                T(0.2);
                NB = 2,
                field_method = :cn,
                cn_workspace = newws,
            )
            @test same(sh, ps, restored, particles)
            @test nparticles(particles) < nparticles(oldps)
        end
        @test all(
            a === b for (a, b) in
            zip(original_arrays, (sh.B..., sh.E..., ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
        )
        @test nparticles(ps) < nparticles(oldps)
        # Assemble the frozen post-push operator from the public result; the
        # independent reference starts from the pre-step B, over both substeps.
        A, b = independent_operator(sh)
        half = Float64(T(0.2) / 2)
        reference = before
        for _ = 1:2
            reference = (I - half * A / 2) \ ((I + half * A / 2) * reference + half * b)
        end
        @test norm(flatten(sh.B) - reference, Inf) <= 3000eps(T) * max(norm(reference, Inf), 1)
        mktemp() do path, io
            close(io)
            checkpoint_shock3d(path, sh, ps)
            restored, particles = restore_shock3d(path)
            newws = ShockCNWorkspace(restored, particles)
            for (dt, nb) in ((T(0.03), 1), (T(0.07), 3), (T(0.02), 2))
                step_shock3d!(sh, ps, dt; NB = nb, field_method = :cn, cn_workspace = ws)
                step_shock3d!(
                    restored,
                    particles,
                    dt;
                    NB = nb,
                    field_method = :cn,
                    cn_workspace = newws,
                )
                @test same(sh, ps, restored, particles)
            end
        end
        for opts in (
            (; cn_rtol = -1),
            (; cn_rtol = Inf),
            (; cn_rtol = 1),
            (; cn_rtol = 0, cn_atol = 0),
            (; cn_atol = NaN),
            (; cn_maxiter = 0),
        )
            oldsh, oldps = deepcopy(sh), deepcopy(ps)
            @test_throws ArgumentError step_shock3d!(sh, ps, T(0.01); field_method = :cn, opts...)
            @test same(sh, ps, oldsh, oldps)
        end
        @test_throws ArgumentError step_shock3d!(sh, ps, T(0.01); cn_workspace = ws)
        @test_throws ArgumentError step_shock3d!(sh, ps, T(0.01); cn_maxiter = 3)
        @test_throws ArgumentError step_shock3d!(
            sh,
            ps,
            T(0.01);
            field_method = :cn,
            cn_workspace = 1,
        )
        @test_throws ArgumentError step_shock3d!(
            ws.trial,
            ws.particles,
            T(0.01);
            field_method = :cn,
            cn_workspace = ws,
        )
        oldsh, oldps = deepcopy(sh), deepcopy(ps)
        ws.solver.x = ws.rhs
        @test_throws ArgumentError step_shock3d!(
            sh,
            ps,
            T(0.01);
            field_method = :cn,
            cn_workspace = ws,
        )
        @test same(sh, ps, oldsh, oldps)
        ws = ShockCNWorkspace(sh, ps)
        sh.η = T(0.03)
        other, more = deepcopy(sh), deepcopy(ps)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        step_shock3d!(other, more, T(0.01); field_method = :cn)
        @test same(sh, ps, other, more)
        sh.dy *= 2
        @test_throws ArgumentError step_shock3d!(
            sh,
            ps,
            T(0.01);
            field_method = :cn,
            cn_workspace = ws,
        )
    end
end

@testset "shock CN rejects unsupported and aliased storage" begin
    for n in (0, 2)
        sh = setup(Float64, :diffusion)
        base = ParticleSet{3,Float64}(n)
        for x in base.x
            fill!(x, 0.5)
        end
        viewed = ParticleSet{3,Float64}(
            map(x -> view(x, :), base.x),
            base.v,
            base.weight,
            base.id,
            base.tag,
            base.q,
            base.m,
        )
        @test_throws ArgumentError ShockCNWorkspace(sh, viewed)
        @test_throws ArgumentError step_shock3d!(sh, viewed, 0.1; field_method = :cn)
        aliased = ParticleSet{3,Float64}(
            (base.x[1], base.x[1], base.x[3]),
            base.v,
            base.weight,
            base.id,
            base.tag,
            base.q,
            base.m,
        )
        @test_throws ArgumentError step_shock3d!(sh, aliased, 0.1; field_method = :cn)
    end
end

@testset "shock CN late arithmetic rejection and workspace recovery" begin
    for T in (Float32, Float64)
        sh, ps = public_fixture(T)
        fill!(ps.weight, T(100))
        sh.Te = floatmax(T)
        ws = ShockCNWorkspace(sh, ps)
        oldsh, oldps = deepcopy(sh), deepcopy(ps)
        @test_throws DomainError step_shock3d!(
            sh,
            ps,
            T(0.01);
            NB = 2,
            field_method = :cn,
            cn_workspace = ws,
        )
        @test same(sh, ps, oldsh, oldps)
        # Both field substeps have completed before final pressure/E rejects.
        @test ws.trial.B != oldsh.B
        sh.Te = zero(T)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        @test all(isfinite, flatten(sh.B))
        @test all(isfinite, flatten(sh.E))
        large = setup(T, :diffusion)
        fill!(large.B[3], floatmax(T) / 2)
        empty = ParticleSet{3,T}(0)
        snap = deepcopy(large)
        @test_throws ShockCNConvergenceError step_shock3d!(large, empty, T(1); field_method = :cn)
        @test same(large, empty, snap, empty)
    end
end

struct CNCustomVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(a::CNCustomVector) = size(a.data)
Base.getindex(a::CNCustomVector, i::Int) = a.data[i]
Base.setindex!(a::CNCustomVector, v, i::Int) = (a.data[i] = v)
@testset "shock CN interval and custom storage boundaries" begin
    for T in (Float32, Float64)
        sh = setup(T, :diffusion)
        ps = ParticleSet{3,T}(0)
        snap = deepcopy(sh)
        for nb in (1, 2, 3)
            @test_throws ArgumentError step_shock3d!(
                sh,
                ps,
                nextfloat(zero(T));
                NB = nb,
                field_method = :cn,
            )
            @test same(sh, ps, snap, ps)
        end
        for n in (0, 2)
            p = ParticleSet{3,T}(n)
            for x in p.x
                fill!(x, T(0.5))
            end
            q = ParticleSet{3,T}(map(CNCustomVector, p.x), p.v, p.weight, p.id, p.tag, p.q, p.m)
            @test_throws ArgumentError step_shock3d!(sh, q, T(0.1); field_method = :cn)
        end
    end
end

@testset "shock CN nonfinite norm and tolerance fail closed" begin
    for T in (Float32, Float64)
        # Singleton Fourier axes retain finite individual constant components;
        # their combined Euclidean norm overflows even though no entry does.
        sh = PerpShock3D(3, 1, 1, T(2), T(1), T(1); η = 0, τ = 0, Te = 0, B0 = 0)
        ps = ParticleSet{3,T}(0)
        for b in sh.B
            fill!(b, floatmax(T) / 2)
        end
        snap = deepcopy(sh)
        err = try
            step_shock3d!(sh, ps, zero(T); NB = 1, field_method = :cn)
            nothing
        catch ex
            ex
        end
        @test err isa ShockCNConvergenceError
        @test err.iterations == 0
        @test isinf(err.tolerance)
        @test same(sh, ps, snap, ps)
        for b in sh.B
            fill!(b, floatmax(T) * sqrt(eps(T)))
        end
        snap = deepcopy(sh)
        @test_throws ShockCNConvergenceError step_shock3d!(
            sh,
            ps,
            zero(T);
            NB = 1,
            field_method = :cn,
            cn_atol = floatmax(T),
            cn_rtol = T(0.5),
        )
        @test same(sh, ps, snap, ps)
    end
end

@testset "shock CN particle growth and later substep rollback" begin
    for T in (Float32, Float64)
        sh, ps = public_fixture(T)
        ws = ShockCNWorkspace(sh, ps)
        step_shock3d!(sh, ps, T(0.2); field_method = :cn, cn_workspace = ws)
        _, incoming = public_fixture(T)
        append_particles!(ps, incoming)
        init_shock3d!(sh, ps)
        other, particles = deepcopy(sh), deepcopy(ps)
        @test nparticles(ps) > nparticles(ws.particles)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        step_shock3d!(other, particles, T(0.01); field_method = :cn)
        @test same(sh, ps, other, particles)

        # Negative SAT damping is a supported configuration. Its scalar CN
        # amplification is three here: the first substep is representable,
        # while the second solution exceeds the floating-point range.
        sh = PerpShock3D(3, 1, 1, T(2), T(1), T(1); Te = 0, η = 0, τ = -0.5, B0 = 0)
        sh.B[3][3, 1, 1] = floatmax(T) / 8
        ps = ParticleSet{3,T}(0)
        ws = ShockCNWorkspace(sh, ps)
        before = deepcopy(sh)
        @test_throws ShockCNConvergenceError step_shock3d!(
            sh,
            ps,
            T(2);
            NB = 2,
            field_method = :cn,
            cn_workspace = ws,
        )
        @test same(sh, ps, before, ps)
        @test ws.trial.B[3][3, 1, 1] == 3before.B[3][3, 1, 1]
    end
end


@testset "shock CN solver history ownership and integer controls" begin
    for T in (Float32, Float64)
        for name in (:residuals, :Aresiduals, :qvals, :Acond),
            target in (:geometry, :particles, :operator, :solver, :empty)

            sh, ps = public_fixture(T)
            ws = ShockCNWorkspace(sh, ps)
            storage =
                target === :geometry ? sh.x :
                target === :particles ? ps.weight :
                target === :operator ? ws.trial.x :
                target === :solver ? ws.rhs :
                getfield(ws.solver.stats, name === :residuals ? :Aresiduals : :residuals)
            setfield!(ws.solver.stats, name, storage)
            before, particles = deepcopy(sh), deepcopy(ps)
            @test_throws ArgumentError step_shock3d!(
                sh,
                ps,
                T(0.01);
                field_method = :cn,
                cn_workspace = ws,
            )
            @test same(sh, ps, before, particles)
        end
        reference, particles = public_fixture(T)
        step_shock3d!(reference, particles, T(0.01); NB = 2, field_method = :cn)
        for nb in (Int8(2), Int128(2), big(2))
            sh, ps = public_fixture(T)
            step_shock3d!(sh, ps, T(0.01); NB = nb, field_method = :cn)
            @test same(sh, ps, reference, particles)
        end
        sh, ps = public_fixture(T)
        before, particles = deepcopy(sh), deepcopy(ps)
        @test_throws ArgumentError step_shock3d!(
            sh,
            ps,
            T(0.01);
            NB = big(typemax(Int)) + 1,
            field_method = :cn,
        )
        @test same(sh, ps, before, particles)
    end
end


@testset "shock CN stationary and tiny-step convergence" begin
    for T in (Float32, Float64),
        memory in (1, 3),
        dt in (zero(T), eps(T), T(0.125)),
        sat in (false, true)

        sh = PerpShock3D(
            3,
            1,
            1,
            T(2),
            T(1),
            T(1);
            η = T(0.1),
            τ = sat ? T(0.2) : zero(T),
            B0 = T(0.5),
            Te = 0,
        )
        if !sat
            fill!(sh.B[1], T(0.25))
            fill!(sh.B[2], T(-0.125))
        end
        ps = ParticleSet{3,T}(0)
        init_shock3d!(sh, ps)
        ws = ShockCNWorkspace(sh, ps; memory)
        beforeB, beforeE = deepcopy(sh.B), deepcopy(sh.E)
        step_shock3d!(sh, ps, dt; NB = 1, field_method = :cn, cn_workspace = ws)
        @test sh.B == beforeB
        @test sh.E == beforeE
    end
    for T in (Float32, Float64)
        compare_step(setup(T, :diffusion), eps(T))
        sh, ps = public_fixture(T)
        ws = ShockCNWorkspace(sh, ps)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        ps = ParticleSet{3,T}(0)
        sh.η = sh.τ = zero(T)
        init_shock3d!(sh, ps)
        beforeB = deepcopy(sh.B)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        @test sh.B == beforeB
        sh.η = T(0.1)
        other, particles = deepcopy(sh), deepcopy(ps)
        step_shock3d!(sh, ps, T(0.01); field_method = :cn, cn_workspace = ws)
        step_shock3d!(other, particles, T(0.01); field_method = :cn)
        @test same(sh, ps, other, particles)
    end
end

end # module
