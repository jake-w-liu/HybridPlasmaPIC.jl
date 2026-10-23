"A shock CN solve failed its true residual check; the public step was not committed."
struct ShockCNConvergenceError{T} <: Exception
    residual::T
    tolerance::T
    iterations::Int
end

function Base.showerror(io::IO, err::ShockCNConvergenceError)
    print(
        io,
        "shock CN solve failed after ",
        err.iterations,
        " iterations: residual ",
        err.residual,
        " exceeds finite tolerance ",
        err.tolerance,
    )
end

const _ShockCNParticles{T} =
    ParticleSet{3,T,Vector{T},Vector{T},Vector{T},Vector{UInt64},Vector{UInt32}}

"""
    ShockCNWorkspace(sh, ps=ParticleSet{3,T}(0); memory=20)

Reusable CPU trial state and restarted GMRES storage for `step_shock3d!` with
`field_method=:cn`. Particle components must be ordinary `Vector`s, including
empty populations. Reuse with the same grid dimensions and precision; current
physical state and mutable parameters are copied at each step. Particle count
may change between calls. Workspace contents are scratch, may change on failure,
and are not checkpoint state. Use a separate workspace per concurrent task.
"""
struct ShockCNWorkspace{T<:Union{Float32,Float64}}
    trial::PerpShock3D{T}
    particles::_ShockCNParticles{T}
    solver::Krylov.GmresWorkspace{T,T,Vector{T}}
    rhs::Vector{T}
    guess::Vector{T}
    residual::Vector{T}
end

_shock_cn_arrays(sh::PerpShock3D) = (
    sh.B...,
    sh.n,
    sh.u...,
    sh.pe,
    sh.E...,
    sh.J...,
    sh.Eind...,
    sh.k1...,
    sh.k2...,
    sh.k3...,
    sh.k4...,
    sh.Btmp...,
    sh.s,
    sh.wsum,
    sh.mx,
    sh.my,
    sh.mz,
)
_shock_cn_geometry(sh::PerpShock3D) = (sh.x, sh.y, sh.z, sh.sbp.H)

function _shock_cn_controls(::Type{T}, rtol, atol, maxiter) where {T}
    r = _require_finite_nonnegative_real("cn_rtol", rtol, T)
    a = _require_finite_nonnegative_real("cn_atol", atol, T)
    r < one(T) || throw(ArgumentError("cn_rtol must be below one"))
    r > zero(T) || a > zero(T) || throw(ArgumentError("at least one CN tolerance must be positive"))
    return r, a, _require_positive_intlike("cn_maxiter", maxiter)
end

function _validate_shock_cn_particles(ps::ParticleSet{3,T}; check_values::Bool = true) where {T}
    ps isa _ShockCNParticles{T} ||
        throw(ArgumentError("CN requires ordinary CPU Vector particle storage"))
    _check_particle_storage_axes(ps)
    _particle_components_mightalias(ps, _particle_backing_mightalias) &&
        throw(ArgumentError("CN particle components must not share backing storage"))
    if check_values
        for a in (ps.x..., ps.v..., ps.weight)
            all(isfinite, a) || throw(ArgumentError("CN particle values must be finite"))
        end
    end
    return nothing
end

function _validate_shock_cn_state(sh::PerpShock3D{T}; check_values::Bool = true) where {T}
    T in (Float32, Float64) || throw(ArgumentError("CN supports Float32 and Float64"))
    sh.nx >= 3 && sh.ny >= 1 && sh.nz >= 1 || throw(ArgumentError("invalid CN grid dimensions"))
    dims = (sh.nx, sh.ny, sh.nz)
    for a in _shock_cn_arrays(sh)
        size(a) == dims || throw(DimensionMismatch("CN shock arrays must match $dims"))
    end
    for (a, n) in zip(_shock_cn_geometry(sh), (sh.nx, sh.ny, sh.nz, sh.nx))
        length(a) == n || throw(DimensionMismatch("CN grid geometry length mismatch"))
        all(isfinite, a) || throw(ArgumentError("CN grid geometry must be finite"))
    end
    for (name, value) in (("Lx", sh.Lx), ("Ly", sh.Ly), ("Lz", sh.Lz), ("nfloor", sh.nfloor))
        _require_finite_positive_real(name, value, T)
    end
    _require_finite_nonnegative_real("Te", sh.Te, T)
    _require_finite_nonnegative_real("η", sh.η, T)
    _require_valid_gamma(sh.γe, T)
    _require_finite_real("τ", sh.τ, T)
    _require_finite_real("B0", sh.B0, T)
    dx = sh.Lx / (sh.nx - 1)
    sh.sbp.n == sh.nx && sh.sbp.dx == dx && sh.dy == sh.Ly / sh.ny && sh.dz == sh.Lz / sh.nz ||
        throw(ArgumentError("CN grid spacings and dimensions are inconsistent"))
    isfinite(inv(dx)) && dx / 2 > zero(T) && sh.dy > zero(T) && sh.dz > zero(T) ||
        throw(ArgumentError("CN grid spacings must be representable"))
    for i in eachindex(sh.sbp.H)
        sh.sbp.H[i] == (i == 1 || i == sh.nx ? dx / 2 : dx) ||
            throw(ArgumentError("CN SBP norm is inconsistent with the grid"))
    end
    arrays = (_shock_cn_arrays(sh)..., _shock_cn_geometry(sh)...)
    for i in eachindex(arrays), j = 1:(i-1)
        Base.mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("CN shock arrays must not share backing storage"))
    end
    if check_values
        for a in (sh.B..., sh.E..., sh.u...)
            all(isfinite, a) || throw(ArgumentError("CN fields and velocity must be finite"))
        end
    end
    return nothing
end

function ShockCNWorkspace(
    sh::PerpShock3D{T},
    ps::ParticleSet{3,T} = ParticleSet{3,T}(0);
    memory::Integer = 20,
) where {T}
    _validate_shock_cn_state(sh)
    _validate_shock_cn_particles(ps)
    m = _require_positive_intlike("CN memory", memory)
    n = 3length(sh.n)
    solver = Krylov.GmresWorkspace(n, n, Vector{T}; memory = min(m, n))
    # Restart and the explicit initial guess use Δx; size it at construction.
    resize!(solver.Δx, n)
    return ShockCNWorkspace{T}(
        deepcopy(sh),
        deepcopy(ps),
        solver,
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
    )
end

function _validate_shock_cn_workspace(ws, sh::PerpShock3D, ps::ParticleSet)
    throw(ArgumentError("cn_workspace must have the same precision as the shock"))
end

function _validate_shock_cn_workspace(
    ws::ShockCNWorkspace{T},
    sh::PerpShock3D{T},
    ps::ParticleSet{3,T},
) where {T}
    _validate_shock_cn_state(sh)
    _validate_shock_cn_particles(ps)
    size(ws.trial.n) == size(sh.n) ||
        throw(DimensionMismatch("CN workspace grid dimensions differ"))
    _validate_shock_cn_state(ws.trial; check_values = false)
    _validate_shock_cn_particles(ws.particles; check_values = false)
    n = 3length(sh.n)
    s = ws.solver
    k = length(s.V)
    s.m == n &&
        s.n == n &&
        1 <= k <= n &&
        length(s.c) == k &&
        length(s.s) == k &&
        length(s.z) == k &&
        length(s.R) == k * (k + 1) ÷ 2 &&
        isempty(s.p) &&
        isempty(s.q) || throw(DimensionMismatch("CN Krylov workspace dimensions differ"))
    buffers = (ws.rhs, ws.guess, ws.residual, s.x, s.w, s.Δx, s.V...)
    all(a -> length(a) == n, buffers) || throw(DimensionMismatch("CN solver vector lengths differ"))
    # Krylov resets every statistics vector even with history collection off.
    # Their backing storage must be independent before any solver call.
    storage = (
        buffers...,
        s.c,
        s.s,
        s.z,
        s.R,
        s.stats.residuals,
        s.stats.Aresiduals,
        s.stats.qvals,
        s.stats.Acond,
        _shock_cn_arrays(ws.trial)...,
        _shock_cn_geometry(ws.trial)...,
        _shock_cn_arrays(sh)...,
        _shock_cn_geometry(sh)...,
    )
    for i in eachindex(storage)
        for j = 1:(i-1)
            (storage[i] === storage[j] || Base.mightalias(storage[i], storage[j])) &&
                throw(ArgumentError("CN workspace and shock storage must be independent"))
        end
        _mightalias_particle_storage(storage[i], ps, _particle_backing_mightalias) &&
            throw(ArgumentError("CN workspace/shock must not alias caller particles"))
        _mightalias_particle_storage(storage[i], ws.particles, _particle_backing_mightalias) &&
            throw(ArgumentError("CN workspace/shock must not alias trial particles"))
    end
    _particle_sets_mightalias(ps, ws.particles, _particle_backing_mightalias) &&
        throw(ArgumentError("CN trial and caller particles must have independent storage"))
    return nothing
end

function _copy_shock_cn_state!(dst::PerpShock3D{T}, src::PerpShock3D{T}) where {T}
    for (to, from) in zip(_shock_cn_arrays(dst), _shock_cn_arrays(src))
        copyto!(to, from)
    end
    for (to, from) in zip(_shock_cn_geometry(dst), _shock_cn_geometry(src))
        copyto!(to, from)
    end
    dst.sbp = SBP1D{T}(src.sbp.n, src.sbp.dx, dst.sbp.H)
    dst.nx, dst.ny, dst.nz = src.nx, src.ny, src.nz
    dst.Lx, dst.Ly, dst.Lz = src.Lx, src.Ly, src.Lz
    dst.dy, dst.dz = src.dy, src.dz
    dst.Te, dst.γe, dst.η = src.Te, src.γe, src.η
    dst.τ, dst.B0, dst.nfloor = src.τ, src.B0, src.nfloor
    return dst
end

function _copy_shock_cn_particles!(dst::_ShockCNParticles, src::_ShockCNParticles)
    n = length(src.weight)
    for (to, from) in zip((dst.x..., dst.v..., dst.weight), (src.x..., src.v..., src.weight))
        resize!(to, n)
        copyto!(to, from)
    end
    resize!(dst.id, n)
    resize!(dst.tag, n)
    copyto!(dst.id, src.id)
    copyto!(dst.tag, src.tag)
    dst.q, dst.m = src.q, src.m
    return dst
end

struct _ShockCNOperator{T<:Union{Float32,Float64}}
    shock::PerpShock3D{T}
    half_dt::T
end
Base.size(A::_ShockCNOperator) = (3length(A.shock.n), 3length(A.shock.n))
Base.eltype(::_ShockCNOperator{T}) where {T} = T

function LinearAlgebra.mul!(out::Vector{T}, A::_ShockCNOperator{T}, x::Vector{T}) where {T}
    T === Float32 ||
        T === Float64 ||
        throw(ArgumentError("CN operator requires a concrete Float32 or Float64 precision"))
    sh = A.shock
    n = length(sh.n)
    for c = 1:3
        copyto!(sh.Btmp[c], 1, x, (c - 1) * n + 1, n)
    end
    _b_rhs3d!(sh.k2, sh.Btmp, sh, zero(T))
    @inbounds for c = 1:3, i = 1:n
        q = (c - 1) * n + i
        out[q] = x[q] - A.half_dt * sh.k2[c][i]
    end
    return out
end

function _cn_b3d!(
    sh::PerpShock3D{T},
    h::T;
    workspace::ShockCNWorkspace{T} = ShockCNWorkspace(sh),
    rtol::Real = 100eps(T),
    atol::Real = 0,
    maxiter::Integer = 200,
) where {T<:Union{Float32,Float64}}
    rtolT, atolT, limit = _shock_cn_controls(T, rtol, atol, maxiter)
    isfinite(h) && h >= zero(T) ||
        throw(ArgumentError("CN interval must be finite and nonnegative"))
    ws = workspace
    n = length(sh.n)
    length(ws.rhs) == 3n || throw(DimensionMismatch("CN workspace grid dimensions differ"))
    _b_rhs3d!(sh.k1, sh.B, sh)
    half = h / 2
    h > zero(T) &&
        half == zero(T) &&
        throw(ArgumentError("CN half interval is not representable; reduce subcycling"))
    @inbounds for c = 1:3, i = 1:n
        q = (c - 1) * n + i
        ws.guess[q] = sh.B[c][i]
        ws.rhs[q] = sh.B[c][i] + half * sh.k1[c][i]
    end
    # f(B)=L*B+b; the other half of the constant SAT forcing belongs in rhs.
    forcing = half * (sh.τ / sh.sbp.H[sh.nx]) * sh.B0
    @inbounds for k = 1:sh.nz, j = 1:sh.ny
        q = 2n + sh.nx + (j - 1) * sh.nx + (k - 1) * sh.nx * sh.ny
        ws.rhs[q] += forcing
    end
    # Storage is a CPU Vector{Float32/64}; use its scaled BLAS Euclidean norm
    # directly, including overflow/nonfinite handling checked below.
    rhsnorm = BLAS.nrm2(ws.rhs)
    tolerance = atolT + rtolT * rhsnorm
    isfinite(rhsnorm) && isfinite(tolerance) && all(isfinite, ws.guess) ||
        throw(ShockCNConvergenceError(T(Inf), tolerance, 0))
    A = _ShockCNOperator(sh, half)
    # Test the initial state against the same true-residual contract. Besides
    # avoiding unnecessary iteration for stationary/tiny steps, this prevents
    # Krylov 0.10.11's zero-residual restarted warm start from adding x₀ twice.
    mul!(ws.residual, A, ws.guess)
    ws.residual .-= ws.rhs
    initial_residual = BLAS.nrm2(ws.residual)
    isfinite(initial_residual) || throw(ShockCNConvergenceError(initial_residual, tolerance, 0))
    initial_residual <= tolerance && return sh
    # Solve to a tighter estimated residual, then independently check the true
    # residual. Explicit B_old as the guess removes hidden restart history.
    Krylov.gmres!(
        ws.solver,
        A,
        ws.rhs,
        ws.guess;
        restart = true,
        reorthogonalization = true,
        atol = tolerance / 2,
        rtol = zero(T),
        itmax = limit,
    )
    solution = Krylov.solution(ws.solver)
    mul!(ws.residual, A, solution)
    ws.residual .-= ws.rhs
    residual = BLAS.nrm2(ws.residual)
    isfinite(residual) && all(isfinite, solution) && residual <= tolerance ||
        throw(ShockCNConvergenceError(residual, tolerance, ws.solver.stats.niter))
    for c = 1:3
        copyto!(sh.B[c], 1, solution, (c - 1) * n + 1, n)
    end
    return sh
end
