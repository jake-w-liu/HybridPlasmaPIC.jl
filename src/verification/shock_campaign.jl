# shock_campaign.jl — 3-D shock campaign + restart (Phase 13 items):
#   • run selected M_A=4,6 3-D cases,
#   • cross-seed robustness (a 3-D effect is not inferred from one noisy run),
#   • numerical-variation robustness (conclusion unchanged under resolution/dt),
#   • controlled 1-D vs 3-D comparison at matched physical parameters,
#   • a production case that completes init → checkpoint/restart → analysis.
#
# step_shock3d! uses no RNG after the load, so a checkpoint→restart continues the
# run BITWISE-identically; production_3d_case asserts that explicitly.

"Serialize a 3-D shock state (`sh`,`ps`) for restart. Returns `path`."
function checkpoint_shock3d(path::AbstractString, sh::PerpShock3D, ps::ParticleSet{3})
    serialize(path, (sh, ps))
    return path
end

function _validate_shock3d_restart_state(state)
    state isa Tuple ||
        throw(ArgumentError("shock3d restart file does not contain a checkpoint tuple"))
    length(state) == 2 ||
        throw(ArgumentError("shock3d restart tuple must have exactly 2 entries, got $(length(state))"))
    sh, ps = state
    sh isa PerpShock3D ||
        throw(ArgumentError("shock3d restart entry 1 must be a PerpShock3D, got $(typeof(sh))"))
    ps isa ParticleSet{3} ||
        throw(ArgumentError("shock3d restart entry 2 must be a 3-D ParticleSet, got $(typeof(ps))"))
    return sh, ps
end

"Restore a `(sh, ps)` 3-D shock state written by [`checkpoint_shock3d`](@ref)."
function restore_shock3d(path::AbstractString)
    return _validate_shock3d_restart_state(deserialize(path))
end

@inline _state_value_equal(a, b) = a == b

@inline function _state_value_equal(a::Tuple, b::Tuple)
    length(a) == length(b) || return false
    @inbounds for i = 1:length(a)
        _state_value_equal(a[i], b[i]) || return false
    end
    return true
end

@inline _state_value_equal(a::SBP1D, b::SBP1D) = _typed_field_state_equal(a, b)

function _typed_field_state_equal(a, b)
    T = typeof(a)
    T === typeof(b) || return false
    for name in fieldnames(T)
        _state_value_equal(getfield(a, name), getfield(b, name)) || return false
    end
    return true
end

_particle_state_equal(a::ParticleSet{D}, b::ParticleSet{D}) where {D} =
    _typed_field_state_equal(a, b)

_shock3d_state_equal(a::PerpShock3D, b::PerpShock3D) = _typed_field_state_equal(a, b)

function _shock3d_restart_bitmatch(
    sh_a::PerpShock3D,
    ps_a::ParticleSet{3},
    sh_b::PerpShock3D,
    ps_b::ParticleSet{3},
)
    return _shock3d_state_equal(sh_a, sh_b) && _particle_state_equal(ps_a, ps_b)
end

# transverse-averaged downstream compression / frozen-in of a 3-D shock state
function _shock3d_diag(sh::PerpShock3D{T}) where {T}
    ny, nz = sh.ny, sh.nz
    nbar = dropdims(sum(sh.n; dims = (2, 3)); dims = (2, 3)) ./ (ny * nz)
    Bzbar = dropdims(sum(sh.B[3]; dims = (2, 3)); dims = (2, 3)) ./ (ny * nz)
    ipk = argmax(Bzbar)
    thr = (sh.B0 + Bzbar[ipk]) / 2
    ifr = sh.nx
    for i = ipk:sh.nx
        if Bzbar[i] < thr
            ifr = i
            break
        end
    end
    ilo = findfirst(>(T(3)), sh.x)
    ilo === nothing && (ilo = 1)
    slab = ilo:max(ilo, ifr - 1)
    n2 = sum(@view nbar[slab]) / length(slab)
    Bz2 = sum(@view Bzbar[slab]) / length(slab)
    return (; n2, Bz2, frozen_ratio = (Bz2 / sh.B0) / n2)
end

"""
    production_3d_case(; MA=4.0, nx=40, ny=8, nz=8, Lx=70, Ly=10, Lz=10, nppc=8,
                         nsteps_pre=150, nsteps_post=150, dt=0.03, seed=1)
        -> (; pass, restart_bitmatch, n2, frozen_ratio, maxdivB)

A full production 3-D pipeline: initialize, run `nsteps_pre` steps, checkpoint to
disk, restart from the checkpoint and run `nsteps_post` more, then analyze. The
restart is verified BITWISE-identical against an uninterrupted continuation
(`step_shock3d!` is deterministic), and the downstream is checked for a real
compressive shock. `pass` is the conjunction of both.
"""
function production_3d_case(;
    MA::Real = 4.0,
    nx::Integer = 40,
    ny::Integer = 8,
    nz::Integer = 8,
    Lx::Real = 70.0,
    Ly::Real = 10.0,
    Lz::Real = 10.0,
    nppc::Integer = 8,
    nsteps_pre::Integer = 150,
    nsteps_post::Integer = 150,
    dt::Real = 0.03,
    seed::Integer = 1,
)
    nsteps_pre >= 0 || throw(ArgumentError("nsteps_pre must be non-negative"))
    nsteps_post >= 0 || throw(ArgumentError("nsteps_post must be non-negative"))
    T = Float64
    _require_valid_positive_shock_ma(MA, T)
    sh, ps = _load_shock3d(; MA, nx, ny, nz, Lx, Ly, Lz, nppc, seed)
    for _ = 1:nsteps_pre
        step_shock3d!(sh, ps, T(dt); NB = 2)
    end
    # checkpoint, then make a deepcopy as the uninterrupted reference
    path = tempname()
    checkpoint_shock3d(path, sh, ps)
    shc = deepcopy(sh)
    psc = deepcopy(ps)
    for _ = 1:nsteps_post
        step_shock3d!(shc, psc, T(dt); NB = 2)
    end
    # restart from disk and continue
    sh2, ps2 = restore_shock3d(path)
    for _ = 1:nsteps_post
        step_shock3d!(sh2, ps2, T(dt); NB = 2)
    end
    rm(path; force = true)

    bitmatch = _shock3d_restart_bitmatch(sh2, ps2, shc, psc)
    d = _shock3d_diag(sh2)
    maxdivB = maximum(abs, magnetic_divergence3d(sh2))
    shock_ok = isfinite(d.n2) && d.n2 > 1.5 && all(isfinite, sh2.B[3])
    return (;
        pass = bitmatch && shock_ok,
        restart_bitmatch = bitmatch,
        d.n2,
        d.frozen_ratio,
        maxdivB,
    )
end

"""
    shock_campaign_3d(; MAs=(4.0, 6.0), seeds=(1, 2), kwargs...) -> Vector{NamedTuple}

Run the selected `MAs` (e.g. `M_A = 4, 6`) 3-D perpendicular-shock cases, each
over several `seeds`, and return one summary NamedTuple per `MA`:
`(; MA, n2_mean, n2_std, frozen_mean, robust)`. `robust` flags that the
cross-seed scatter is small relative to the mean (`n2_std < 0.15·n2_mean`) — the
quantitative statement that the result is not an artifact of a single noisy
realization. `kwargs` (grid size, steps, …) are forwarded to
[`run_perp_shock3d`](@ref).
"""
function shock_campaign_3d(; MAs = (4.0, 6.0), seeds = (1, 2), kwargs...)
    isempty(seeds) && throw(ArgumentError("seeds must be non-empty"))
    out = NamedTuple[]
    for MA in MAs
        n2s = Float64[]
        frs = Float64[]
        for s in seeds
            r = run_perp_shock3d(; MA = MA, seed = s, kwargs...)
            push!(n2s, r.n2)
            push!(frs, r.frozen_ratio)
        end
        μ = sum(n2s) / length(n2s)
        σ = length(n2s) > 1 ? sqrt(sum(abs2, n2s .- μ) / (length(n2s) - 1)) : 0.0
        push!(
            out,
            (;
                MA,
                n2_mean = μ,
                n2_std = σ,
                frozen_mean = sum(frs) / length(frs),
                robust = isfinite(μ) && σ < 0.15 * μ,
            ),
        )
    end
    return out
end

"""
    compare_dims_shock(; MA=3.0, kwargs...) -> (; oned, threed, frozen_consistent, both_compress)

Run the 1-D ([`run_perp_shock`](@ref)) and 3-D ([`run_perp_shock3d`](@ref))
perpendicular shocks at the SAME controlled physical parameters (`MA`, `Te`,
`vthi`, `η`) and compare. `frozen_consistent` checks both reproduce flux freezing
((Bz/B0)/n ≈ 1 to 8%); `both_compress` checks both give `n₂ > 1.5`. This is the
controlled 1-D/3-D comparison (matched physics, differing dimensionality).
"""
function compare_dims_shock(;
    MA::Real = 3.0,
    Te::Real = 0.125,
    vthi::Real = 0.35,
    η::Real = 0.02,
    oned_kwargs = (; N = 256, nsteps = 500),
    threed_kwargs = (;
        nx = 40,
        ny = 8,
        nz = 8,
        Lx = 70.0,
        Ly = 10.0,
        Lz = 10.0,
        nppc = 8,
        nsteps = 500,
        dt = 0.03,
    ),
)
    o = run_perp_shock(; MA = MA, Te = Te, vthi = vthi, η = η, oned_kwargs...)
    t = run_perp_shock3d(; MA = MA, Te = Te, vthi = vthi, η = η, threed_kwargs...)
    frozen_consistent = abs(o.frozen_ratio - 1) < 0.08 && abs(t.frozen_ratio - 1) < 0.08
    both_compress = o.n2 > 1.5 && t.n2 > 1.5
    return (; oned = o, threed = t, frozen_consistent, both_compress)
end
