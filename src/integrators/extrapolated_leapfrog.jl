# integrator.jl — explicit hybrid leapfrog with a subcycled magnetic field.
#
# Documented time levels (the cycle below advances all of these by one dt):
#   particle position x        : integer n          (→ n+1)
#   particle velocity v        : half    n−1/2      (→ n+1/2)
#   magnetic field   B         : integer n          (→ n+1)
#   electric field   E (carried): integer n         (recomputed → n+1)
#   moments for the B-subcycle : n+1/2 (midpoint positions x^{n+1/2}, v^{n+1/2})
#
# The loaded velocity is the physical v^0, so the FIRST step primes the leapfrog once
# (v^{-1/2} = v^0 − (dt/2)·a^0) — without it the loaded v^0 acts as v^{-1/2} and the run is only
# 1st-order in dt. This scheme is 2nd-order in dt (verified by a convergence study, incl. on
# magnetized problems); getting there needs both the priming and the step-4 u_i re-centering below.
#
# One step (CAM/CL-style; the B subcycle is RK4 for stability of the stiff
# whistler branch, with the n+1/2 ion moments held frozen across the subcycle):
#   1. Gather E^n, B^n to particles; Boris push  v^{n−1/2}→v^{n+1/2}, x^n→x^{n+1}.
#      Record midpoint positions x^{n+1/2}; wrap positions (periodic box).
#   2. Deposit frozen moments n^{n+1/2}, u_i^{n+1/2} from (x^{n+1/2}, v^{n+1/2}).
#   3. Subcycle B: integrate ∂B/∂t = −∇×E(B; frozen moments) over dt with N_B RK4
#      substeps  →  B^{n+1}.
#   4. Recompute carried E from n^{n+1}, u_i^{n+1/2}, B^{n+1}; then RE-CENTER u_i to n+1 via a
#      predictor half-kick and recompute → 2nd-order E^{n+1} (u_i^{n+1/2} alone would be 1st-order
#      in the −u_i×B convection term). See _recenter_carried_E!.
# The Boris update is centered at n; the B update is centered at n+1/2; the scheme is 2nd-order.

"""
    HybridStepper(g, model, shape, N::Integer)

Workspace and configuration for time-stepping a single proton species of `N`
particles on grid `g` under `model` with deposition/gather `shape`. Set the
initial field in `stepper.fields.B`, then call [`init!`](@ref) before stepping.
"""
mutable struct HybridStepper{D,T,SH<:ShapeFunction,M<:HybridModel}
    g::FourierGrid{D,T}
    model::M
    shape::SH
    fields::HybridFields{D,T}          # canonical B, carried E, plus Ohm scratch
    fn::Array{T,D}                     # frozen n^{n+1/2}
    fui::NTuple{3,Array{T,D}}          # frozen u_i^{n+1/2}
    Escr::NTuple{3,Array{T,D}}         # E scratch for the B-subcycle
    k1::NTuple{3,Array{T,D}}
    k2::NTuple{3,Array{T,D}}
    k3::NTuple{3,Array{T,D}}
    k4::NTuple{3,Array{T,D}}
    Btmp::NTuple{3,Array{T,D}}
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    vpred::NTuple{3,Vector{T}}         # predicted v^{n+1} for the carried-E re-centering (2nd order)
    xmid::NTuple{D,Vector{T}}
    work::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
end

function HybridStepper(
    g::FourierGrid{D,T},
    model::M,
    shape::SH,
    N::Integer,
) where {D,T,M<:HybridModel,SH<:ShapeFunction}
    nc = g.n
    Np = _particle_length(N)
    HybridStepper{D,T,SH,M}(
        g,
        model,
        shape,
        HybridFields{D,T}(nc; anisotropic = is_anisotropic(model.closure)),
        zeros(T, nc),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, Np), Val(3)),   # Ep
        ntuple(_ -> zeros(T, Np), Val(3)),   # Bp
        ntuple(_ -> zeros(T, Np), Val(3)),   # vpred
        ntuple(_ -> zeros(T, Np), Val(D)),   # xmid
        zeros(T, Np),
        Ref(zero(T)),
        Ref(0),
    )
end

function _resize_hybrid_particle_workspaces!(st, n::Integer)
    n >= 0 || throw(ArgumentError("particle workspace length must be nonnegative, got $n"))
    N = Int(n)
    D = length(st.xmid)
    T = eltype(st.work)
    length(st.work) == N &&
        all(length(st.Ep[c]) == N for c = 1:3) &&
        all(length(st.Bp[c]) == N for c = 1:3) &&
        all(length(st.vpred[c]) == N for c = 1:3) &&
        all(length(st.xmid[d]) == N for d = 1:D) &&
        return st

    st.Ep = ntuple(_ -> Vector{T}(undef, N), 3)
    st.Bp = ntuple(_ -> Vector{T}(undef, N), 3)
    st.vpred = ntuple(_ -> Vector{T}(undef, N), 3)
    st.xmid = ntuple(_ -> Vector{T}(undef, N), D)
    st.work = Vector{T}(undef, N)
    return st
end

# moments (density into nout, bulk velocity into uout) with density floor
function _moments!(
    nout::Array{T,D},
    uout::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor::T,
    work::Vector{T},
) where {D,T}
    density!(nout, ps, g, shape)
    momentum!(uout, ps, g, shape; work)
    @inbounds for I in eachindex(nout)
        inv = one(T) / max(nout[I], nfloor)
        uout[1][I] *= inv
        uout[2][I] *= inv
        uout[3][I] *= inv
    end
    return nout
end

"""
    init!(stepper, ps)

Compute the carried electric field E⁰ from the initial particle moments and
`stepper.fields.B`. Call once after loading particles and setting B.
"""
function init!(st::HybridStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))
    nf = T(st.model.nfloor)
    _moments!(st.fields.n, st.fields.ui, ps, st.g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    st.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    st.step[] = 0
    return st
end

# F = −∇×E(Btrial; frozen moments) into out
@inline function _bfield_rhs!(
    out::NTuple{3,<:Array{T,D}},
    Btrial::NTuple{3,<:Array{T,D}},
    st::HybridStepper{D,T},
) where {D,T}
    f = st.fields
    # frozen pe/∇pe/ninv (scalar) or the frozen ∇·P_e force (anisotropic/CGL) were
    # precomputed once per step; only the B-dependent part is evaluated here each subcycle
    # stage. The `is_anisotropic` test is compile-time-resolved from the closure type, so
    # scalar closures keep the exact original call below.
    if is_anisotropic(st.model.closure)
        # the gyrotropic force depends on the field DIRECTION b, so it must track the trial
        # B within the subcycle (a step-frozen b would misalign as B_⊥ grows and drive a
        # numerical instability); recompute ∇·P_e(n^{n+1/2}, B_trial) each stage.
        anisotropic_pressure_force!(f.pforce, st.fn, Btrial, st.model.closure, st.g)
        _ohm_Efield_aniso!(
            st.Escr,
            st.fui,
            Btrial,
            f.J,
            f.lapJ,
            f.pforce,
            f.ninv,
            T(st.model.η),
            T(st.model.ηH),
            st.g,
        )
    else
        _ohm_Efield!(
            st.Escr,
            st.fui,
            Btrial,
            f.J,
            f.lapJ,
            f.gradp,
            f.ninv,
            T(st.model.η),
            T(st.model.ηH),
            st.g,
        )
    end
    _apply_electron_inertia!(st.Escr, T(st.model.de2), st.g)   # E ← E/(1+d_e²k²) if de2>0
    faraday_rhs!(out, st.Escr, st.g)
    return out
end

# one RK4 substep of size h on canonical B (st.fields.B), frozen moments
function _rk4_B!(st::HybridStepper{D,T}, h::T) where {D,T}
    B = st.fields.B
    _bfield_rhs!(st.k1, B, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k1[c]
    end
    _bfield_rhs!(st.k2, st.Btmp, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + (h / 2) * st.k2[c]
    end
    _bfield_rhs!(st.k3, st.Btmp, st)
    for c = 1:3
        @. st.Btmp[c] = B[c] + h * st.k3[c]
    end
    _bfield_rhs!(st.k4, st.Btmp, st)
    for c = 1:3
        @. B[c] += (h / 6) * (st.k1[c] + 2 * st.k2[c] + 2 * st.k3[c] + st.k4[c])
    end
    return B
end

function _validated_step_dt(
    ::Type{T},
    dt::Real,
    NB::Integer;
    min_NB::Integer,
    name::AbstractString,
) where {T}
    NB >= min_NB ||
        throw(ArgumentError("$name requires NB >= $min_NB magnetic subcycles (got $NB)"))
    return _validated_nonnegative_dt(T, dt; name)
end

"""
    step!(stepper, ps, dt; NB=1)

Advance the plasma one timestep `dt` with `NB` magnetic subcycles (periodic box).
"""
# One-time leapfrog priming. The loaded velocity is the physical v^0, but the leapfrog Boris push
# expects v^{-1/2}; push v back a half step v^{-1/2} = v^0 − (dt/2)·a^0 (forward-Euler, with
# a^0 = qm(E^0 + v^0×B^0) from the just-gathered fields) so the first Boris kick is centred at
# t=0. Applied once (first step after init!) — a one-time O(dt²) IC correction that restores the
# documented 2nd-order accuracy for a real init!+step! run (verified: end-to-end rate 1 → 2).
# `rel=true` (relativistic push, works in momentum u=γv) backs up MOMENTUM: u^{-1/2}=γ^0 v^0 − h·a^0,
# v^{-1/2}=u^{-1/2}/γ(u); the default non-relativistic path is byte-identical for the hybrid callers.
@inline function _prime_leapfrog!(v, Ep, Bp, qm, h, np, rel::Bool = false, c = one(eltype(v[1])))
    vx, vy, vz = v
    ex, ey, ez = Ep
    bx, by, bz = Bp
    if rel
        c2 = c * c
        o = one(c)
        @inbounds for p = 1:np
            v0x, v0y, v0z = vx[p], vy[p], vz[p]
            β2 = (v0x * v0x + v0y * v0y + v0z * v0z) / c2
            β2 >= o && (β2 = o - eps(c))
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
        @inbounds for p = 1:np
            ux, uy, uz = vx[p], vy[p], vz[p]
            vx[p] = ux - h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
            vy[p] = uy - h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
            vz[p] = uz - h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
        end
    end
    return v
end

# Forward half-kick predictor v^{n+1/2} → v^{n+1} into `vp`, from particle-gathered E, B.
@inline function _predict_half_kick!(vp, v, Ep, Bp, qm, h, np)
    vx, vy, vz = v
    vpx, vpy, vpz = vp
    ex, ey, ez = Ep
    bx, by, bz = Bp
    @inbounds for p = 1:np
        ux, uy, uz = vx[p], vy[p], vz[p]
        vpx[p] = ux + h * qm * (ex[p] + (uy * bz[p] - uz * by[p]))
        vpy[p] = uy + h * qm * (ey[p] + (uz * bx[p] - ux * bz[p]))
        vpz[p] = uz + h * qm * (ez[p] + (ux * by[p] - uy * bx[p]))
    end
    return vp
end

# Re-center the carried E to integer level n+1. The E just computed in step 4 deposits u_i from
# the half-step velocity v^{n+1/2}, so the −u_i×B convection term lags by dt/2 — an O(dt) force
# error that makes the whole scheme 1st-order on magnetized problems. Predict v^{n+1} with a
# forward half-kick from the fresh E^{n+1}, B^{n+1}, re-deposit u_i^{n+1}, and recompute the
# carried E, restoring the documented CAM/CL 2nd-order accuracy (verified: magnetized temporal
# convergence rate 1 → 2). Unmagnetized behavior is unchanged (u_i drops out of Ohm's law when
# B=0). Generic over the stepper (HybridStepper/CAMCLStepper share Ep/Bp/vpred/work). No alloc.
function _recenter_carried_E!(st, ps::ParticleSet{D,T}, dtT::T) where {D,T}
    qm = _validated_qm(ps)
    g = st.g
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)   # E^{n+1}
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)   # B^{n+1}
    _predict_half_kick!(st.vpred, ps.v, st.Ep, st.Bp, qm, dtT / 2, length(ps.weight))
    pspred = ParticleSet{D,T}(ps.x, st.vpred, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    _moments!(st.fields.n, st.fields.ui, pspred, g, st.shape, T(st.model.nfloor), st.work)
    ohms_law!(st.fields, st.model, g)
    return st
end

function step!(st::HybridStepper{D,T}, ps::ParticleSet{D,T}, dt::Real; NB::Integer = 1) where {D,T}
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step!")
    iszero(dtT) && return st        # dt=0 is a true no-op: advance nothing, and (critically) do
    #                                 NOT consume the step==0 one-time leapfrog-priming guard.
    qm = _validated_qm(ps)
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))
    g = st.g
    h = dtT / 2
    nf = T(st.model.nfloor)
    lo = ntuple(_ -> zero(T), D)
    hi = g.L

    # 1. push particles with carried E^n, B^n
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
    # prime the leapfrog once: loaded v is physical v^0 → v^{-1/2} for 2nd-order (see above).
    st.step[] == 0 && _prime_leapfrog!(ps.v, st.Ep, st.Bp, qm, h, length(ps.weight))
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = boris_kick(
            vx[p],
            vy[p],
            vz[p],
            st.Ep[1][p],
            st.Ep[2][p],
            st.Ep[3][p],
            st.Bp[1][p],
            st.Bp[2][p],
            st.Bp[3][p],
            qm,
            dtT,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            st.xmid[d][p] = ps.x[d][p] + h * ps.v[d][p]
            ps.x[d][p] += dtT * ps.v[d][p]
        end
    end
    apply_periodic!(ps, lo, hi)

    # 2. frozen n+1/2 moments from midpoint positions + v^{n+1/2}
    psmid = ParticleSet{D,T}(st.xmid, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    apply_periodic!(psmid, lo, hi)
    _moments!(st.fn, st.fui, psmid, g, st.shape, nf, st.work)
    # frozen-moment Ohm terms computed ONCE; reused every subcycle. Scalar closures freeze
    # pe/∇pe/1/n; the anisotropic (CGL) closure freezes 1/n and the pressure-tensor force
    # ∇·P_e(n^{n+1/2}, B^n) instead (compile-time-resolved from the closure type).
    f = st.fields
    if is_anisotropic(st.model.closure)
        # only 1/n is frozen; the gyrotropic force ∇·P_e is recomputed per subcycle stage
        # (it depends on the trial B direction — see _bfield_rhs!).
        _ohm_ninv!(f.ninv, st.fn, T(st.model.nfloor), f.floor_count)
    else
        _ohm_prep!(
            f.pe,
            f.gradp,
            f.ninv,
            st.fn,
            st.model.closure,
            T(st.model.nfloor),
            f.floor_count,
            g,
        )
    end

    # 3. subcycle B: n → n+1
    hb = dtT / NB
    for _ = 1:NB
        _rk4_B!(st, hb)
    end

    # 4. recompute carried E = E^{n+1} (n^{n+1}, u_i^{n+1/2}, B^{n+1}), then re-center u_i to n+1
    #    via a predictor half-kick so the carried E is 2nd-order accurate (see _recenter_carried_E!).
    _moments!(st.fields.n, st.fields.ui, ps, g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    _recenter_carried_E!(st, ps, dtT)

    st.time[] += dtT
    st.step[] += 1
    return st
end

# ---------------------------------------------------------------- diagnostics

"Total ion kinetic energy Σ_p ½ m w |v|²."
function kinetic_energy(ps::ParticleSet{D,T}) where {D,T}
    s = zero(T)
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        s += ps.weight[p] * (vx[p]^2 + vy[p]^2 + vz[p]^2)
    end
    return T(0.5) * ps.m * s
end

"Magnetic energy ∫ ½|B|² dV (normalized units)."
function magnetic_energy(B::NTuple{3,<:Array{T,D}}, g::FourierGrid{D,T}) where {D,T}
    _require_grid_tuple(:B, B, g)
    s = zero(T)
    @inbounds for c = 1:3, I in eachindex(B[c])
        s += B[c][I]^2
    end
    return T(0.5) * s * prod(g.dx)
end

"Electron internal energy ∫ p_e/(γ−1) dV (polytropic γ≠1; isothermal has no closed invariant)."
function electron_internal_energy(
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_array(:n, n, g)
    # a gyrotropic (CGL) closure has no scalar pressure ⇒ no scalar internal energy from this
    # (density-only) interface; degrade to NaN like isothermal rather than throwing downstream.
    is_anisotropic(closure) && return T(NaN)
    γ = T(closure_gamma(closure))
    γ == one(T) && return T(NaN)
    pe = similar(n)
    electron_pressure!(pe, n, closure)
    return sum(pe) / (γ - one(T)) * prod(g.dx)
end

"""
    mode_amplitude(field, g, m::NTuple{D,Int})

Complex Fourier amplitude Σ_g field_g · exp(−i k·x_g) of integer mode `m`
(k_d = 2π m_d / L_d), for tracking a single wave mode over time.
"""
function mode_amplitude(field::Array{T,D}, g::FourierGrid{D,T}, m::NTuple{D,Int}) where {D,T}
    _require_grid_array(:field, field, g)
    acc = zero(Complex{T})
    kk = ntuple(d -> T(2π) * m[d] / g.L[d], D)
    @inbounds for I in CartesianIndices(field)
        t = Tuple(I)
        phase = zero(T)
        for d = 1:D
            phase += kk[d] * (t[d] - 1) * g.dx[d]
        end
        acc += field[I] * cis(-phase)
    end
    return acc
end
