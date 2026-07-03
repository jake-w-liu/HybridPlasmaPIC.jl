# camcl.jl — CAM-CL hybrid integrator (Matthews, J. Comput. Phys. 112, 102 (1994)).
#
# A SECOND integrator for the SAME massless-electron hybrid model that
# `HybridStepper`/`step!` advances (generalized Ohm's law
#   E = −u_i×B + (J×B)/n − ∇p_e/n + η J − ηH ∇²J,  ∂B/∂t = −∇×E),
# on a periodic FourierGrid. It differs from `HybridStepper` in two ways that
# define the Matthews scheme:
#
#   • CAM  (Current Advance Method): the ion current/moments are advanced to the
#     HALF step n+1/2 so that the electric field used to push B is time-centred.
#     Here the half-step moments are deposited from the midpoint positions
#     x^{n+1/2} = x^n + (dt/2) v^{n+1/2} together with the freshly Boris-kicked
#     velocities v^{n+1/2} — the standard CAM "advanced current" — and are then
#     held FROZEN across the whole magnetic subcycle.
#
#   • CL  (Cyclic LeapFrog): the magnetic field is subcycled with a
#     time-reversible, non-dissipative cyclic leapfrog rather than RK4. Two
#     staggered copies of B (offset by one substep) are advanced alternately by
#     a full leapfrog kick each, "cycling" which copy leads; a final averaging of
#     the two copies removes the odd-even (leapfrog) splitting error and lands
#     both at n+1. CL is the magnetic mover Matthews pairs with CAM because it is
#     cheap (one Ohm/curl evaluation per substep, vs four for RK4) and conserves
#     magnetic energy far better than a dissipative scheme over many substeps.
#
# Documented time levels over one step! advancing all by dt:
#   particle position x   : integer n            (→ n+1)
#   particle velocity v   : half    n−1/2        (→ n+1/2)
#   ion moments (n, u_i)  : half    n+1/2  (FROZEN through the B subcycle)
#   magnetic field   B    : integer n            (→ n+1)
#   electric field   E    : integer n  (carried) (recomputed → n+1)
# The Boris kick is centred at n; the B advance is centred at n+1/2.

"""
    CAMCLStepper(g, model, shape, N::Integer)

Workspace and configuration for the CAM-CL hybrid integrator advancing a single
proton species of `N` particles on grid `g` under `model` with deposition/gather
`shape`. Set the initial field in `stepper.fields.B`, then call [`init_camcl!`](@ref)
before stepping with [`step_camcl!`](@ref). Advances the same hybrid model as
[`HybridStepper`](@ref) using Matthews' Current-Advance-Method + Cyclic-LeapFrog.
"""
mutable struct CAMCLStepper{D,T,SH<:ShapeFunction,M<:HybridModel}
    g::FourierGrid{D,T}
    model::M
    shape::SH
    fields::HybridFields{D,T}          # canonical B, carried E, plus Ohm scratch
    fn::Array{T,D}                     # frozen n^{n+1/2}
    fui::NTuple{3,Array{T,D}}          # frozen u_i^{n+1/2}
    Escr::NTuple{3,Array{T,D}}         # E scratch for the B subcycle
    Blead::NTuple{3,Array{T,D}}        # leading CL copy of B
    Blag::NTuple{3,Array{T,D}}         # lagging CL copy of B
    rhs::NTuple{3,Array{T,D}}          # −∇×E scratch
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    vpred::NTuple{3,Vector{T}}         # predicted v^{n+1} for the carried-E re-centering (2nd order)
    xmid::NTuple{D,Vector{T}}
    work::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
end

function CAMCLStepper(
    g::FourierGrid{D,T},
    model::M,
    shape::SH,
    N::Integer,
) where {D,T,M<:HybridModel,SH<:ShapeFunction}
    # the CL subcycle freezes the scalar ∇p_e (f.gradp), which is incompatible with the
    # B-dependent gyrotropic force — use HybridStepper for anisotropic (CGL) closures.
    is_anisotropic(model.closure) && throw(
        ArgumentError(
            "CAMCLStepper does not support anisotropic (CGL) closures; use HybridStepper",
        ),
    )
    nc = g.n
    Np = _particle_length(N)
    CAMCLStepper{D,T,SH,M}(
        g,
        model,
        shape,
        HybridFields{D,T}(nc),
        zeros(T, nc),
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

"""
    init_camcl!(stepper, ps)

Compute the carried electric field E⁰ from the initial particle moments and
`stepper.fields.B`. Call once after loading particles and setting B, before
[`step_camcl!`](@ref). Analogous to [`init!`](@ref) for `HybridStepper`.
"""
function init_camcl!(st::CAMCLStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    nf = T(st.model.nfloor)
    _moments!(st.fields.n, st.fields.ui, ps, st.g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    st.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    st.step[] = 0
    return st
end

# F = −∇×E(Btrial; frozen n+1/2 moments) into out, reusing Ohm scratch.
@inline function _camcl_rhs!(
    out::NTuple{3,<:Array{T,D}},
    Btrial::NTuple{3,<:Array{T,D}},
    st::CAMCLStepper{D,T},
) where {D,T}
    f = st.fields
    # frozen pe/∇pe/ninv precomputed once per step (via _ohm_prep! in step_camcl!).
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
    _apply_electron_inertia!(st.Escr, T(st.model.de2), st.g)   # E ← E/(1+d_e²k²) if de2>0
    faraday_rhs!(out, st.Escr, st.g)
    return out
end

# Cyclic leapfrog of the magnetic field from n to n+1 over NB substeps of size h
# (Matthews 1994, §3.2, "cyclic leapfrog"). CL integrates the first-order field
# ODE  dB/dt = F(B) ≡ −∇×E(B; frozen moments)  with a leapfrog by carrying two
# copies of B staggered by one substep h, so that the RHS is always evaluated at
# the time level CENTRED for the copy being kicked (the source of leapfrog's
# second-order accuracy and its time-reversibility / non-dissipation).
#
# Time levels of the two copies (in units of h, both starting at B^n = level 0):
#   • Half-step bootstrap puts the half-copy C at level 1/2 (integer copy A @ 0),
#     staggering them by h/2.
#   • Staggered leapfrog (NB steps of h): each kick is time-centred on the OTHER
#     copy — A_{k+1}=A_k+h·F(C_{k+1/2}), C_{k+3/2}=C_{k+1/2}+h·F(A_{k+1}) — giving
#     leapfrog's second-order accuracy and time-reversibility. After NB steps A is
#     at level NB and C at level NB−1/2.
#   • A final half-step synchronizes C to level NB; the two copies (both now at
#     n+1) are AVERAGED, which cancels the odd/even (computational-mode) leapfrog
#     splitting and lands the result centred EXACTLY at n+1 — Matthews' reason for
#     the final average and for CL conserving magnetic energy far better than a
#     dissipative scheme. (The prior full-step-bootstrap + (NB−1)×2h-kick variant
#     averaged levels NB and NB−1, landing at n+1−h/2: an O(h) centring error.)
#
# Stability: like any leapfrog, CL has a CFL bound tighter than RK4's; the
# stiff whistler branch needs enough subcycles (NB) that ω_whistler·h ≲ 1.
# Advances B from integer level n to integer level n+1, time-centred at n+1/2
# with the frozen CAM moments — the level structure step_camcl! documents.
function _cl_subcycle_B!(st::CAMCLStepper{D,T}, dtT::T, NB::Integer) where {D,T}
    NB >= 2 || throw(ArgumentError("CAM-CL needs NB ≥ 2 magnetic subcycles (got $NB)"))
    B = st.fields.B
    h = dtT / NB

    A = st.Blead        # integer-level copy
    C = st.Blag         # half-level copy
    for c = 1:3
        copyto!(A[c], B[c])
        copyto!(C[c], B[c])
    end
    # half-step bootstrap: C → level 1/2 using F at level 0.
    _camcl_rhs!(st.rhs, A, st)
    for c = 1:3
        @. C[c] += (h / 2) * st.rhs[c]
    end
    # staggered leapfrog, NB steps of h, each kick time-centred:
    #   A_{k+1} = A_k + h·F(C_{k+1/2});   C_{k+3/2} = C_{k+1/2} + h·F(A_{k+1}).
    @inbounds for k = 0:(NB-1)
        _camcl_rhs!(st.rhs, C, st)              # F at half-level k+1/2
        for c = 1:3
            @. A[c] += h * st.rhs[c]            # A: k → k+1
        end
        if k < NB - 1
            _camcl_rhs!(st.rhs, A, st)          # F at integer level k+1
            for c = 1:3
                @. C[c] += h * st.rhs[c]        # C: k+1/2 → k+3/2
            end
        end
    end
    # A @ NB, C @ NB−1/2: synchronize C to NB (final half-step), then average the
    # two copies — cancels the odd/even leapfrog computational mode and lands
    # centred EXACTLY at n+1 (verified exact for a constant RHS).
    _camcl_rhs!(st.rhs, A, st)
    for c = 1:3
        @. C[c] += (h / 2) * st.rhs[c]
        @. B[c] = T(0.5) * (A[c] + C[c])
    end
    return B
end

"""
    step_camcl!(stepper, ps, dt; NB=2)

Advance the plasma one timestep `dt` using the CAM-CL scheme with `NB` cyclic-
leapfrog magnetic subcycles (periodic box; `NB ≥ 2`, since the cyclic leapfrog
needs the bootstrap plus at least one leapfrog kick). Mirrors [`step!`](@ref)
for `HybridStepper`: Boris kick centred at n, CAM half-step ion current, CL
magnetic subcycle, then recompute the carried E at n+1.
"""
function step_camcl!(
    st::CAMCLStepper{D,T},
    ps::ParticleSet{D,T},
    dt::Real;
    NB::Integer = 2,
) where {D,T}
    dtT = _validated_step_dt(T, dt, NB; min_NB = 2, name = "step_camcl!")
    g = st.g
    h = dtT / 2
    nf = T(st.model.nfloor)
    lo = ntuple(_ -> zero(T), D)
    hi = g.L

    # 1. Boris push with carried E^n, B^n:  v^{n−1/2}→v^{n+1/2}, x^n→x^{n+1}.
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
    qm = ps.q / ps.m
    # prime the leapfrog once: loaded v is physical v^0 → v^{-1/2} for 2nd-order accuracy.
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
            st.xmid[d][p] = ps.x[d][p] + h * ps.v[d][p]   # x^{n+1/2}
            ps.x[d][p] += dtT * ps.v[d][p]                # x^{n+1}
        end
    end
    apply_periodic!(ps, lo, hi)

    # 2. CAM: advanced (half-step) ion current/moments from x^{n+1/2}, v^{n+1/2}.
    #    Frozen across the magnetic subcycle so E is time-centred at n+1/2.
    psmid = ParticleSet{D,T}(st.xmid, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    apply_periodic!(psmid, lo, hi)
    _moments!(st.fn, st.fui, psmid, g, st.shape, nf, st.work)
    # frozen-moment Ohm terms (pe, ∇pe, 1/n) computed ONCE; reused every subcycle.
    f = st.fields
    _ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, st.model.closure, T(st.model.nfloor), f.floor_count, g)

    # 3. CL: cyclic-leapfrog subcycle of B from n → n+1.
    _cl_subcycle_B!(st, dtT, NB)

    # 4. recompute carried E = E^{n+1}, then re-center u_i to integer level n+1 via a predictor
    #    half-kick so the carried E is 2nd-order accurate (consistent with HybridStepper).
    _moments!(st.fields.n, st.fields.ui, ps, g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    _recenter_carried_E!(st, ps, dtT)

    st.time[] += dtT
    st.step[] += 1
    return st
end
