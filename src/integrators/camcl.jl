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
#     copies of B, staggered by h/2 inside the subcycle (kick-drift-kick form),
#     PERSIST across particle steps and are advanced alternately by leapfrog
#     kicks; the field the rest of the step consumes is their (output-only)
#     average at n+1, which removes the odd-even (leapfrog) splitting error —
#     see _cl_subcycle_B! for why persistence is essential to CL's neutral
#     stability. CL is the magnetic mover Matthews pairs with CAM because it is
#     cheap (one Ohm/curl evaluation per copy per substep, vs four for RK4) and
#     conserves magnetic energy far better than a dissipative scheme over many
#     substeps.
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
    Blead::NTuple{3,Array{T,D}}        # persistent CL copy, integer substep levels
    Blag::NTuple{3,Array{T,D}}         # persistent CL copy, half levels inside the subcycle
    Bout::NTuple{3,Array{T,D}}         # last CL output B (detects external edits of fields.B)
    rhs::NTuple{3,Array{T,D}}          # −∇×E scratch
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    vpred::NTuple{3,Vector{T}}         # predicted v^{n+1} for the carried-E re-centering (2nd order)
    xmid::NTuple{D,Vector{T}}
    work::Vector{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    cl_h::Base.RefValue{T}             # substep size h the persistent CL pair was advanced with
    cl_live::Base.RefValue{Bool}       # persistent CL pair valid to continue from?
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
        ntuple(_ -> zeros(T, nc), Val(3)),   # fui
        ntuple(_ -> zeros(T, nc), Val(3)),   # Escr
        ntuple(_ -> zeros(T, nc), Val(3)),   # Blead
        ntuple(_ -> zeros(T, nc), Val(3)),   # Blag
        ntuple(_ -> zeros(T, nc), Val(3)),   # Bout
        ntuple(_ -> zeros(T, nc), Val(3)),   # rhs
        ntuple(_ -> zeros(T, Np), Val(3)),   # Ep
        ntuple(_ -> zeros(T, Np), Val(3)),   # Bp
        ntuple(_ -> zeros(T, Np), Val(3)),   # vpred
        ntuple(_ -> zeros(T, Np), Val(D)),   # xmid
        zeros(T, Np),
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),                        # cl_h
        Ref(false),                          # cl_live
    )
end

"""
    init_camcl!(stepper, ps)

Compute the carried electric field E⁰ from the initial particle moments and
`stepper.fields.B`. Call once after loading particles and setting B, before
[`step_camcl!`](@ref). Analogous to [`init!`](@ref) for `HybridStepper`.
"""
function init_camcl!(st::CAMCLStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))
    nf = T(st.model.nfloor)
    _moments!(st.fields.n, st.fields.ui, ps, st.g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    st.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    st.step[] = 0
    st.cl_live[] = false                # re-bootstrap the persistent CL pair from fields.B
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
# copies of B staggered by h/2, so that the RHS is always evaluated at the time
# level CENTRED for the copy being kicked (the source of leapfrog's second-order
# accuracy and its time-reversibility / non-dissipation).
#
# The two copies PERSIST across particle steps (A = st.Blead, C = st.Blag; both
# land on the integer level n+1 at the end of every step, staggered by h/2 only
# INSIDE the subcycle — a kick-drift-kick / Verlet arrangement of the staggered
# leapfrog). Persistence is what makes the pair neutrally stable: the per-step
# map (opening half C-kick, NB A-kicks interleaved with NB−1 C-kicks, closing
# half C-kick, all unit-determinant shears) has unit-modulus eigenvalues for
# every whistler mode with ω·h ≤ 2 (verified on the exact 2×2 mode map and in
# code, see below). The previous variant re-initialized BOTH copies from the
# single averaged B every particle step; that extra projection makes the
# composition weakly unstable at ALL ω·h > 0 (measured per-step |m| = 1.008 at
# ω·h = 1 and 1.11 at ω·h = 1.6 — production runs at recommended_dt(:camcl) blew
# up in ~700–2000 steps).
#
# Second-order accuracy across the step seam: the closing half C-kick of step n
# uses F(A^{n+1}) with step-n moments (frozen at n+1/2) and the opening half
# C-kick of step n+1 uses F(A^{n+1}) with step-(n+1) moments (frozen at n+3/2);
# their sum is one full kick centred exactly on the boundary n+1 evaluated with
# the AVERAGE of the two adjacent frozen-moment fields — the moment-lag errors
# cancel to O(dt²). (Keeping C permanently staggered across the seam instead
# biases its kicks by h/2 against the moment centring and degrades the whole
# stepper to 1st order — measured convergence ratios ~2 instead of ~4.)
#
# Output: the field the rest of the step consumes is Matthews' copy average
# B^{n+1} = ½·(A + C), both at n+1, which cancels the odd/even (computational-
# mode) leapfrog splitting. The average is OUTPUT-ONLY: it is never fed back
# into the persistent pair, so it cannot re-create the restart instability; the
# (bounded, |eig| = 1) computational mode merely rides along in (A, C).
#
# Re-bootstrap (A = C = B, one Matthews averaging event) happens only when the
# pair cannot be continued: first step after init_camcl!, a change of substep
# size h, an external edit of fields.B (detected against the stored last output
# Bout), or copy drift max|A − C| > 10% of max|A|. The drift trigger is
# Matthews' "average the copies occasionally": under real-axis damping
# (η/ηH > 0) the leapfrog computational mode grows as e^{+νh} per substep
# (det = 1: whatever damps the physical mode grows the parasitic one), and the
# on-demand re-sync arrests it; in non-dissipative runs it never fires and CL
# stays exactly neutral.
#
# Stability: neutrally stable for ω_whistler·h < 2 (measured on the stiffest
# n=64 mode over 10⁴ steps: bounded non-normal beat envelopes with per-step
# multiplier ≤ 1+1e-5, saturated by mid-run; non-finite at ω·h = 2.05), tighter
# than RK4's 2.8 but two Ohm evaluations per substep instead of four. Advances B
# from integer level n to integer level n+1, time-centred at n+1/2 with the
# frozen CAM moments — the level structure step_camcl! documents.
function _cl_subcycle_B!(st::CAMCLStepper{D,T}, dtT::T, NB::Integer) where {D,T}
    NB >= 2 || throw(ArgumentError("CAM-CL needs NB ≥ 2 magnetic subcycles (got $NB)"))
    B = st.fields.B
    h = dtT / NB

    A = st.Blead        # persistent copy, integer substep levels
    C = st.Blag         # persistent copy, half levels inside the subcycle
    if !(st.cl_live[] && st.cl_h[] == h && all(B[c] == st.Bout[c] for c = 1:3))
        # (re)bootstrap: both copies from B at level 0 (a Matthews averaging
        # event — rare by construction, see the drift trigger above).
        for c = 1:3
            copyto!(A[c], B[c])
            copyto!(C[c], B[c])
        end
    end
    # opening half-kick: C → level 1/2 using F(A @ 0) — with the previous step's
    # closing half-kick this forms one boundary-centred full kick (see above).
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
    # A @ NB, C @ NB−1/2: closing half-kick synchronizes C to NB inside the
    # persistent state; output B = ½(A + C) (Matthews' average, cancels the
    # odd/even computational mode) WITHOUT feeding it back into the pair. Also
    # measure the copies' drift to decide whether to re-sync next step.
    _camcl_rhs!(st.rhs, A, st)
    drift = zero(T)
    scale = zero(T)
    @inbounds for c = 1:3
        Ac, Cc, Rc, Bc, Oc = A[c], C[c], st.rhs[c], B[c], st.Bout[c]
        for i in eachindex(Ac, Cc, Rc, Bc, Oc)
            cs = Cc[i] + (h / 2) * Rc[i]
            Cc[i] = cs
            d = abs(Ac[i] - cs)
            drift = ifelse(d > drift, d, drift)
            a = abs(Ac[i])
            scale = ifelse(a > scale, a, scale)
            b = T(0.5) * (Ac[i] + cs)
            Bc[i] = b
            Oc[i] = b
        end
    end
    st.cl_h[] = h
    st.cl_live[] = drift <= T(0.1) * scale      # else: re-sync (re-bootstrap) next step
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
    iszero(dtT) && return st        # dt=0 no-op: do not consume the one-time priming guard.
    qm = _validated_qm(ps)
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))
    g = st.g
    h = dtT / 2
    nf = T(st.model.nfloor)
    lo = ntuple(_ -> zero(T), D)
    hi = g.L

    # 1. Boris push with carried E^n, B^n:  v^{n−1/2}→v^{n+1/2}, x^n→x^{n+1}.
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
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
