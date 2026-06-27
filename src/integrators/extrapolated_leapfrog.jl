# integrator.jl — explicit hybrid leapfrog with a subcycled magnetic field.
#
# Documented time levels (the cycle below advances all of these by one dt):
#   particle position x        : integer n          (→ n+1)
#   particle velocity v        : half    n−1/2      (→ n+1/2)
#   magnetic field   B         : integer n          (→ n+1)
#   electric field   E (carried): integer n         (recomputed → n+1)
#   moments for the B-subcycle : n+1/2 (midpoint positions x^{n+1/2}, v^{n+1/2})
#
# One step (CAM/CL-style; the B subcycle is RK4 for stability of the stiff
# whistler branch, with the n+1/2 ion moments held frozen across the subcycle):
#   1. Gather E^n, B^n to particles; Boris push  v^{n−1/2}→v^{n+1/2}, x^n→x^{n+1}.
#      Record midpoint positions x^{n+1/2}; wrap positions (periodic box).
#   2. Deposit frozen moments n^{n+1/2}, u_i^{n+1/2} from (x^{n+1/2}, v^{n+1/2}).
#   3. Subcycle B: integrate ∂B/∂t = −∇×E(B; frozen moments) over dt with N_B RK4
#      substeps  →  B^{n+1}.
#   4. Recompute carried E from n^{n+1}, u_i^{n+1/2}, B^{n+1}  →  E^{n+1}.
# The Boris update is centered at n; the B update is centered at n+1/2.

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
        HybridFields{D,T}(nc),
        zeros(T, nc),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, nc), Val(3)),
        ntuple(_ -> zeros(T, Np), Val(3)),
        ntuple(_ -> zeros(T, Np), Val(3)),
        ntuple(_ -> zeros(T, Np), Val(D)),
        zeros(T, Np),
        Ref(zero(T)),
        Ref(0),
    )
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
    nf = T(st.model.nfloor)
    _moments!(st.fields.n, st.fields.ui, ps, st.g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)
    return st
end

# F = −∇×E(Btrial; frozen moments) into out
@inline function _bfield_rhs!(
    out::NTuple{3,<:Array{T,D}},
    Btrial::NTuple{3,<:Array{T,D}},
    st::HybridStepper{D,T},
) where {D,T}
    f = st.fields
    # frozen pe/∇pe/ninv were precomputed once per step (via _ohm_prep! in step!);
    # only the B-dependent part is evaluated here, each subcycle stage.
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
    dtT = T(dt)
    isfinite(dtT) || throw(ArgumentError("$name requires finite dt (got $dt)"))
    dtT >= zero(T) || throw(ArgumentError("$name requires nonnegative dt (got $dt)"))
    return dtT
end

"""
    step!(stepper, ps, dt; NB=1)

Advance the plasma one timestep `dt` with `NB` magnetic subcycles (periodic box).
"""
function step!(st::HybridStepper{D,T}, ps::ParticleSet{D,T}, dt::Real; NB::Integer = 1) where {D,T}
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step!")
    g = st.g
    h = dtT / 2
    nf = T(st.model.nfloor)
    lo = ntuple(_ -> zero(T), D)
    hi = g.L

    # 1. push particles with carried E^n, B^n
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
    qm = ps.q / ps.m
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
    # frozen-moment Ohm terms (pe, ∇pe, 1/n) computed ONCE; reused every subcycle.
    f = st.fields
    _ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, st.model.closure, T(st.model.nfloor), f.floor_count, g)

    # 3. subcycle B: n → n+1
    hb = dtT / NB
    for _ = 1:NB
        _rk4_B!(st, hb)
    end

    # 4. recompute carried E = E^{n+1}
    _moments!(st.fields.n, st.fields.ui, ps, g, st.shape, nf, st.work)
    ohms_law!(st.fields, st.model, st.g)

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
