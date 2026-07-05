# hall_mhd.jl — periodic single-fluid Hall-MHD with massless fluid electrons.
#
# Normalized equations on the same FourierGrid operators as the hybrid PIC model:
#   ∂n/∂t = -∇·(n u)
#   ∂u/∂t = -u·∇u + (J×B - ∇p_i - ∇p_e)/n
#   J = ∇×B
#   E = -u×B + (J×B)/n - ∇p_e/n + ηJ - ηH∇²J
#   ∂B/∂t = -∇×E
#
# Ions are a fluid with isothermal pressure p_i = Ti n. Electrons use the same
# ElectronClosure hierarchy as HybridModel. Spatial dimension D is parametric;
# u and B always carry three vector components.

"""
    HallMHDModel(closure; Ti=0.0, η=0.0, ηH=0.0, nfloor=1e-6)

Parameters for periodic Hall-MHD. `closure` gives the electron pressure,
`Ti` is the isothermal ion temperature (`p_i = Ti*n`), `η` and `ηH` are
resistivity and hyperresistivity in Ohm's law, and `nfloor` protects `1/n`.
"""
struct HallMHDModel{C<:ElectronClosure}
    closure::C
    Ti::Float64
    η::Float64
    ηH::Float64
    nfloor::Float64
end

function HallMHDModel(closure::ElectronClosure; Ti = 0.0, η = 0.0, ηH = 0.0, nfloor = 1e-6)
    TiT = _require_finite_nonnegative_real("Ti", Ti, Float64)
    ηT = _require_finite_nonnegative_real("η", η, Float64)
    ηHT = _require_finite_nonnegative_real("ηH", ηH, Float64)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, Float64)
    return HallMHDModel(closure, TiT, ηT, ηHT, nfloorT)
end

"""
    HallMHDState(g, model)

State and workspaces for a periodic Hall-MHD solve on `g`. Initialize
`state.fields.n`, `state.fields.ui`, and `state.fields.B`, then call
[`hall_mhd_ohms_law!`](@ref), [`hall_mhd_rhs!`](@ref), or
[`step_hall_mhd!`](@ref).
"""
mutable struct HallMHDState{D,T,M<:HallMHDModel}
    g::FourierGrid{D,T}
    model::M
    fields::HybridFields{D,T}
    stage::HybridFields{D,T}
    rhs_n::Array{T,D}
    rhs_u::NTuple{3,Array{T,D}}
    rhs_B::NTuple{3,Array{T,D}}
    n0::Array{T,D}
    u0::NTuple{3,Array{T,D}}
    B0::NTuple{3,Array{T,D}}
    ntmp::Array{T,D}
    utmp::NTuple{3,Array{T,D}}
    Btmp::NTuple{3,Array{T,D}}
    kn::NTuple{4,Array{T,D}}
    ku::NTuple{4,NTuple{3,Array{T,D}}}
    kB::NTuple{4,NTuple{3,Array{T,D}}}
    flux_n::NTuple{D,Array{T,D}}
    gradn::NTuple{D,Array{T,D}}
    gradu::NTuple{D,Array{T,D}}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
end

function HallMHDState(g::FourierGrid{D,T}, model::M) where {D,T,M<:HallMHDModel}
    nc = g.n
    z() = zeros(T, nc)
    vec3() = ntuple(_ -> z(), 3)
    kvec3() = ntuple(_ -> vec3(), 4)
    return HallMHDState{D,T,M}(
        g,
        model,
        HybridFields{D,T}(nc; anisotropic = is_anisotropic(model.closure)),
        HybridFields{D,T}(nc; anisotropic = is_anisotropic(model.closure)),
        z(),
        vec3(),
        vec3(),
        z(),
        vec3(),
        vec3(),
        z(),
        vec3(),
        vec3(),
        ntuple(_ -> z(), 4),
        kvec3(),
        kvec3(),
        ntuple(_ -> z(), D),
        ntuple(_ -> z(), D),
        ntuple(_ -> z(), D),
        Ref(zero(T)),
        Ref(0),
    )
end

function _copy_hall_state!(
    n_dst::Array{T,D},
    u_dst::NTuple{3,<:Array{T,D}},
    B_dst::NTuple{3,<:Array{T,D}},
    n_src::Array{T,D},
    u_src::NTuple{3,<:Array{T,D}},
    B_src::NTuple{3,<:Array{T,D}},
) where {D,T}
    copyto!(n_dst, n_src)
    for c = 1:3
        copyto!(u_dst[c], u_src[c])
        copyto!(B_dst[c], B_src[c])
    end
    return nothing
end

function _hall_stage!(
    nout::Array{T,D},
    uout::NTuple{3,<:Array{T,D}},
    Bout::NTuple{3,<:Array{T,D}},
    n0::Array{T,D},
    u0::NTuple{3,<:Array{T,D}},
    B0::NTuple{3,<:Array{T,D}},
    kn::Array{T,D},
    ku::NTuple{3,<:Array{T,D}},
    kB::NTuple{3,<:Array{T,D}},
    a::T,
) where {D,T}
    @. nout = n0 + a * kn
    for c = 1:3
        @. uout[c] = u0[c] + a * ku[c]
        @. Bout[c] = B0[c] + a * kB[c]
    end
    return nothing
end

function _hall_final_stage!(
    nout::Array{T,D},
    uout::NTuple{3,<:Array{T,D}},
    Bout::NTuple{3,<:Array{T,D}},
    st::HallMHDState{D,T},
    dt::T,
) where {D,T}
    k1n, k2n, k3n, k4n = st.kn
    @. nout = st.n0 + (dt / 6) * (k1n + 2 * k2n + 2 * k3n + k4n)
    for c = 1:3
        k1u, k2u, k3u, k4u = (st.ku[1][c], st.ku[2][c], st.ku[3][c], st.ku[4][c])
        k1B, k2B, k3B, k4B = (st.kB[1][c], st.kB[2][c], st.kB[3][c], st.kB[4][c])
        @. uout[c] = st.u0[c] + (dt / 6) * (k1u + 2 * k2u + 2 * k3u + k4u)
        @. Bout[c] = st.B0[c] + (dt / 6) * (k1B + 2 * k2B + 2 * k3B + k4B)
    end
    return nothing
end

function _validate_hall_candidate!(st::HallMHDState{D,T}, n, u, B) where {D,T}
    nfloor = T(st.model.nfloor)
    @inbounds for x in n
        isfinite(x) || throw(ArgumentError("Hall-MHD density went non-finite"))
        x > nfloor || throw(ArgumentError("Hall-MHD density fell below nfloor"))
    end
    for c = 1:3
        _require_all_finite("Hall-MHD velocity component $c", u[c], "step_hall_mhd!")
        _require_all_finite("Hall-MHD magnetic component $c", B[c], "step_hall_mhd!")
    end
    return nothing
end

function _hall_ohms_law_fields!(
    f::HybridFields{D,T},
    model::HallMHDModel,
    g::FourierGrid{D,T},
) where {D,T}
    if is_anisotropic(model.closure)
        _ohm_ninv!(f.ninv, f.n, T(model.nfloor), f.floor_count)
        anisotropic_pressure_force!(f.pforce, f.n, f.B, model.closure, g)
        _ohm_Efield_aniso!(
            f.E,
            f.ui,
            f.B,
            f.J,
            f.lapJ,
            f.pforce,
            f.ninv,
            T(model.η),
            T(model.ηH),
            g,
        )
    else
        _ohm_E!(
            f.E,
            f.n,
            f.ui,
            f.B,
            f.J,
            f.lapJ,
            f.pe,
            f.gradp,
            f.ninv,
            model.closure,
            T(model.η),
            T(model.ηH),
            T(model.nfloor),
            f.floor_count,
            g,
        )
    end
    return f
end

"""
    hall_mhd_ohms_law!(state)

Update `state.fields.E`, `J`, electron pressure, pressure gradient, reciprocal
density, and density-floor count from the current Hall-MHD `n`, `u`, and `B`.
"""
function hall_mhd_ohms_law!(st::HallMHDState{D,T}) where {D,T}
    _hall_ohms_law_fields!(st.fields, st.model, st.g)
    return st
end

function _stage_ohms_law!(
    st::HallMHDState{D,T},
    n::Array{T,D},
    u::NTuple{3,<:Array{T,D}},
    B::NTuple{3,<:Array{T,D}},
) where {D,T}
    _copy_hall_state!(st.stage.n, st.stage.ui, st.stage.B, n, u, B)
    _hall_ohms_law_fields!(st.stage, st.model, st.g)
    return st.stage
end

function _hall_mhd_rhs!(
    dn::Array{T,D},
    du::NTuple{3,<:Array{T,D}},
    dB::NTuple{3,<:Array{T,D}},
    st::HallMHDState{D,T},
    n::Array{T,D},
    u::NTuple{3,<:Array{T,D}},
    B::NTuple{3,<:Array{T,D}},
) where {D,T}
    g = st.g
    model = st.model
    stage = _stage_ohms_law!(st, n, u, B)

    for d = 1:D
        @. st.flux_n[d] = n * u[d]
    end
    divergence!(dn, st.flux_n, g)
    dn .*= -one(T)

    faraday_rhs!(dB, stage.E, g)

    gradient!(st.gradn, n, g)
    Ti = T(model.Ti)
    nfloor = T(model.nfloor)
    aniso = is_anisotropic(model.closure)
    Bx, By, Bz = B
    Jx, Jy, Jz = stage.J
    for c = 1:3
        for d = 1:D
            deriv!(st.gradu[d], u[c], g, d)
        end
        @inbounds for I in eachindex(n)
            adv = zero(T)
            for d = 1:D
                adv += u[d][I] * st.gradu[d][I]
            end
            lorentz = if c == 1
                Jy[I] * Bz[I] - Jz[I] * By[I]
            elseif c == 2
                Jz[I] * Bx[I] - Jx[I] * Bz[I]
            else
                Jx[I] * By[I] - Jy[I] * Bx[I]
            end
            # total pressure force: ion ∇p_i = Ti∇n on the spatial components, plus
            # the electron pressure the stage Ohm's law computed (scalar ∇p_e on the
            # spatial components, CGL ∇·P_e on all three) — the same force it puts
            # into E as −∇p_e/n, which must also act on the quasineutral bulk fluid.
            pressure = c <= D ? Ti * st.gradn[c][I] : zero(T)
            if aniso
                pressure += stage.pforce[c][I]
            elseif c <= D
                pressure += stage.gradp[c][I]
            end
            du[c][I] = -adv + (lorentz - pressure) / max(n[I], nfloor)
        end
    end
    return (; dn, du, dB)
end

"""
    hall_mhd_rhs!(state) -> (; dn, du, dB)

Compute the Hall-MHD right-hand side from the current state into reusable
workspaces and return them.
"""
function hall_mhd_rhs!(st::HallMHDState{D,T}) where {D,T}
    return _hall_mhd_rhs!(st.rhs_n, st.rhs_u, st.rhs_B, st, st.fields.n, st.fields.ui, st.fields.B)
end

"""
    step_hall_mhd!(state, dt; project_b=true)

Advance one explicit RK4 step for periodic Hall-MHD. The candidate final state is
checked for finite values and positive density before it replaces the current
state. `project_b=true` projects the final magnetic field onto its divergence-free
subspace.
"""
function step_hall_mhd!(st::HallMHDState{D,T}, dt::Real; project_b::Bool = true) where {D,T}
    dtT = _validated_nonnegative_dt(T, dt; name = "step_hall_mhd!")
    iszero(dtT) && return st

    _copy_hall_state!(st.n0, st.u0, st.B0, st.fields.n, st.fields.ui, st.fields.B)

    _hall_mhd_rhs!(st.kn[1], st.ku[1], st.kB[1], st, st.n0, st.u0, st.B0)
    _hall_stage!(
        st.ntmp,
        st.utmp,
        st.Btmp,
        st.n0,
        st.u0,
        st.B0,
        st.kn[1],
        st.ku[1],
        st.kB[1],
        dtT / 2,
    )
    _hall_mhd_rhs!(st.kn[2], st.ku[2], st.kB[2], st, st.ntmp, st.utmp, st.Btmp)
    _hall_stage!(
        st.ntmp,
        st.utmp,
        st.Btmp,
        st.n0,
        st.u0,
        st.B0,
        st.kn[2],
        st.ku[2],
        st.kB[2],
        dtT / 2,
    )
    _hall_mhd_rhs!(st.kn[3], st.ku[3], st.kB[3], st, st.ntmp, st.utmp, st.Btmp)
    _hall_stage!(st.ntmp, st.utmp, st.Btmp, st.n0, st.u0, st.B0, st.kn[3], st.ku[3], st.kB[3], dtT)
    _hall_mhd_rhs!(st.kn[4], st.ku[4], st.kB[4], st, st.ntmp, st.utmp, st.Btmp)

    _hall_final_stage!(st.ntmp, st.utmp, st.Btmp, st, dtT)
    _validate_hall_candidate!(st, st.ntmp, st.utmp, st.Btmp)
    _copy_hall_state!(st.fields.n, st.fields.ui, st.fields.B, st.ntmp, st.utmp, st.Btmp)
    project_b && project_b!(st.fields, st.g)
    hall_mhd_ohms_law!(st)
    st.time[] += dtT
    st.step[] += 1
    return st
end
