# OhmsLaw.jl — generalized Ohm's law + Faraday + div control (from hybrid.jl)

"""
    ohms_law!(f, model, g)

Fill `f.E` from generalized Ohm's law using `f.n`, `f.ui`, `f.B`. Also updates
`f.J = ∇×B`, `f.pe`, `f.gradp`, `f.ninv`, and the density-floor counter.
"""
ohms_law!(f::HybridFields{D,T}, model::HybridModel, g::FourierGrid{D,T}) where {D,T} = _ohm_E!(
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

# FROZEN-moment part of Ohm's law (depends only on n, NOT B): electron pressure,
# its gradient, and 1/n with the density floor. During a B-subcycle n is frozen,
# so this is computed ONCE per step (via _ohm_prep!) and reused across every
# subcycle RHS evaluation instead of being recomputed each time — the gradient!
# (D inverse FFTs) was the dominant redundant cost. Writes pe, gradp, ninv, fc.
function _ohm_prep!(
    pe::Array{T,D},
    gradp::NTuple{D,<:Array{T,D}},
    ninv::Array{T,D},
    n::Array{T,D},
    closure::ElectronClosure,
    nfloor::T,
    fc::Base.RefValue{Int},
    g::FourierGrid{D,T},
) where {D,T}
    size(n) == g.n || throw(DimensionMismatch("n size $(size(n)) does not match grid size $(g.n)"))
    size(pe) == g.n ||
        throw(DimensionMismatch("pe size $(size(pe)) does not match grid size $(g.n)"))
    size(ninv) == g.n ||
        throw(DimensionMismatch("ninv size $(size(ninv)) does not match grid size $(g.n)"))
    for d = 1:D
        size(gradp[d]) == g.n || throw(
            DimensionMismatch("gradp[$d] size $(size(gradp[d])) does not match grid size $(g.n)"),
        )
    end
    electron_pressure!(pe, n, closure)
    gradient!(gradp, pe, g)                              # ∇p_e (D spatial comps)
    cnt = 0
    @inbounds for I in eachindex(n)
        nv = n[I]
        if nv < nfloor
            cnt += 1
            nv = nfloor
        end
        ninv[I] = one(T) / nv
    end
    fc[] = cnt
    return ninv
end

# B-DEPENDENT part of Ohm's law: E from an explicit (trial) B, using the
# pre-computed `gradp` and `ninv` from _ohm_prep!. Recomputes J = ∇×B (and ∇²J if
# ηH≠0) each call — these genuinely depend on B. Writes E, J, lapJ.
function _ohm_Efield!(
    E::NTuple{3,<:Array{T,D}},
    ui::NTuple{3,<:Array{T,D}},
    B::NTuple{3,<:Array{T,D}},
    J::NTuple{3,<:Array{T,D}},
    lapJ::NTuple{3,<:Array{T,D}},
    gradp::NTuple{D,<:Array{T,D}},
    ninv::Array{T,D},
    η::T,
    ηH::T,
    g::FourierGrid{D,T},
) where {D,T}
    size(ninv) == g.n ||
        throw(DimensionMismatch("ninv size $(size(ninv)) does not match grid size $(g.n)"))
    for c = 1:3
        size(E[c]) == g.n ||
            throw(DimensionMismatch("E[$c] size $(size(E[c])) does not match grid size $(g.n)"))
        size(ui[c]) == g.n ||
            throw(DimensionMismatch("ui[$c] size $(size(ui[c])) does not match grid size $(g.n)"))
        size(lapJ[c]) == g.n || throw(
            DimensionMismatch("lapJ[$c] size $(size(lapJ[c])) does not match grid size $(g.n)"),
        )
    end
    for d = 1:D
        size(gradp[d]) == g.n || throw(
            DimensionMismatch("gradp[$d] size $(size(gradp[d])) does not match grid size $(g.n)"),
        )
    end
    curl!(J, B, g)                                       # J = ∇×B
    if ηH != 0                                           # hyperresistivity ∇²J
        for c = 1:3
            laplacian!(lapJ[c], J[c], g)
        end
    end
    Bx, By, Bz = B
    ux, uy, uz = ui
    Jx, Jy, Jz = J
    Ex, Ey, Ez = E
    gp = gradp
    LJx, LJy, LJz = lapJ
    @inbounds for I in eachindex(ninv)
        bx = Bx[I]
        by = By[I]
        bz = Bz[I]
        vx = ux[I]
        vy = uy[I]
        vz = uz[I]
        jx = Jx[I]
        jy = Jy[I]
        jz = Jz[I]
        inv = ninv[I]
        # −(u_i×B) + (J×B)/n + η J − ηH ∇²J
        e1 = -(vy * bz - vz * by) + (jy * bz - jz * by) * inv + η * jx - ηH * LJx[I]
        e2 = -(vz * bx - vx * bz) + (jz * bx - jx * bz) * inv + η * jy - ηH * LJy[I]
        e3 = -(vx * by - vy * bx) + (jx * by - jy * bx) * inv + η * jz - ηH * LJz[I]
        # −∇p_e/n on the spatial components only
        e1 -= gp[1][I] * inv
        if D >= 2
            e2 -= gp[2][I] * inv
        end
        if D >= 3
            e3 -= gp[3][I] * inv
        end
        Ex[I] = e1
        Ey[I] = e2
        Ez[I] = e3
    end
    return E
end

# Full Ohm's law (prep + field) in one call — used for the carried E (init! and
# end-of-step), which recomputes everything from the full-step density n.
function _ohm_E!(
    E::NTuple{3,<:Array{T,D}},
    n::Array{T,D},
    ui::NTuple{3,<:Array{T,D}},
    B::NTuple{3,<:Array{T,D}},
    J::NTuple{3,<:Array{T,D}},
    lapJ::NTuple{3,<:Array{T,D}},
    pe::Array{T,D},
    gradp::NTuple{D,<:Array{T,D}},
    ninv::Array{T,D},
    closure::ElectronClosure,
    η::T,
    ηH::T,
    nfloor::T,
    fc::Base.RefValue{Int},
    g::FourierGrid{D,T},
) where {D,T}
    _ohm_prep!(pe, gradp, ninv, n, closure, nfloor, fc, g)
    _ohm_Efield!(E, ui, B, J, lapJ, gradp, ninv, η, ηH, g)
    return E
end

# ---------------------------------------------------------------- Faraday + divergence

"""
    faraday_rhs!(dB, E, g)

Compute ∂B/∂t = −∇×E into `dB`.
"""
function faraday_rhs!(
    dB::NTuple{3,<:Array{T,D}},
    E::NTuple{3,<:Array{T,D}},
    g::FourierGrid{D,T},
) where {D,T}
    curl!(dB, E, g)
    for c = 1:3
        dB[c] .*= -one(T)
    end
    return dB
end

"Discrete ∇·B into `out`; returns its L2 norm."
function magnetic_divergence!(
    out::Array{T,D},
    f::HybridFields{D,T},
    g::FourierGrid{D,T},
) where {D,T}
    divergence!(out, f.B, g)
    s = zero(T)
    @inbounds for x in out
        s += x * x
    end
    return sqrt(s)
end

"Project B onto its divergence-free part in place (preserving the k=0 mean field)."
project_b!(f::HybridFields{D,T}, g::FourierGrid{D,T}) where {D,T} = project_divfree!(f.B, g)

# ---------------------------------------------------------------- electron velocity

"""
    electron_velocity!(ue, ui, J, n; nfloor=1e-6)

Electron bulk velocity `u_e = u_i − J/n` (§6.4), with `J = ∇×B` the normalized
current and `n` the quasineutral density. `ue`, `ui`, `J` are 3-tuples of grid
arrays; `n` a grid array. The density floor protects the `1/n` division. `ue` may
alias `ui`. Returns `ue`.
"""
function electron_velocity!(
    ue::NTuple{3,<:Array{T,D}},
    ui::NTuple{3,<:Array{T,D}},
    J::NTuple{3,<:Array{T,D}},
    n::Array{T,D};
    nfloor = 1e-6,
) where {D,T}
    nf = _require_finite_positive_real("nfloor", nfloor, T)
    size(n) == size(ue[1]) ||
        throw(DimensionMismatch("n size $(size(n)) does not match ue size $(size(ue[1]))"))
    @inbounds for c = 1:3
        uec, uic, Jc = ue[c], ui[c], J[c]
        for I in eachindex(uec)
            uec[I] = uic[I] - Jc[I] / max(n[I], nf)
        end
    end
    return ue
end

"""
    electron_velocity!(ue, f::HybridFields, g; nfloor=1e-6)

Compute `J = ∇×B` into `f.J`, then the electron velocity `u_e = u_i − J/n` into
`ue` from the field state `f`. Returns `ue`.
"""
function electron_velocity!(
    ue::NTuple{3,<:Array{T,D}},
    f::HybridFields{D,T},
    g::FourierGrid{D,T};
    nfloor = 1e-6,
) where {D,T}
    curl!(f.J, f.B, g)
    return electron_velocity!(ue, f.ui, f.J, f.n; nfloor)
end

# ---------------------------------------------------------------- multi-species
