# EnergyEquation.jl — electron pressure evolution (from hybrid.jl)

"""
    advance_electron_pressure!(pe, ue, dt, γe, g)

One explicit step of the adiabatic electron-pressure equation
`∂p_e/∂t = −u_e·∇p_e − γ_e p_e ∇·u_e` on a periodic grid (the energy-equation
closure; source terms −∇·q_e, resistive heating are optional and omitted here).
"""
function advance_electron_pressure!(
    pe::Array{T,D},
    ue::NTuple{3,<:Array{T,D}},
    dt,
    γe,
    g::FourierGrid{D,T},
) where {D,T}
    gradpe = ntuple(_ -> similar(pe), D)
    gradient!(gradpe, pe, g)
    divu = similar(pe)
    divergence!(divu, ue, g)
    dtT = T(dt)
    γ = T(γe)
    @inbounds for I in eachindex(pe)
        adv = zero(T)
        for d = 1:D
            adv += ue[d][I] * gradpe[d][I]
        end
        pe[I] += dtT * (-adv - γ * pe[I] * divu[I])
    end
    return pe
end
