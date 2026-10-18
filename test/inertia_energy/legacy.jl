# Frozen public owner from 99760f0; canonical migration checks preserve this baseline.
function original_budget(
    ps::H._HybridPopulations{D,T},
    B::NTuple{3,<:Array{T,D}},
    n::Array{T,D},
    closure::ElectronClosure,
    g::FourierGrid{D,T};
    de2::Real = 0.0,
    pe = nothing,
) where {D,T}
    H._require_grid_tuple(:B, B, g)
    H._require_grid_array(:n, n, g)
    de2T = H._require_finite_nonnegative_real("de2", de2, T)
    H._evolves_pressure(closure) &&
        !iszero(de2T) &&
        throw(ArgumentError("EvolvingElectrons energy requires de2 = 0"))
    ek = kinetic_energy(ps)
    em = magnetic_energy(B, g)
    ei = electron_internal_energy(n, B, closure, g; pe)
    ej = zero(T)
    if de2T > zero(T)
        J = ntuple(_ -> similar(B[1]), 3)
        curl!(J, B, g)
        s = zero(T)
        @inbounds for I in eachindex(J[1], J[2], J[3])
            s += J[1][I]^2 + J[2][I]^2 + J[3][I]^2
        end
        ej = T(0.5) * de2T * s * prod(g.dx)
    end
    return (
        kinetic = ek,
        magnetic = em,
        electron_internal = ei,
        electron_inertia = ej,
        total = ek + em + ei + ej,
    )
end
