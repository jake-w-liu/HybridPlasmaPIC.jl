# Arithmetic compatibility oracle from 7104c61; production has one owner.
module Legacy
using HybridPlasmaPIC
function boundary_energy_flux(
    sh::PerpShock{T};
    ps::Union{Nothing,ParticleSet{1,T}} = nothing,
) where {T}
    N = sh.s.n
    γe = sh.γe
    # electron enthalpy flux  h_e = γe/(γe−1) · p_e · u_x  (internal energy + p·u work)
    enth(i) = (γe / (γe - one(T))) * sh.pe[i] * sh.ux[i]
    # magnetic / Poynting flux S_x = E_y B_z at each boundary
    mag_inflow = sh.Ey[N] * sh.Bz[N]
    mag_wall = sh.Ey[1] * sh.Bz[1]
    enth_inflow = enth(N)
    enth_wall = enth(1)
    if ps === nothing
        # bulk-only ion kinetic flux ½ n |u|² u_x (ion thermal flux needs the ions)
        bulk(i) = T(0.5) * sh.n[i] * (sh.ux[i]^2 + sh.uy[i]^2) * sh.ux[i]
        kin_inflow = bulk(N)
        kin_wall = bulk(1)
    else
        # full kinetic-ion energy flux Σ_p w ½|v|² v_x (bulk + thermal) at the boundary
        Fe = _ion_energy_flux(sh, ps)
        kin_inflow = Fe[N]
        kin_wall = Fe[1]
    end
    magnetic = (mag_inflow, mag_wall)
    kinetic = (kin_inflow, kin_wall)
    enthalpy = (enth_inflow, enth_wall)
    total = (mag_inflow + kin_inflow + enth_inflow, mag_wall + kin_wall + enth_wall)
    return (; magnetic, kinetic, enthalpy, total)
end

# Full kinetic-ion energy-flux density Σ_p w·½|v|²·v_x deposited on the SBP grid
# (the same CIC stencil + H-norm as deposit_moments!), returned as a length-N
# vector. CIC's partition of unity makes Σ_i H[i]·F[i] = Σ_p w·½|v|²·v_x exactly.
function _ion_energy_flux(sh::PerpShock{T}, ps::ParticleSet{1,T}) where {T}
    N = sh.s.n
    dx = sh.s.dx
    acc = zeros(T, N)
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    w = ps.weight
    HybridPlasmaPIC._require_finite_real_sequence("particle positions", xp)
    HybridPlasmaPIC._require_finite_real_sequence("particle velocities vx", vx)
    HybridPlasmaPIC._require_finite_real_sequence("particle velocities vy", vy)
    HybridPlasmaPIC._require_finite_real_sequence("particle velocities vz", vz)
    HybridPlasmaPIC._require_finite_real_sequence("particle weights", w)
    @inbounds for p in eachindex(w)
        a, b, wa, wb = HybridPlasmaPIC._cic_sbp(xp[p], dx, N)
        keflux = T(0.5) * (vx[p]^2 + vy[p]^2 + vz[p]^2) * vx[p] * w[p]
        acc[a] += wa * keflux
        acc[b] += wb * keflux
    end
    @inbounds for i = 1:N
        acc[i] /= sh.s.H[i]
    end
    return acc
end


end
