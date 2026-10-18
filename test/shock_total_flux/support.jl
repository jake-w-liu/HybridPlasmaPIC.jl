using HybridPlasmaPIC, Test, Random
include("../shock_energy_flux/legacy.jl")
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function rounded(::Type{T}, q::Q) where {T}
    setprecision(BigFloat, 8192) do
        lo = setrounding(BigFloat, RoundDown) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        hi = setrounding(BigFloat, RoundUp) do
            BigFloat(numerator(q)) / BigFloat(denominator(q))
        end
        @test Q(lo) <= q <= Q(hi)
        a, b = T(lo, RoundNearest), T(hi, RoundNearest)
        @test isequal(a, b)
        a
    end
end
function reference(sh, ps)
    ion = fill(zero(Q), sh.s.n)
    if ps !== nothing
        for p in eachindex(ps.weight)
            position = Q(ps.x[1][p]) / Q(sh.s.dx)
            cell = floor(Int, position)
            fraction = position - cell
            a, b = clamp(cell + 1, 1, sh.s.n), clamp(cell + 2, 1, sh.s.n)
            value = Q(ps.weight[p]) * sum(Q(v[p])^2 for v in ps.v) * Q(ps.v[1][p]) / 2
            ion[a] += (1 - fraction) * value
            ion[b] += fraction * value
        end
    end
    ntuple(
        j -> begin
            i = j == 1 ? sh.s.n : 1
            kinetic =
                ps === nothing ?
                Q(sh.n[i]) * (Q(sh.ux[i])^2 + Q(sh.uy[i])^2) * Q(sh.ux[i]) / 2 :
                ion[i] / Q(sh.s.H[i])
            Q(sh.Ey[i]) * Q(sh.Bz[i]) +
            kinetic +
            Q(sh.γe) * Q(sh.pe[i]) * Q(sh.ux[i]) / (Q(sh.γe) - 1)
        end,
        2,
    )
end
old(sh::S, ps::P) where {S,P} = @timed Legacy.boundary_energy_flux(sh; ps = ps)
new(sh::S, ps::P) where {S,P} = @timed H.boundary_energy_flux(sh; ps = ps)
node_old(sh::S, ps::P) where {S,P} = @allocated Legacy._ion_energy_flux(sh, ps)
node_new(sh::S, ps::P) where {S,P} = @allocated H._ion_energy_flux(sh, ps)
