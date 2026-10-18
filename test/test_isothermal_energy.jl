module IsothermalEnergyTests
using HybridPlasmaPIC, Test, Random
include("isothermal_energy/oracles.jl")
include("isothermal_energy/cold.jl")
include("isothermal_energy/nearunit.jl")
include("isothermal_energy/dyadic.jl")
include("isothermal_energy/general.jl")
include("isothermal_energy/public.jl")
end
