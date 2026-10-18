module PolytropicEnergyChecks
using Test, HybridPlasmaPIC
@testset "Complete polytropic energy" begin
    include("polytropic_energy/conversion.jl")
    include("polytropic_energy/strict.jl")
    include("polytropic_energy/equality.jl")
    include("polytropic_energy/cold.jl")
    include("polytropic_energy/ordinary.jl")
    include("polytropic_energy/gamma_one.jl")
    include("polytropic_energy/fixed_log.jl")
    include("polytropic_energy/resources.jl")
    include("polytropic_energy/public.jl")
end
end
