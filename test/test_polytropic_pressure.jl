module PolytropicPressureTests
using HybridPlasmaPIC, Test
include("isothermal_energy/oracles.jl")
include("polytropic_pressure/cold.jl")
include("polytropic_pressure/ordinary.jl")
include("polytropic_pressure/float32.jl")
include("polytropic_pressure/edges.jl")
include("polytropic_pressure/precision.jl")
include("polytropic_pressure/arrays.jl")
include("polytropic_pressure/reuse.jl")
include("polytropic_pressure/shocks.jl")
include("polytropic_pressure/resources.jl")
include("polytropic_pressure/reuse_resources.jl")
include("polytropic_pressure/binomial.jl")
end
