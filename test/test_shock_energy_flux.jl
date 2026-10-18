# Complete boundary-flux arithmetic, ordinary compatibility and ownership.
for file in ("range.jl", "properties.jl", "bounds.jl", "resources.jl", "contracts.jl")
    include(joinpath("shock_energy_flux", file))
end
