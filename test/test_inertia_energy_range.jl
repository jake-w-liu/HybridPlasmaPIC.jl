for file in ("range.jl", "public.jl", "properties.jl", "quality.jl", "zero_resources.jl")
    include(joinpath("inertia_energy", file))
end
