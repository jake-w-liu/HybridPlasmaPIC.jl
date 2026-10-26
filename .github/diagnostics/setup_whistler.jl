using Pkg, TOML

length(ARGS) == 2 || error("Expected source checkout and disposable environment paths")
source, environment = abspath.(ARGS)
project = TOML.parsefile(joinpath(source, "Project.toml"))
dependencies = merge(project["deps"], project["extras"])
dependencies[project["name"]] = project["uuid"]
compat = Dict(k => v for (k, v) in project["compat"] if k == "julia" || haskey(dependencies, k))
# Both failed CI jobs resolved this registered provider. Preserve that boundary.
compat["SpectralOperators"] = "=0.1.6"
mkpath(environment)
open(joinpath(environment, "Project.toml"), "w") do io
    TOML.print(
        io,
        Dict(
            "deps" => dependencies,
            "compat" => compat,
            "sources" => Dict(project["name"] => Dict("path" => source)),
        ),
    )
end
Pkg.activate(environment)
Pkg.instantiate()
Pkg.status(; mode = Pkg.PKGMODE_MANIFEST)
