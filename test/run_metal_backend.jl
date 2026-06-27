#!/usr/bin/env julia

# Run the optional Metal backend tests in a temporary project.
#
# Metal is a weak dependency, so the default `Pkg.test()` environment does not
# install it. This runner keeps the ordinary CPU/CI test target portable while
# giving macOS workstations a reproducible command for extension coverage:
#
#   julia test/run_metal_backend.jl

import Pkg

const TEST_DIR = @__DIR__
const PKG_DIR = dirname(TEST_DIR)

tmp_project = mktempdir(; prefix = "hybridplasmapic-metal-test-")

try
    Pkg.activate(tmp_project)
    Pkg.add(Pkg.PackageSpec(url = "https://github.com/jake-w-liu/SpectralOperators.jl.git"))
    Pkg.develop(Pkg.PackageSpec(path = PKG_DIR))
    Pkg.add(Pkg.PackageSpec(name = "Metal"))
    Pkg.instantiate()

    include(joinpath(TEST_DIR, "test_gpu_backend.jl"))
finally
    rm(tmp_project; force = true, recursive = true)
end
