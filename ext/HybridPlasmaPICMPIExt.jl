module HybridPlasmaPICMPIExt

import HybridPlasmaPIC
import MPI

HybridPlasmaPIC.extension_dependency_module(::Val{:mpi}) = MPI

end # module
