# Aqua package-quality checks (§21.1): no method ambiguities or piracy, no
# undefined exports, no stale/under-constrained dependencies.
#
# unbound_args is disabled: the only flagged method is the FourierGrid(NTuple{D})
# constructor, "unbound" solely for the impossible D=0 (empty-tuple) case — the
# well-known Aqua false positive for NTuple{D} APIs.

using HybridPlasmaPIC, Test, Aqua

@testset "Aqua quality" begin
    Aqua.test_all(HybridPlasmaPIC; unbound_args = false, persistent_tasks = false)
end
