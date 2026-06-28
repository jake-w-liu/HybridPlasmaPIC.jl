# Aqua package-quality checks (§21.1): no method ambiguities or piracy, no
# undefined exports, no stale/under-constrained dependencies.
#
# unbound_args is disabled: the only flagged method is the FourierGrid(NTuple{D})
# constructor, "unbound" solely for the impossible D=0 (empty-tuple) case — the
# well-known Aqua false positive for NTuple{D} APIs.

using HybridPlasmaPIC, Test

const AQUA_AVAILABLE = try
    @eval import Aqua
    true
catch err
    @warn "Aqua unavailable; skipping quality checks" exception = err
    false
end

@testset "Aqua quality" begin
    if AQUA_AVAILABLE
        import Aqua
        Aqua.test_all(HybridPlasmaPIC; unbound_args = false, persistent_tasks = false)
    else
        @test_skip "Aqua not available in this environment"
    end
end
