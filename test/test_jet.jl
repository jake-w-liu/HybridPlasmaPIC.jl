# JET static-analysis smoke test (§21.2): run JET's whole-package abstract
# interpretation over HybridPlasmaPIC to surface obvious type instabilities and
# undefined-variable / dispatch errors at load time.
#
# This is REPORT-ONLY: JET routinely flags benign issues that originate in
# upstream dependencies (FFTW, Base, the standard library) and that are not
# under our control. Failing CI on those would make the test a flaky gate
# rather than a useful signal, so we run the analysis, print its summary, and
# assert only that the analysis itself completed. The printed report is the
# artifact a developer inspects when tightening type stability.
#
# If JET cannot be loaded (e.g. a Julia version where the pinned JET is
# unavailable), the test is skipped rather than erroring.

using HybridPlasmaPIC
using Test

const JET_AVAILABLE = try
    @eval import JET
    true
catch err
    @info "JET unavailable; skipping static analysis" exception = err
    false
end

@testset "JET static analysis" begin
    if JET_AVAILABLE
        # report_package runs JET's optimization/error analysis over every
        # method reachable from the package's bindings. toplevel_logger=nothing
        # silences the progress chatter so the captured CI log stays compact.
        rep = JET.report_package(HybridPlasmaPIC; toplevel_logger = nothing)
        reports = JET.get_reports(rep)
        @info "JET report-only summary" n_reports = length(reports)
        # Surface the findings in the log without gating CI on upstream noise.
        if !isempty(reports)
            show(stdout, MIME"text/plain"(), rep)
            println()
        end
        # The analysis ran to completion — that is what we verify.
        @test true
    else
        @test_skip "JET not available in this environment"
    end
end
