# JET static-analysis smoke test (§21.2): run JET's whole-package abstract
# interpretation over HybridPlasmaPIC to surface obvious type instabilities and
# undefined-variable / dispatch errors at load time.
#
# Package-wide abstract reports can originate in our code or dependencies.
# Print every report for triage and bound the count to catch regressions.
# Concrete-call analysis is needed before classifying a report as a runtime
# defect, an inference limitation, or an upstream issue.
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
        # Preserve the full findings so the bounded count does not hide their cause.
        if !isempty(reports)
            show(stdout, MIME"text/plain"(), rep)
            println()
        end
        # The ceiling bounds package-wide findings, including dependency reports;
        # it is not a claim that every report below the ceiling is harmless.
        JET_REPORT_CEILING = 30
        @test length(reports) <= JET_REPORT_CEILING
    else
        @test_skip "JET not available in this environment"
    end
end
