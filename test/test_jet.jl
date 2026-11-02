# JET static-analysis smoke test (§21.2): run JET's whole-package abstract
# interpretation over HybridPlasmaPIC to surface obvious type instabilities and
# undefined-variable / dispatch errors at load time.
#
# Package-wide abstract reports can originate in our code or dependencies.
# Print every report for triage and bound the count to catch regressions.
# Concrete-call analysis is needed before classifying a report as a runtime
# defect, an inference limitation, or an upstream issue.
#
# JET is a required test dependency on supported Julia versions. Import and
# analysis failures must fail this gate instead of silently skipping it.

using HybridPlasmaPIC
using Test
import JET

@testset "JET static analysis" begin
    # report_package runs JET's optimization/error analysis over every
    # method reachable from the package's bindings. toplevel_logger=nothing
    # silences the progress chatter so the captured CI log stays compact.
    rep = JET.report_package(HybridPlasmaPIC; toplevel_logger = nothing)
    reports = JET.get_reports(rep)
    @info "JET package analysis summary" n_reports = length(reports)
    # Preserve the full findings so the bounded count does not hide their cause.
    if !isempty(reports)
        show(stdout, MIME"text/plain"(), rep)
        println()
    end
    # The ceiling bounds package-wide findings, including dependency reports;
    # it is not a claim that every report below the ceiling is harmless.
    JET_REPORT_CEILING = 30
    @test length(reports) <= JET_REPORT_CEILING
    # Retain the validated byte-vector payload type through nested decoding,
    # including when optional packages add unrelated read methods.
    for path_type in (String, SubString{String})
        loader = JET.report_call(load_run, Tuple{path_type})
        if !isempty(JET.get_reports(loader))
            show(stdout, MIME"text/plain"(), loader)
            println()
        end
        @test isempty(JET.get_reports(loader))
    end
end

@testset "Specialized EMPIC1D collection concrete compiler calls" begin
    for T in (Float32, Float64)
        ps = ParticleSet{1,T}(2; q = 0)
        ps.x[1] .= T[0.25, 0.75]
        species = [ps]
        state = EMPIC1D(FourierGrid((4,), (one(T),)), species)
        set_empic1d_schedule!(state; n_sub = 2, species_stride = [2])
        init_empic!(state, species)
        saved = HybridPlasmaPIC._checkpoint_state(state, species)
        calls = (
            (set_empic1d_schedule!, Tuple{typeof(state)}),
            (init_empic!, Tuple{typeof(state),typeof(species)}),
            (step_empic!, Tuple{typeof(state),typeof(species),T}),
            (electromagnetic_diagnostics, Tuple{typeof(state),typeof(species)}),
            (HybridPlasmaPIC._checkpoint_state, Tuple{typeof(state),typeof(species)}),
            (
                HybridPlasmaPIC._restore_checkpoint_state!,
                Tuple{typeof(state),typeof(species),typeof(saved)},
            ),
        )
        for (f, types) in calls
            report = JET.report_call(f, types)
            if !isempty(JET.get_reports(report))
                show(stdout, MIME"text/plain"(), report)
                println()
            end
            @test isempty(JET.get_reports(report))
        end
    end
end
