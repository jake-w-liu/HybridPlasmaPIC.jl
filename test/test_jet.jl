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
    # Reports located in HybridPlasmaPIC itself are always a defect: gate them
    # at zero. The total below is a measured dependency/Base baseline.
    located = JET.report_package(
        HybridPlasmaPIC;
        toplevel_logger = nothing,
        target_modules = (HybridPlasmaPIC,),
    )
    located_reports = JET.get_reports(located)
    @info "JET package-located reports" n_reports = length(located_reports)
    if !isempty(located_reports)
        show(stdout, MIME"text/plain"(), located)
        println()
    end
    @test isempty(located_reports)
    # Measured 2026-11-03 with JET 0.13.1 once package-located reports reached
    # zero: 86 on Julia 1.12.7 and 79 on 1.13.3. All remaining reports are
    # located in dependencies or Base, mostly on solver paths through
    # Krylov/MPI/Serialization; the cap is the larger measurement.
    JET_TOTAL_BASELINE = 86
    @test length(reports) <= JET_TOTAL_BASELINE
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
