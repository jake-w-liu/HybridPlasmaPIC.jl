#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "plot_validation.jl"))

@testset "validation plot case selection" begin
    options = _parse_plot_args(
        ["--artifact-dir", "tmp-validation-artifacts", "--case", "case_b"];
        default_artifact_dir = "unused",
        allow_cases = true,
    )
    @test basename(options.artifact_dir) == "tmp-validation-artifacts"
    @test options.selected == ["case_b"]

    plots = ["case_a" => identity, "case_b" => identity, "case_c" => identity]
    @test _selected_plots(plots, String[]) == plots
    @test _selected_plots(plots, ["case_b"]) == ["case_b" => identity]
    @test_throws ArgumentError _selected_plots(plots, ["case_missing"])
    @test_throws ArgumentError _selected_plots(plots, ["case_b", "case_b"])
    @test_throws ArgumentError _parse_plot_args(
        ["--case", "case_b"];
        default_artifact_dir = "unused",
    )
end

@testset "validation summary plot selection" begin
    mktempdir() do dir
        artifact_dir = joinpath(dir, "artifacts")
        mkpath(artifact_dir)
        open(joinpath(artifact_dir, "validation_summary.csv"), "w") do io
            println(io, "id,category,reference_kind,reference,metric,measured,expected,error_kind,error,tolerance,status,artifact,notes")
            println(io, "case_a,cat,analytic,ref,m_a,0,0,absolute,0.1,1,pass,,")
            println(io, "case_b,cat,analytic,ref,m_b,0,0,absolute,0.2,1,pass,,")
        end

        output = _summary_plot(artifact_dir, ["case_b"])
        @test output == joinpath(artifact_dir, "validation_summary.pdf")
        @test isfile(output)

        label_key = read(joinpath(artifact_dir, "validation_summary_plot_labels.csv"), String)
        @test occursin("case_b", label_key)
        @test !occursin("case_a", label_key)
    end
end
