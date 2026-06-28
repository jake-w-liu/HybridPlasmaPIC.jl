#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_04_particle_push_collision_boundaries(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "04_particle_push_collision_boundaries.csv",
        "04_particle_push_collision_boundaries.pdf";
        title = "Particle push, boundary, and collision validation",
    )
end

VALIDATION_PLOT = plot_04_particle_push_collision_boundaries

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
