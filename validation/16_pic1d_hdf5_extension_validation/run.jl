#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_16_pic1d_hdf5_extension_validation(artifact_dir::AbstractString)
    id = "16_pic1d_hdf5_extension_validation"
    results = ValidationResult[]

    n = 32
    l = 2π
    g = FourierGrid((n,), (l,))
    es = Electrostatic1D(g, 0; n0 = 1.0)
    amp = 0.2
    mode = 2
    k = 2π * mode / l
    expected_E = [amp / k * sin(k * (i - 1) * g.dx[1]) for i = 1:n]
    for i = 1:n
        es.ne[i] = es.n0 - amp * cos(k * (i - 1) * g.dx[1])
    end
    poisson_E!(es)
    espic_error = maximum(abs, es.E .- expected_E)
    field_energy_error = abs(field_energy(es) - 0.5 * sum(abs2, es.E) * prod(g.dx))

    gem = FourierGrid((16,), (2π,))
    N = 16 * 16
    electrons = ParticleSet{1,Float64}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(electrons, 0.0, 2π)
    set_density_weight!(electrons, 1.0, gem)
    for p = 1:N
        electrons.v[1][p] = 0.01 * sin(electrons.x[1][p])
    end
    em = EMPIC1D(gem, N; n0 = 1.0, c = 4.0, shape = CIC())
    init_empic!(em, electrons)
    for _ = 1:20
        step_empic!(em, electrons, 0.01)
    end
    empic_residual = charge_conservation_residual(em, 0.01)
    empic_finite_error = isfinite(em_field_energy(em) + kinetic_energy(electrons)) ? 0.0 : 1.0

    ions = ParticleSet{1,Float64}(N; q = 1.0, m = 25.0)
    load_lattice_1d!(ions, 0.0, 2π)
    set_density_weight!(ions, 1.0, gem)
    mobile_e = ParticleSet{1,Float64}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(mobile_e, 0.0, 2π)
    set_density_weight!(mobile_e, 1.0, gem)
    mobile = EMPIC1D(gem, N; n0 = 1.0, c = 4.0, shape = CIC(), mobile = true, mi = 25.0, n_sub = 2)
    init_empic!(mobile, mobile_e, ions)
    step_empic!(mobile, mobile_e, ions, 0.02)
    mobile_residual = charge_conservation_residual(mobile, 0.02)

    artifact = joinpath(artifact_dir, "16_pic1d_hdf5_extension_validation.csv")
    rows = Any[
        ("electrostatic1d_poisson_max_abs_error", espic_error, 0.0, "absolute", espic_error, 1e-12),
        (
            "electrostatic1d_field_energy_abs_error",
            field_energy(es),
            0.5 * sum(abs2, es.E) * prod(g.dx),
            "absolute",
            field_energy_error,
            0.0,
        ),
        (
            "empic1d_charge_conservation_residual",
            empic_residual,
            0.0,
            "absolute",
            empic_residual,
            1e-8,
        ),
        (
            "empic1d_finite_energy_error",
            empic_finite_error,
            0.0,
            "absolute",
            empic_finite_error,
            0.0,
        ),
        (
            "empic1d_mobile_subcycle_charge_residual",
            mobile_residual,
            0.0,
            "absolute",
            mobile_residual,
            1e-8,
        ),
    ]

    try
        @eval import HDF5
        h5_error = Base.invokelatest(_hdf5_roundtrip_error)
        push!(
            rows,
            ("hdf5_extension_roundtrip_max_abs_error", h5_error, 0.0, "absolute", h5_error, 0.0),
        )
    catch err
        push!(
            results,
            _skip_result(
                id = id,
                category = "hdf5_extension",
                reference_kind = "external_open_source",
                reference = "HDF5.jl package extension",
                metric = "hdf5_extension_roundtrip_max_abs_error",
                notes = "Skipped because HDF5 extension is not loadable: $(typeof(err))",
            ),
        )
    end
    _write_metric_csv(artifact, rows)
    append!(
        results,
        _metric_rows_to_results(
            id = id,
            category = "pic_and_extensions",
            reference_kind = "analytic_or_roundtrip",
            reference = "1D PIC analytic checks, EMPIC1D charge conservation, and HDF5 dense-array roundtrip",
            rows = rows,
            artifact = artifact,
        ),
    )
    return results
end

function _hdf5_roundtrip_error()
    h5path = tempname() * ".h5"
    try
        A = reshape(collect(Float64, 1:12), 3, 4)
        write_field_hdf5(h5path, "density", A)
        B = read_field_hdf5(h5path, "density")
        return maximum(abs, B .- A)
    finally
        rm(h5path; force = true)
    end
end


VALIDATION_CASE = ValidationCase(
    id = "16_pic1d_hdf5_extension_validation",
    default = true,
    description = "Electrostatic1D, EMPIC1D, mobile/subcycled PIC, and HDF5 extension.",
    runner = case_16_pic1d_hdf5_extension_validation,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
