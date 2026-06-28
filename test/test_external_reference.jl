using HybridPlasmaPIC, Test

@testset "published external hybrid-code reference metadata" begin
    @test published_hybrid_reference_ids() == (:preisser2020_65deg_Bavg_y,)

    meta = published_hybrid_reference_metadata()
    @test meta.id == :preisser2020_65deg_Bavg_y
    @test meta.doi == "10.5281/zenodo.3697360"
    @test meta.doi_url == "https://doi.org/10.5281/zenodo.3697360"
    @test meta.license == "CC-BY-4.0"
    @test meta.file == "Fig2_65deg_1perc_5perc_10perc_Bavg_y.h5"
    @test meta.file_checksum == "md5:2ea4f239d7221cd52607705f383161a1"
    @test occursin("2D local hybrid simulations", meta.notes)
    @test_throws ArgumentError published_hybrid_reference_metadata(:unknown_reference)
end

@testset "published external hybrid-code Bavg_y summaries" begin
    r1 = published_hybrid_reference(; alpha_fraction = 0.01)
    @test r1.dataset == "65deg_1perc_Bavg_y"
    @test r1.thetaBn_deg == 65.0
    @test r1.alpha_fraction == 0.01
    @test r1.nsamples == 1000
    @test r1.Bavg_y_min == 0.9968850596311658
    @test r1.Bavg_y_max == 3.698209909020219
    @test r1.Bavg_y_mean == 1.924506582194626

    r5 = published_hybrid_reference(; alpha_fraction = 0.05)
    @test r5.dataset == "65deg_5perc_Bavg_y"
    @test r5.Bavg_y_min == 0.997214995499391
    @test r5.Bavg_y_max == 3.500705042674052
    @test r5.Bavg_y_mean == 1.9245920883357377

    r10 = published_hybrid_reference(; alpha_fraction = 0.10)
    @test r10.dataset == "65deg_10perc_Bavg_y"
    @test r10.Bavg_y_min == 0.9969985157435545
    @test r10.Bavg_y_max == 3.5037309774874927
    @test r10.Bavg_y_mean == 1.9378039653740708

    @test_throws ArgumentError published_hybrid_reference(; id = :unknown_reference)
    @test_throws ArgumentError published_hybrid_reference(; alpha_fraction = 0.02)
    @test_throws ArgumentError published_hybrid_reference(; alpha_fraction = NaN)
end

@testset "published external hybrid-code comparison" begin
    reference = published_hybrid_reference(; alpha_fraction = 0.05)
    pass = compare_to_published_hybrid_reference(reference; alpha_fraction = 0.05, rtol = 0.0)
    @test pass.pass
    @test pass.metadata.doi == "10.5281/zenodo.3697360"
    @test pass.comparison.maxrelerr == 0.0

    perturbed = merge(reference, (; Bavg_y_mean = reference.Bavg_y_mean * 1.01))
    @test compare_to_published_hybrid_reference(perturbed; alpha_fraction = 0.05, rtol = 0.02).pass
    @test !compare_to_published_hybrid_reference(
        perturbed;
        alpha_fraction = 0.05,
        rtol = 0.001,
    ).pass

    missing = (; Bavg_y_mean = reference.Bavg_y_mean)
    failed = compare_to_published_hybrid_reference(missing; alpha_fraction = 0.05)
    @test !failed.pass
    @test any(d -> d[1] == :thetaBn_deg && !d[5], failed.comparison.details)
end
