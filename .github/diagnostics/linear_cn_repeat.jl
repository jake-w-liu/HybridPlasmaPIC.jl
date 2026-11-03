# Repeat the permanent test_linear_cn_range.jl in one Ubuntu process. Each
# include performs 24 packet record comparisons; 20 includes give 480.
using Test

const REPEATS = 20
results = @testset "linear CN repeat" begin
    for i = 1:REPEATS
        @testset "repeat $i" begin
            include(joinpath(@__DIR__, "..", "..", "test", "test_linear_cn_range.jl"))
        end
    end
end
