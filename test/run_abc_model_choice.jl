
using Test

@testset "ABC model choice tests" begin
    @test include("abc_model_choice/toy_example.jl")
    @test include("abc_model_choice/toy_example_ma.jl")
end

