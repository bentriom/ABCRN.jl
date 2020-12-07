
using Test

@testset "Cosmos tests" begin
    test_ER_1D = include("cosmos/distance_F/ER_1D.jl")
    test_ER_R5 = include("cosmos/distance_G/ER_R5.jl")
    test_ER_R6 = include("cosmos/distance_G_F/ER_R6.jl")
    @test test_ER_1D
    @test test_ER_R5
    @test test_ER_R6
end

