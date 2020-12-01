
using Test

@testset "Cosmos tests" begin
    @test include("cosmos/distance_F/ER_1D.jl")
    @test include("cosmos/distance_G/ER_R5.jl")
end

