
using Test

@testset "Unit tests" begin
    @test include("unit/load_model.jl")
    @test include("unit/load_module.jl")
end

