
using Test

@testset "Unit tests" begin
    @test include("unit/load_model.jl")
    @test include("unit/load_module.jl")
    @test include("unit/simulate_sir.jl")
    @test include("unit/simulate_sir_bounded.jl")
    @test include("unit/simulate_er.jl")
end

