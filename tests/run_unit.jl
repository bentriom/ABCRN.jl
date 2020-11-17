
using Test

@testset "Unit tests" begin
    @test include("unit/load_model.jl")
    @test include("unit/load_model_bench.jl")
    @test include("unit/load_module.jl")
    @test include("unit/simulate_sir.jl")
    @test include("unit/simulate_sir_bounded.jl")
    @test include("unit/simulate_er.jl")
    @test include("unit/change_obs_var_sir.jl")
    @test include("unit/change_obs_var_sir_2.jl")
    @test include("unit/getindex_access_trajectory.jl")
end

