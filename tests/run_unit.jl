
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
    @test include("unit/is_always_bounded_sir.jl")
    @test include("unit/length_obs_var.jl")
    @test include("unit/dist_lp_var.jl")
    @test include("unit/dist_lp.jl")
    @test include("unit/l_dist_lp.jl")
    @test include("unit/check_trajectory_consistency.jl")
    @test include("unit/check_model_consistency.jl")
    @test include("unit/set_param.jl")
    @test include("unit/side_effects_1.jl")
end

