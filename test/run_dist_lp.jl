
using Test

@testset "Distance Lp tests" begin
    @test include("dist_lp/dist_l1_case_1.jl")
    @test include("dist_lp/dist_case_2.jl")
    @test include("dist_lp/dist_case_3.jl")
    @test include("dist_lp/dist_sim_sir.jl")
end

