
using Test

@testset "Automata tests" begin
    @test include("automata/absorbing_x0_p.jl")
    @test include("automata/accept_R5.jl")
    @test include("automata/euclidean_distance.jl")
    @test include("automata/euclidean_distance_single.jl")
    @test include("automata/read_trajectory_last_state_F.jl")
    @test include("automata/read_trajectory_last_state_G.jl")
    @test include("automata/sync_simulate_last_state_F.jl")
    @test include("automata/sync_simulate_last_state_G.jl")
    @test include("automata/two_automata.jl")
end

