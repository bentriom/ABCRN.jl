
using Test

@testset "Automata tests" begin
    @test include("automata/read_trajectory_last_state_F.jl")
    @test include("automata/sync_simulate_last_state_F.jl")
end

