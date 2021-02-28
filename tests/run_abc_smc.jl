
using Test

@testset "ABC SMC and automaton-ABC tests" begin
    @test include("automaton_abc/R1.jl")
    @test include("automaton_abc/distributed_R1.jl")
    @test include("automaton_abc/abc_euclidean_distance.jl")
    #@test test_distributed_R1
end

