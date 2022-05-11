
using Test
import MarkovProcesses: get_module_path

path_bench  = get_module_path() * "/bench/"
function include_arg!(path_file::String, arg::String)
    global ARGS = [arg]
    include(path_bench * "array_memory_order/sim.jl")
    return true
end

@testset "Benchmark" begin
    @test include_arg!(path_bench * "array_memory_order/sim.jl", "SIR")
    @test include_arg!(path_bench * "pygmalion/sim_sir.jl", "SIR")
end

