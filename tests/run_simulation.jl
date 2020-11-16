
using Test
import MarkovProcesses: get_module_path

str_dir_pics = get_module_path() * "/tests/simulation/res_pics"
if !isdir(str_dir_pics) mkdir(str_dir_pics) end

@testset "Simulation tests" begin
    @test include("simulation/sim_sir.jl")
    @test include("simulation/sim_sir_bounded.jl")
    @test include("simulation/sim_sir_col_buffer_bounded.jl")
    @test include("simulation/sim_sir_row_buffer_bounded.jl")
    @test include("simulation/sim_er.jl")
end

