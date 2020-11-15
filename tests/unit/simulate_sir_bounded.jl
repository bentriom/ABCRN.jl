
using MarkovProcesses

load_model("SIR")
SIR.time_bound = 100.0

σ = simulate(SIR)

d1 = Dict("S" => 1, "I"=> 2, "R" => 3)
d2 = Dict("I" => 1)

bool_test = SIR.g == ["I"] && SIR._g_idx == [2] && 
            SIR.map_var_idx == d1 && 
            SIR._map_obs_var_idx == d2 && σ["times"][end] <= SIR.time_bound

return bool_test

