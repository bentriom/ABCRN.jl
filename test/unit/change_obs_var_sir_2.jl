
using MarkovProcesses

load_model("SIR")

σ = simulate(SIR)
set_observed_var!(SIR, [:R, :S])

d1 = Dict(:S => 1, :I => 2, :R => 3)
d2 = Dict(:R => 1, :S => 2)

bool_test = SIR.g == [:R, :S] && SIR._g_idx == [3,1] && 
            SIR.map_var_idx == d1 && 
            SIR._map_obs_var_idx == d2

return bool_test

