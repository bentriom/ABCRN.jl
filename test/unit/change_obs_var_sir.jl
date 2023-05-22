
using ABCRN

load_model("SIR")

σ = simulate(SIR)
set_observed_var!(SIR, [:I, :R])

d1 = Dict(:S => 1, :I => 2, :R => 3)
d2 = Dict(:I => 1, :R => 2)

bool_test = SIR.g == [:I, :R] && SIR._g_idx == [2,3] && 
            SIR.map_var_idx == d1 && 
            SIR._map_obs_var_idx == d2

return bool_test

