
using ABCRN

load_model("ER")

σ = simulate(ER)

d1 = Dict(:E => 1, :S => 2, :ES => 3, :P => 4)
d2 = Dict(:P => 1)

bool_test = ER.g == [:P] && ER._g_idx == [4] && 
            ER.map_var_idx == d1 && 
            ER._map_obs_var_idx == d2

return bool_test

