
using MarkovProcesses

load_model("SIR")
σ = simulate(SIR)
test_1 = length_obs_var(σ) == 1

set_observed_var!(SIR, [:I, :R])
σ = simulate(SIR)
test_2 = length_obs_var(σ) == 2

load_model("ER")
σ = simulate(ER)
test_3 = length_obs_var(σ) == 1

set_observed_var!(ER, [:P, :S, :ES])
σ = simulate(ER)
test_4 = length_obs_var(σ) == 3

test = test_1 && test_2 && test_3 && test_4

return test

