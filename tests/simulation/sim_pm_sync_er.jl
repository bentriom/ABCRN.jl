
using MarkovProcesses

load_plots()
load_model("ER")
load_automaton("automaton_F")

A_F = create_automaton_F(ER, 0.0, 100.0, 7.0, 8.0, :P)
pm_sync_ER = ParametricModel(ER*A_F, (:k1, Uniform(0.0,100.0)), (:k2, Uniform(0.0,100.0)))

prior_p = [0.2, 40.0]

σ = simulate(pm_sync_ER, prior_p)
plot(σ; filename = get_module_path() * "/tests/simulation/res_pics/sim_pm_sync_er_long.svg")

return true

