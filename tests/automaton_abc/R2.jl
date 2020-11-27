
using MarkovProcesses
using Plots

load_model("ER")
load_automaton("automaton_F")
A_F_R2 = create_automaton_F(ER, 50.0, 75.0, 0.05, 0.075, "P")
sync_ER = A_F_R2*ER 
pm_sync_ER = ParametricModel(sync_ER, ("k3", Uniform(0.0, 100.0)))

r = automaton_abc(pm_sync_ER; nbr_particles = 1000)
@show r.nbr_sim

histogram(r.mat_p_end', weights = r.weights, normalize = :density)
png("R2_hist.png")

