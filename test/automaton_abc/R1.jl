
using BiochemNetABC

load_model("ER")
load_automaton("automaton_F")
A_F_R1 = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P)
sync_ER = A_F_R1*ER 
pm_sync_ER = ParametricModel(sync_ER, (:k3, Uniform(0.0, 100.0)))
nbr_pa = 502

r = automaton_abc(pm_sync_ER; nbr_particles = nbr_pa)

test = size(r.mat_p_end)[1] == pm_sync_ER.df &&
       size(r.mat_p_end)[2] == nbr_pa &&
       length(r.vec_dist) == nbr_pa &&
       length(r.weights) == nbr_pa &&
       r.epsilon == 0.0

return test

#using Plots
#histogram(r.mat_p_end', weights = r.weights, normalize = :density)
#savefig("R1_hist.svg")

