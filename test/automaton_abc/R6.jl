
using ABCRN

load_model("ER")
observe_all!(ER)
load_automaton("automaton_G_and_F")
A_G_F_R6 = create_automaton_G_and_F(ER, 50.0, 100.0, 0.0, 0.8, :E, 30.0, 100.0, 0.8, 0.9, :P)
sync_ER = A_G_F_R6 * ER
pm_sync_ER = ParametricModel(sync_ER, (:k1, Uniform(0.0,100.0)), (:k2, Uniform(0.0,100.0)))

r = automaton_abc(pm_sync_ER; alpha=0.2)

using Plots
histogram2d(r.mat_p_end[1,:], r.mat_p_end[2,:], bins=50)
savefig("R6_hist.svg")

