
using MarkovProcesses

load_model("ER")
load_automaton("automaton_F")
A_F = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, "P")
sync_ER = A_F*ER 
pm_sync_ER = ParametricModel(sync_ER, ("k3", Uniform(0.0, 100.0)))

automaton_abc(pm_sync_ER)

