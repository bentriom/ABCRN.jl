
using MarkovProcesses

load_model("SIR")
load_model("ER")

new_SIR = deepcopy(SIR)
sync_SIR = new_SIR * create_automaton_F(new_SIR, 0.0, Inf, 100.0, 110.0, :I)
new_ER = deepcopy(ER)
sync_ER = new_ER * create_automaton_F(new_ER, 0.0, 100.0, 4.0, 5.0, :P)

simulate(sync_SIR)
simulate(sync_ER)

return true

