
using BiochemNetABC

load_model("SIR")
load_model("ER")
load_automaton("automaton_F")
load_automaton("automaton_G")
load_automaton("automaton_G_and_F")

new_SIR = deepcopy(SIR)
new_ER = deepcopy(ER)
observe_all!(new_SIR)
observe_all!(new_ER)

sync_SIR = new_SIR * create_automaton_F(new_SIR, 0.0, 0.0, 100.0, 110.0, :I)
sync_ER = new_ER * create_automaton_F(new_ER, 0.0, 100.0, 4.0, 5.0, :P)
simulate(sync_SIR)
simulate(sync_ER)

sync_SIR = new_SIR * create_automaton_G(new_SIR, 1.0, Inf, 0.0, 100.0, :I)
sync_ER = new_ER * create_automaton_G(new_ER, 50.0, 100.0, 0.0, 5.0, :E)
simulate(sync_SIR)
simulate(sync_ER)

sync_SIR = new_SIR * create_automaton_G_and_F(new_SIR, 1.0, Inf, 0.0, 100.0, :I, 0.0, Inf, 100.0, 110.0, :I)
sync_ER = new_ER * create_automaton_G_and_F(new_ER, 50.0, 100.0, 0.0, 5.0, :E, 40.0, 100.0, 5.0, 6.0, :P)
simulate(sync_SIR)
simulate(sync_ER)

return true

