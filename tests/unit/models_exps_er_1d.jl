
using MarkovProcesses

load_model("ER")
observe_all!(ER)
load_automaton("automaton_F")
load_automaton("automaton_G")

A_F_R1 = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, "P") 
A_F_R2 = create_automaton_F(ER, 50.0, 75.0, 0.05, 0.075, "P") 
A_F_R3 = create_automaton_F(ER, 25.0, 50.0, 0.05, 0.075, "P")
A_G_R5 = create_automaton_G(ER, 50.0, 100.0, 0.0, 0.8, "E") 

return true

