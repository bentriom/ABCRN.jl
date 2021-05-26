
using MarkovProcesses

load_model("ER")
observe_all!(ER)
load_automaton("automaton_F")
load_automaton("automaton_G")
load_automaton("automaton_G_and_F")

A_F_R4 = create_automaton_F(ER, 5.0, 15.0, 8.0, 10.0, :P)
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
A_G_R5 = create_automaton_G(ER, x1, x2, t1, t2, :E) 
x3, x4, t3, t4 = 30.0, 100.0, 0.8, 0.9
A_G_F_R6 = create_automaton_G_and_F(ER, x1, x2, t1, t2, :E,
                                    x3, x4, t3, t4, :P)  

return true

