
using BiochemNetABC

load_model("SIR")
load_model("ER")
load_automaton("automaton_F")
load_automaton("automaton_G")

t1, t2, x1, x2 = 100.0, 120.0, 1.0, 100.0
A_F = create_automaton_F(SIR, x1, x2, t1, t2, :I)

t1, t2, x1, x2 = 0.0, 0.8, 50.0, 100.0
A_G = create_automaton_G(ER, x1, x2, t1, t2, :P)

return true

