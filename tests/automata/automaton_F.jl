
using MarkovProcesses

load_model("SIR")
load_automaton("automaton_F")
SIR.time_bound = 120.0
x1, x2, t1, t2 = 0.0, Inf, 100.0, 120.0 
A_F = create_automaton_F(SIR, x1, x2, t1, t2, "I") # <: LHA

#sync_SIR = SIR * A_F
#σ, state_lha = simulate(sync_SIR)  

