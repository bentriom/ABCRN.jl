
using BenchmarkTools
using ABCRN

load_model("ER")
load_automaton("automaton_F")
A_F = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P)

S0 = init_state(A_F, ER.x0, ER.t0)
# nplus1 values
xnplus1 = [99, 99, 1, 0] 
tnplus1 = 0.01
tr_nplus1 = :R1
ptr_loc, ptr_time = [S0.loc], [S0.time] 
values = S0.values
edge_candidates = Vector{EdgeAutomatonF}(undef, 2)

b = @benchmark next_state!(A_F, ptr_loc, values, ptr_time, xnplus1, tnplus1, tr_nplus1, ER.x0, ER.p, edge_candidates)
@show minimum(b), mean(b), maximum(b)
@btime next_state!(A_F, ptr_loc, values, ptr_time, xnplus1, tnplus1, tr_nplus1, ER.x0, ER.p, edge_candidates)

