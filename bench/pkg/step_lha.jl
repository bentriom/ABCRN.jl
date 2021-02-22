
using BenchmarkTools
using MarkovProcesses

load_model("ER")
p_sim = ER.p
load_automaton("automaton_F")
A_F = create_automaton_F(ER, 50.0, 75.0, 0.05, 0.075, :P)
S0 = init_state(A_F, ER.x0, 0.0)
xn = ER.x0
tn = ER.t0
Snplus1 = copy(S0)
xnplus1 = zeros(Int, ER.dim_state)
l_t = Float64[0.0]
l_tr = Transition[nothing]
getfield(Main, ER.f!)(xnplus1, l_t, l_tr, xn, tn, ER.p)
tnplus1 = l_t[1]
tr_nplus1 = l_tr[1]
ptr_loc, ptr_time = [S0.loc], [S0.time] 
values = S0.values
edge_candidates = Vector{EdgeAutomatonF}(undef, 2)

b_ns = @benchmark next_state!($(A_F), $(ptr_loc), $(values), $(ptr_time), $(xnplus1), $(tnplus1), $(Meta.quot(tr_nplus1)), $(xn), $(ER.p), edge_candidates)
@show minimum(b_ns), mean(b_ns), maximum(b_ns)
@btime next_state!($(A_F), $(ptr_loc), $(values), $(ptr_time), $(xnplus1), $(tnplus1), $(Meta.quot(tr_nplus1)), $(xn), $(ER.p), edge_candidates)

#next_state!(Snplus1, A_F, xnplus1, tnplus1, tr_nplus1, S0, xn, p_sim)

