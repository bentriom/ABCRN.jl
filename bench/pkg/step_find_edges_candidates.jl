
using Profile
import Statistics: mean
using BenchmarkTools
using MarkovProcesses

load_model("ER")
load_automaton("automaton_F")
A_F = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P)

EdgeTransition = Union{Nothing,Vector{Transition}}
struct newEdge
    tr::EdgeTransition
    sym_cc::Symbol
    sym_us::Symbol
end
newEdge2 = Tuple{EdgeTransition,Symbol,Symbol}
newEdge3 = Tuple{EdgeTransition,Function,Function}

S0 = init_state(A_F, ER.x0, ER.t0)
Sn = copy(S0)
Sn_values = [100.0, Inf, 0.0]
Sn.values = Sn_values
Sn_time = Sn.time

# nplus1 values
xnplus1 = [99, 99, 1, 0] 
tnplus1 = 0.01
constants = A_F.constants
p_sim = ER.p

edge_candidates = Vector{EdgeAutomatonF}(undef, 2)
edges_from_current_loc = getfield(A_F, :map_edges)[:l0]
Λ = getfield(A_F, :Λ)
nbr_candidates = _find_edge_candidates!(Vector{EdgeAutomatonF}(undef, 2), edges_from_current_loc, Λ, Sn_time, Sn_values, xnplus1, p_sim, false)
Snplus1 = copy(Sn)

edges_from_current_loc = getfield(A_F, :map_edges)[:l1]
b_find = @benchmark _find_edge_candidates!(edge_candidates, edges_from_current_loc, A_F.Λ, Sn_time, Sn_values, xnplus1, p_sim, false)
edge = edge_candidates[1]
b_get_idx = @benchmark _get_edge_index(edge_candidates, nbr_candidates, false, :R1)
println("Find the edge")
@show minimum(b_find), mean(b_find), maximum(b_find)
@show minimum(b_get_idx), mean(b_get_idx), maximum(b_get_idx)

b_cc = @benchmark check_constraints_AutomatonF(edge, Sn_time, Sn_values, xnplus1, p_sim)

#=
edge_type2 = newEdge(edge.transitions, sym_cc_edge, sym_us_edge)
@show typeof(edge_type2)
b_cc2 = @benchmark getfield(Main, getfield(edge_type2, :sym_cc))(Sn_time, Sn_values, xnplus1, p_sim)
edge_type3 = (edge.transitions, sym_cc_edge, sym_us_edge)
@show typeof(edge_type3)
b_cc3 = @benchmark getfield(Main, edge_type3[2])(Sn_time, Sn_values, xnplus1, p_sim)
b_cc3_bis = @benchmark getfield(Main, @inbounds(edge_type3[2]))(Sn_time, Sn_values, xnplus1, p_sim)
edge_type4 = (edge.transitions, edge.check_constraints, edge.update_state!)
@show typeof(edge_type4)
b_cc4 = @benchmark edge_type4[2](Sn_time, Sn_values, xnplus1, p_sim)
b_cc4_bis = @benchmark @inbounds(edge_type4[2])(Sn_time, Sn_values, xnplus1, p_sim)
=#

function mycc(t::Float64, values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64})
    return (t <= 0.025 && (values[1] < 50.0 || values[1] > 75.0))
end
vals_state = Sn.values
time_state = Sn.time
l_sym = Symbol[:mycc]
b_mycc = @benchmark getfield(Main, :mycc)(Sn_time, Sn_values, xnplus1, p_sim)
b_mycc2 = @benchmark getfield(Main, l_sym[1])(Sn_time, Sn_values, xnplus1, p_sim)
#=
new_edge = Edge(nothing, getfield(Main, :mycc), getfield(Main, :mycc))
b_mycc3 = @benchmark getfield(new_edge, :check_constraints)(Sn_time, Sn_values, xnplus1, p_sim)
=#
println("Check constraints")
#@show Symbol(edge.check_constraints)
@show minimum(b_cc), mean(b_cc), maximum(b_cc)
#=
@show minimum(b_cc2), mean(b_cc2), maximum(b_cc2)
@show minimum(b_cc3), mean(b_cc3), maximum(b_cc3)
@show minimum(b_cc3_bis), mean(b_cc3_bis), maximum(b_cc3_bis)
@show minimum(b_cc4), mean(b_cc4), maximum(b_cc4)
@show minimum(b_cc4_bis), mean(b_cc4_bis), maximum(b_cc4_bis)
=#
@show minimum(b_mycc), mean(b_mycc), maximum(b_mycc)
@show minimum(b_mycc2), mean(b_mycc2), maximum(b_mycc2)
#@show minimum(b_mycc3), mean(b_mycc3), maximum(b_mycc3)

println("Update state")
#@show Symbol(edge.update_state!)
b_us = @benchmark update_state_AutomatonF!(edge, Sn_time, Sn_values, xnplus1, p_sim)
@show minimum(b_us), mean(b_us), maximum(b_us)

Profile.clear_malloc_data()
for i = 1:10
    check_constraints_AutomatonF(edge, Sn_time, Sn_values, xnplus1, p_sim)
end

