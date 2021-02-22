
import Statistics: mean
using BenchmarkTools
using MarkovProcesses

load_model("ER")
load_automaton("automaton_F")
A_F = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P)
p_sim = ER.p

S0 = init_state(A_F, ER.x0, ER.t0)
Sn = copy(S0)

# nplus1 values
xnplus1 = [99, 99, 1, 0] 
tnplus1 = 0.01
constants = A_F.constants

Sn[:n] = 100.0
Sn[:d] = Inf
Sn[:isabs] = 0.0
Sn_time = Sn.time
Sn_values = Sn.values
Snplus1 = copy(Sn)
edge_candidates = Vector{Edge}(undef, 2)
edges_from_current_loc = getfield(A_F, :map_edges)[:l0]
nbr_candidates = MarkovProcesses._find_edge_candidates!(edge_candidates, edges_from_current_loc, A_F.Λ, Sn_time, Sn_values, xnplus1, p_sim, false)

edges_from_current_loc = getfield(A_F, :map_edges)[:l1]

function _find_edge_candidates2!(edge_candidates::Vector{Edge},
                                edges_from_current_loc::Dict{Location,Vector{Edge}},
                                Λ::Dict{Location,Function},
                                S_time::Float64, S_values::Vector{Float64},
                                x::Vector{Int}, p::Vector{Float64},
                                only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        if !Λ[target_loc](x) continue end
        for edge in edges_from_current_loc[target_loc]
            cc_func = getfield(edge, :check_constraints)
            if cc_func(S_time, S_values, x, p)
                if getfield(edge, :transitions) == nothing
                    MarkovProcesses._push_edge!(edge_candidates, edge, nbr_candidates)
                    nbr_candidates += 1
                    return nbr_candidates
                else
                    if !only_asynchronous
                        MarkovProcesses._push_edge!(edge_candidates, edge, nbr_candidates)
                        nbr_candidates += 1
                    end
                end
            end
        end
    end
    return nbr_candidates
end

function one_loop_edges(edges_from_current_loc::Dict{Location,Vector{Edge}}, loc::Symbol)
    for edge in edges_from_current_loc[loc]
        cc_func = getfield(edge, :check_constraints)
        tr = getfield(edge, :transitions)
    end
end

function two_loops_edges(edges_from_current_loc::Dict{Location,Vector{Edge}}, 
                         S_time::Float64, S_values::Vector{Float64},
                         x::Vector{Int}, p::Vector{Float64})
    for target_loc in keys(edges_from_current_loc)
        for edge in edges_from_current_loc[target_loc] 
            tr = getfield(edge, :transitions)
            cc_func = getfield(edge, :check_constraints)
            #(S_time, S_values, x, p)
            #    nothing
            #end
        end
    end
end

check_constraints_edge(edge::Edge, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
getfield(edge, :check_constraints)(S_time, S_values, x, p)

function two_loops_edges_full(edge_candidates::Vector{Edge},
                                edges_from_current_loc::Dict{Location,Vector{Edge}},
                                Λ::Dict{Location,Function},
                                S_time::Float64, S_values::Vector{Float64},
                                x::Vector{Int}, p::Vector{Float64},
                                only_asynchronous::Bool)
    for target_loc in keys(edges_from_current_loc)
        for edge in edges_from_current_loc[target_loc] 
            tr = getfield(edge, :transitions)
            if getfield(edge, :check_constraints)(S_time, S_values, x, p)
                nothing
            end
            #ischeck = check_constraints_edge(edge, S_time, S_values, x, p)
            #    nothing
            #end
        end
    end

end

function _find_edge_candidates3!(edge_candidates::Vector{Edge},
                                 edges_from_current_loc::Dict{Location,Vector{Edge}},
                                 Λ::Dict{Location,Function},
                                 Sn_time::Float64, Sn_values::Vector{Float64},
                                 x::Vector{Int}, p::Vector{Float64},
                                 only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        cc_funcs = (e -> getfield(e, :check_constraints)).(edges_from_current_loc[target_loc])
        l_tr = (e -> getfield(e, :transitions)).(edges_from_current_loc[target_loc])
        for i = eachindex(cc_funcs)
            if Λ[target_loc](x) && @inbounds(cc_funcs[i](Sn_time, Sn_values, x, p))
                if l_tr[i] == nothing
                    MarkovProcesses._push_edge!(edge_candidates, edges_from_current_loc[target_loc][i], nbr_candidates)
                    nbr_candidates += 1
                    return nbr_candidates
                else
                    if !only_asynchronous
                        MarkovProcesses._push_edge!(edge_candidates, edges_from_current_loc[target_loc][i], nbr_candidates)
                        nbr_candidates += 1
                    end
                end
            end
        end
    end
    return nbr_candidates
end

edge_candidates = Vector{Edge}(undef, 2)
Λ = A_F.Λ
two_loops_edges_full(edge_candidates, edges_from_current_loc, Λ, Sn_time, Sn_values, xnplus1, p_sim, false)


b_find = @benchmark MarkovProcesses._find_edge_candidates!($(Vector{Edge}(undef, 2)), edges_from_current_loc, Λ, Sn_time, Sn_values, xnplus1, p_sim, false)
@show minimum(b_find), mean(b_find), maximum(b_find)

b_find2 = @benchmark _find_edge_candidates2!($(Vector{Edge}(undef, 2)), edges_from_current_loc, Λ, Sn_time, Sn_values, xnplus1, p_sim, false)
@show minimum(b_find2), mean(b_find2), maximum(b_find2)

b_find3 = @benchmark _find_edge_candidates3!($(Vector{Edge}(undef, 2)), edges_from_current_loc, Λ, Sn_time, Sn_values, xnplus1, p_sim, false)
@show minimum(b_find3), mean(b_find3), maximum(b_find3)

