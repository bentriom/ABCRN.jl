
using StaticArrays
using BenchmarkTools
using ABCRN
using Profile

load_model("ER")
observe_all!(ER)
set_param!(ER, [:k1, :k2], [0.2, 40.0])
set_time_bound!(ER, 0.9)

load_automaton("automaton_G_and_F")
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
x3, x4, t3, t4 = 30.0, 100.0, 0.8, 0.9
A_G_F = create_automaton_G_and_F(ER, x1, x2, t1, t2, :E,
                                 x3, x4, t3, t4, :P)
sync_ER = ER * A_G_F

S = init_state(A_G_F, ER.x0, ER.t0)
time = 0.1
values = S.values
x, p = ER.x0, ER.p
for loc_from in keys(A_G_F.map_edges)
    for loc_to in keys(A_G_F.map_edges[loc_from])
        println("$loc_from => $loc_to")
        for i = eachindex(A_G_F.map_edges[loc_from][loc_to])
            println("edge $i:")
            global global_edge = A_G_F.map_edges[loc_from][loc_to][i]
            #@btime check_constraints_AutomatonGandF(global_edge, time, $(copy(values)), x, p)
            #b = @benchmark check_constraints_AutomatonGandF(global_edge, time, $(copy(values)), x, p)
            #@show mean(b).time, mean(b).memory
        end
    end
end

function my_find_edge_candidates!(edge_candidates::Vector{EdgeAutomatonGandF},
                                  edges_from_current_loc::Dict{Location,Vector{EdgeAutomatonGandF}},
                                  Λ::Dict{Location,Function},
                                  S_time::Float64, S_values::Vector{Float64},
                                  x::Vector{Int}, p::Vector{Float64},
                                  only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        if !Λ[target_loc](x) continue end
        for i = eachindex(edges_from_current_loc[target_loc])
            edge = edges_from_current_loc[target_loc][i]
            test_cc = check_constraints_AutomatonGandF(edge, S_time, S_values, x, p)
            if test_cc
                if edge.transitions == nothing
                    _push_edge!(edge_candidates, edge, nbr_candidates)
                    nbr_candidates += 1
                    return nbr_candidates
                else
                    if !only_asynchronous
                        _push_edge!(edge_candidates, edge, nbr_candidates)
                        nbr_candidates += 1
                    end
                end
            end
        end
    end
    return nbr_candidates
end

function svector_find_edge_candidates!(edge_candidates::Vector{EdgeAutomatonGandF},
                                  edges_from_current_loc::Dict,
                                  Λ::Dict{Location,Function},
                                  S_time::Float64, S_values::Vector{Float64},
                                  x::Vector{Int}, p::Vector{Float64},
                                  only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        if !Λ[target_loc](x) continue end
        for edge in edges_from_current_loc[target_loc]
            test_cc = check_constraints_AutomatonGandF(edge, S_time, S_values, x, p)
            if test_cc
                if edge.transitions == nothing
                    _push_edge!(edge_candidates, edge, nbr_candidates)
                    nbr_candidates += 1
                    return nbr_candidates
                else
                    if !only_asynchronous
                        _push_edge!(edge_candidates, edge, nbr_candidates)
                        nbr_candidates += 1
                    end
                end
            end
        end
    end
    return nbr_candidates
end


edge_candidates = Vector{EdgeAutomatonGandF}(undef, 2)
edges_from_current_loc = getfield(A_G_F, :map_edges)[:l1G]
Λ = A_G_F.Λ

function transform_dict_tuple(dict::Dict{Symbol,Dict{Symbol,Vector{EdgeAutomatonGandF}}})
    new_dict = Dict()
    for from_loc in keys(dict)
        for to_loc in keys(dict[from_loc])
            new_dict[from_loc] = (to_loc, dict[from_loc][to_loc])
        end
    end
    return new_dict
end
function transform_dict_static(dict::Dict{Symbol,Dict{Symbol,Vector{EdgeAutomatonGandF}}})
    new_dict = Dict{Symbol,Dict{Symbol,SVector}}()
    for from_loc in keys(dict)
        new_dict[from_loc] = Dict{Symbol,SVector}()
        for to_loc in keys(dict[from_loc])
            d = dict[from_loc][to_loc]
            new_dict[from_loc][to_loc] = SVector{length(d)}(d) 
        end
    end
    return new_dict
end
new_map_edges = transform_dict_tuple(A_G_F.map_edges)
new_edges_from_current_loc = new_map_edges[:l1G]
vec_edges = new_map_edges[:l1G][2]
#concrete_vec_edges = Vector{typeof(vec_edges[1])}(vec_edges)
svector_map_edges = transform_dict_static(A_G_F.map_edges)
svector_edges_from_current_loc = svector_map_edges[:l1G]
static_vec_edges = SVector{length(vec_edges)}(vec_edges)
tuple_edges = Tuple(vec_edges)

println("Test of the current implementation")
@btime begin
    for target_loc in keys(edges_from_current_loc)
        for edge in edges_from_current_loc[target_loc]
            test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
        end
    end
end
println("Test with pairs")
@btime begin
    for pair_target_loc_edges in edges_from_current_loc
        for edge in pair_target_loc_edges.second
            test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
        end
    end
end
println("New map")
@btime begin
    for target_loc_edges in new_edges_from_current_loc
        for edge in new_edges_from_current_loc[2]
            test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
        end
    end
end
println("Read of a edge vector 1 (collection iteration)")
@btime begin
    for edge in vec_edges
        test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
    end
end
println("Read of a edge vector 2 (eachindex)")
@btime begin
    for i = eachindex(vec_edges)
        test_cc = check_constraints_AutomatonGandF(vec_edges[i], time, values, x, p)
    end
end
println("Read of a edge vector 3 (steprange)")
@btime begin
    for i = 1:length(vec_edges)
        test_cc = check_constraints_AutomatonGandF(vec_edges[i], time, values, x, p)
    end
end
println("Read of a edge vector 4 (steprange 2)")
nb_edges = length(vec_edges)
@btime begin
    for i = 1:nb_edges
        test_cc = check_constraints_AutomatonGandF(vec_edges[i], time, values, x, p)
    end
end
println("Read of an edge vector 5 (static vectors)")
@btime begin
    for i = eachindex(static_vec_edges)
        test_cc = check_constraints_AutomatonGandF(static_vec_edges[i], time, values, x, p)
    end
end
println("Read of an edge vector 6 (tuples)")
@show typeof(tuple_edges)
@btime begin
    for i = eachindex(tuple_edges)
        test_cc = check_constraints_AutomatonGandF(tuple_edges[i], time, values, x, p)
    end
end
#=
println("Read of an edge vector 7 (concrete type)")
@show typeof(vec_edges)
@show typeof(concrete_vec_edges)
@btime begin
    for edge in concrete_vec_edges
        #test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
    end
end
=#
println("Read of an edge vector 8 (concrete type v2)")
@show vec_edges
concrete_2_vec_edges = Vector{Union{Nothing,Vector{Symbol}}}([e.transitions for e in vec_edges])
@btime begin
    for i = eachindex(concrete_2_vec_edges)
        concrete_2_vec_edges[i]
        #test_cc = check_constraints_AutomatonGandF(edge, time, values, x, p)
    end
end

b_find = @benchmark _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, time, values, x, p, false)
@show mean(b_find).time, mean(b_find).memory

b_svector_find = @benchmark svector_find_edge_candidates!(edge_candidates, svector_edges_from_current_loc, Λ, time, values, x, p, false)
@show mean(b_svector_find).time, mean(b_svector_find).memory

Profile.clear_malloc_data()
b_myfind = @benchmark my_find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, time, values, x, p, false)
@show mean(b_myfind).time, mean(b_myfind).memory

#=
The problem is that in each step of these benchmarks, an abstract type is involved
=#

