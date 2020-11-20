
function create_automaton_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, str_obs::String)
    # Locations
    l_loc = ["l0", "l1", "l2", "l3"]

    ## Invariant predicates
    true_inv_predicate = (A::LHA, S:: StateLHA) -> return true 
    Λ_F = Dict("l0" => true_inv_predicate, "l1" => true_inv_predicate,
               "l2" => true_inv_predicate, "l3" => true_inv_predicate)
    
    ## Init and final loc
    l_loc_init = ["l0"]
    l_loc_final = ["l2"]

    #S.n <=> S.l_var[A.map_var_automaton_idx["n"]] 
    #P <=> xn[map_var_model_idx[l_ctes[str_O]] with str_O = "P". On stock str_O dans l_ctes
    
    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}("n" => 1, "d" => 2)

    ## Flow of variables
    l_flow = Dict{VariableAutomaton,Vector{Float64}}("l0" => [0.0,0.0], 
                                                     "l1" => [0.0,0.0], 
                                                     "l2" => [0.0,0.0], 
                                                     "l3" => [0.0,0.0])

    ## Edges
    map_edges = Dict{Tuple{Location,Location}, Vector{Edge}}()

    # l0 loc : we construct  the edges of the form l0 => (..)
    # "cc" as check_constraints
    tuple = ("l0", "l1")
    cc_aut_F_l0l1_1(A::LHA, S::StateLHA) = true
    us_aut_F_l0l1_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S.loc = "l1"; S["d"] = Inf; S["n"] = get_value(A, x, str_obs))
    edge1 = Edge([nothing], cc_aut_F_l0l1_1, us_aut_F_l0l1_1!)
    map_edges[tuple] = [edge1]

    # l1 loc : we construct  the edges of the form l1 => (..)
    tuple = ("l1", "l2")
    cc_aut_F_l1l2_1(A::LHA, S::StateLHA) = 
        (A.l_ctes["x1"] <= S["n"] <= A.l_ctes["x2"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"])
    us_aut_F_l1l2_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S["d"] = 0; S.loc = "l2")
    edge1 = Edge([nothing], cc_aut_F_l1l2_1, us_aut_F_l1l2_1!)

    cc_aut_F_l1l2_2(A::LHA, S::StateLHA) = (S["d"] > 0 && S.time > A.l_ctes["t2"])
    us_aut_F_l1l2_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2")
    edge2 = Edge([nothing], cc_aut_F_l1l2_2, us_aut_F_l1l2_2!)
    
    cc_aut_F_l1l2_3(A::LHA, S::StateLHA) = (S["d"] == 0 && S.time >= A.l_ctes["t1"])
    us_aut_F_l1l2_3!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2")
    edge3 = Edge([nothing], cc_aut_F_l1l2_3, us_aut_F_l1l2_3!)
    
    map_edges[tuple] = [edge1, edge2, edge3]

    tuple = ("l1", "l3")
    cc_aut_F_l1l3_1(A::LHA, S::StateLHA) = (A.l_ctes["x1"] <= S["n"] <= A.l_ctes["x2"])
    us_aut_F_l1l3_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S["d"] = 0; S.loc = "l3")
    edge1 = Edge([nothing], cc_aut_F_l1l3_1, us_aut_F_l1l3_1!)
    
    cc_aut_F_l1l3_2(A::LHA, S::StateLHA) = 
        (A.l_ctes["x1"] > S["n"] || S["n"] > A.l_ctes["x2"]) && (S.time < A.l_ctes["t1"])
    us_aut_F_l1l3_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S["d"] = min(sqrt((S.time - A.l_ctes["t1"])^2 + (S["n"] - A.l_ctes["x2"])^2), 
                      sqrt((S.time - A.l_ctes["t1"])^2 + (S["n"] - A.l_ctes["x1"])^2)); S.loc = "l3")
    edge2 = Edge([nothing], cc_aut_F_l1l3_2, us_aut_F_l1l3_2!)

    cc_aut_F_l1l3_3(A::LHA, S::StateLHA) = 
        (A.l_ctes["x1"] > S["n"] || S["n"] > A.l_ctes["x2"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"])
    us_aut_F_l1l3_3!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S["d"] = min(S["d"], min(abs(S["n"] - A.l_ctes["x1"]), abs(S["n"] - A.l_ctes["x2"]))); S.loc = "l3")
    edge3 = Edge([nothing], cc_aut_F_l1l3_3, us_aut_F_l1l3_3!)
    map_edges[tuple] = [edge1, edge2, edge3]

    # l3 loc
    tuple = ("l3", "l1")
    cc_aut_F_l3l1_1(A::LHA, S::StateLHA) = true
    us_aut_F_l3l1_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S["n"] = get_value(A, x, str_obs); S.loc = "l1")
    edge1 = Edge(["ALL"], cc_aut_F_l3l1_1, us_aut_F_l3l1_1!)
    tuple = ("l3", "l2")
    cc_aut_F_l3l2_1(A::LHA, S::StateLHA) = (S.time >= A.l_ctes["t2"])
    us_aut_F_l3l2_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S["n"] = get_value(A, x, str_obs); S.loc = "l2")
    edge2 = Edge([nothing], cc_aut_F_l3l2_1, us_aut_F_l3l2_1!)
    map_edges[tuple] = [edge1, edge2]

    ## Constants
    l_ctes = Dict{String,Float64}("x1" => x1, "x2" => x2, "t1" => t1, "t2" => t2)

    A = LHA(m.l_transitions, l_loc, Λ_F, l_loc_init, l_loc_final, 
            map_var_automaton_idx, l_flow, map_edges, l_ctes, m.map_var_idx)
    return A
end

export create_automaton_F

