
function create_automaton_G(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, str_obs::String)
    @assert str_obs in m.g
    # Locations
    l_loc = ["l0", "l1", "l2", "l3", "l4"]

    # Invariant predicates
    true_inv_predicate = (A::LHA, S:: StateLHA) -> return true 
    Λ_F = Dict("l0" => true_inv_predicate, "l1" => true_inv_predicate,
               "l2" => true_inv_predicate, "l3" => true_inv_predicate, 
               "l4" => true_inv_predicate)
    
    ## Init and final loc
    l_loc_init = ["l0"]
    l_loc_final = ["l2"]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}("t'" => 1, "in" => 2,
                                                         "n" => 3,  "d" => 4)

    ## Flow of variables
    l_flow = Dict{VariableAutomaton,Vector{Float64}}("l0" => [0.0,0.0,0.0,0.0], 
                                                     "l1" => [0.0,0.0,0.0,0.0], 
                                                     "l2" => [0.0,0.0,0.0,0.0], 
                                                     "l3" => [0.0,0.0,0.0,0.0], 
                                                     "l4" => [1.0,0.0,0.0,0.0])

    ## Edges
    map_edges = Dict{Tuple{Location,Location}, Vector{Edge}}()

    isin(val::Float64) = convert(Bool, val)

    # l0 loc
    tuple = ("l0", "l1")
    cc_aut_G_l0l1_1(A::LHA, S::StateLHA) = true
    us_aut_G_l0l1_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l1"; S["d"] = 0; S["n"] = get_value(A, x, str_obs); S["in"] = true)
    edge1 = Edge([nothing], cc_aut_G_l0l1_1, us_aut_G_l0l1_1!)
    map_edges[tuple] = [edge1]

    # l1 loc
    tuple = ("l1", "l3")
    cc_aut_G_l1l3_1(A::LHA, S::StateLHA) = 
        (S.time < A.l_ctes["t1"] && (S["n"] < A.l_ctes["x1"] || S["n"] > A.l_ctes["x2"]))
    us_aut_G_l1l3_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S.loc = "l3"; S["d"] = min(abs(A.l_ctes["x1"] - S["n"]), abs(A.l_ctes["x2"] - S["n"])); S["in"] = false)
    edge1 = Edge([nothing], cc_aut_G_l1l3_1, us_aut_G_l1l3_1!)
    cc_aut_G_l1l3_2(A::LHA, S::StateLHA) = 
        (S.time < A.l_ctes["t1"] && (A.l_ctes["x1"] <= S["n"] <= A.l_ctes["x2"]))
    us_aut_G_l1l3_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S.loc = "l3"; S["d"] = 0; S["in"] = false)
    edge2 = Edge([nothing], cc_aut_G_l1l3_2, us_aut_G_l1l3_2!)
    cc_aut_G_l1l3_3(A::LHA, S::StateLHA) = 
        (!isin(S["in"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"]) && (A.l_ctes["x1"] <= S["n"] <= A.l_ctes["x2"]))
    us_aut_G_l1l3_3!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l3"; S["d"] = S["d"] * (S.time - A.l_ctes["t1"]); S["t'"] = 0.0)
    edge3 = Edge([nothing], cc_aut_G_l1l3_3, us_aut_G_l1l3_3!)
    cc_aut_G_l1l3_4(A::LHA, S::StateLHA) = 
        (isin(S["in"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"]) && (A.l_ctes["x1"] <= S["n"] <= A.l_ctes["x2"]))
    us_aut_G_l1l3_4!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l3"; S["t'"] = 0.0)
    edge4 = Edge([nothing], cc_aut_G_l1l3_4, us_aut_G_l1l3_4!)
    map_edges[tuple] = [edge1, edge2, edge3, edge4]

    tuple = ("l1", "l4")
    cc_aut_G_l1l4_1(A::LHA, S::StateLHA) = 
        (!isin(S["in"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"]) && (S["n"] < A.l_ctes["x1"] || S["n"] > A.l_ctes["x2"]))
    us_aut_G_l1l4_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l4"; S["d"] += S["d"] * (S.time - A.l_ctes["t1"]))
    edge1 = Edge([nothing], cc_aut_G_l1l4_1, us_aut_G_l1l4_1!)
    cc_aut_G_l1l4_2(A::LHA, S::StateLHA) = 
        (isin(S["in"]) && (A.l_ctes["t1"] <= S.time <= A.l_ctes["t2"]) && (S["n"] < A.l_ctes["x1"] || S["n"] > A.l_ctes["x2"]))
    us_aut_G_l1l4_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l4")
    edge2 = Edge([nothing], cc_aut_G_l1l4_2, us_aut_G_l1l4_2!)
    map_edges[tuple] = [edge1, edge2]

    tuple = ("l1", "l2")
    cc_aut_G_l1l2_1(A::LHA, S::StateLHA) = (isin(S["in"]) && S.time >= A.l_ctes["t2"])
    us_aut_G_l1l2_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2")
    edge1 = Edge([nothing], cc_aut_G_l1l2_1, us_aut_G_l1l2_1!)
    cc_aut_G_l1l2_2(A::LHA, S::StateLHA) = (!isin(S["in"]) && S.time >= A.l_ctes["t2"])
    us_aut_G_l1l2_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2"; S["d"] = S["d"] * (A.l_ctes["t2"] - A.l_ctes["t1"]))
    edge2 = Edge([nothing], cc_aut_G_l1l2_2, us_aut_G_l1l2_2!)
    map_edges[tuple] = [edge1, edge2]

    # l3 loc
    tuple = ("l3", "l1")
    cc_aut_G_l3l1_1(A::LHA, S::StateLHA) = true
    us_aut_G_l3l1_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l1"; S["n"] = get_value(A, x, str_obs))
    edge1 = Edge(["ALL"], cc_aut_G_l3l1_1, us_aut_G_l3l1_1!)
    map_edges[tuple] = [edge1]

    tuple = ("l3", "l2")
    cc_aut_G_l3l2_1(A::LHA, S::StateLHA) = (isin(S["in"]) && S.time >= A.l_ctes["t2"])
    us_aut_G_l3l2_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2")
    edge1 = Edge([nothing], cc_aut_G_l3l2_1, us_aut_G_l3l2_1!)
    cc_aut_G_l3l2_2(A::LHA, S::StateLHA) = (!isin(S["in"]) && S.time >= A.l_ctes["t2"])
    us_aut_G_l3l2_2!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = (S.loc = "l2"; S["d"] = S["d"] * (A.l_ctes["t2"] - A.l_ctes["t1"]))
    edge2 = Edge([nothing], cc_aut_G_l3l2_2, us_aut_G_l3l2_2!)
    map_edges[tuple] = [edge1, edge2]

    # l4 loc
    tuple = ("l4", "l1")
    cc_aut_G_l4l1_1(A::LHA, S::StateLHA) = true
    us_aut_G_l4l1_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S.loc = "l1"; S["d"] += S["t'"] * min(abs(A.l_ctes["x1"] - S["n"]), abs(A.l_ctes["x2"] - S["n"])); 
         S["t'"] = 0.0; S["n"] = get_value(A, x, str_obs); S["in"] = true)
    edge1 = Edge(["ALL"], cc_aut_G_l4l1_1, us_aut_G_l4l1_1!)
    map_edges[tuple] = [edge1]

    tuple = ("l4", "l2")
    cc_aut_G_l4l2_1(A::LHA, S::StateLHA) = (S.time >= A.l_ctes["t2"])
    us_aut_G_l4l2_1!(A::LHA, S::StateLHA, x::SubArray{Int,1}) = 
        (S.loc = "l2"; S["d"] +=  S["t'"] * min(abs(A.l_ctes["x1"] - S["n"]), abs(A.l_ctes["x2"] - S["n"])); S["t'"] = 0.0)
    edge1 = Edge([nothing], cc_aut_G_l4l2_1, us_aut_G_l4l2_1!)
    map_edges[tuple] = [edge1]

    ## Constants
    l_ctes = Dict{String,Float64}("x1" => x1, "x2" => x2, "t1" => t1, "t2" => t2)

    A = LHA(m.l_transitions, l_loc, Λ_F, l_loc_init, l_loc_final, 
            map_var_automaton_idx, l_flow, map_edges, l_ctes, m.map_var_idx)
    return A
end

export create_automaton_G

