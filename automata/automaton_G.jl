
function create_automaton_G(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    @assert sym_obs in m.g
    # Locations
    locations = [:l0, :l1, :l2, :l3, :l4]

    # Invariant predicates
    true_inv_predicate = (x::Vector{Int}) -> return true 
    Λ_F = Dict(:l0 => true_inv_predicate, :l1 => true_inv_predicate,
               :l2 => true_inv_predicate, :l3 => true_inv_predicate, 
               :l4 => true_inv_predicate)
    
    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:l2]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:tprime => 1, :in => 2,
                                                        :n => 3,  :d => 4, :isabs => 5)

    ## Flow of variables
    flow = Dict{Location,Vector{Float64}}(:l0 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l1 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l2 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l3 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l4 => [1.0,0.0,0.0,0.0,0.0])

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end
    istrue(val::Float64) = convert(Bool, val)

    # l0 loc
    tuple = (:l0, :l1)
    cc_aut_G_l0l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l0l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1; 
         S[:d] = 0; 
         S[:n] = get_value(S, x, sym_obs); 
         S[:in] = true; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([nothing], cc_aut_G_l0l1_1, us_aut_G_l0l1_1!)
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    tuple = (:l1, :l3)
    cc_aut_G_l1l3_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        getfield(S, :time) <= constants[:t1] && 
        S[:n] < constants[:x1] || S[:n] > constants[:x2]
    us_aut_G_l1l3_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3; 
         S[:d] = min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:in] = false)
    edge1 = Edge([nothing], cc_aut_G_l1l3_1, us_aut_G_l1l3_1!)

    cc_aut_G_l1l3_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
         !istrue(S[:in]) && 
         (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
         (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1l3_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3; 
         S[:d] = S[:d] * (getfield(S, :time) - constants[:t1]); 
         S[:tprime] = 0.0)
    edge3 = Edge([nothing], cc_aut_G_l1l3_3, us_aut_G_l1l3_3!)
   
    cc_aut_G_l1l3_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) <= constants[:t1]) && 
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1l3_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3; 
         S[:d] = 0; 
         S[:in] = false)
    edge2 = Edge([nothing], cc_aut_G_l1l3_2, us_aut_G_l1l3_2!)

    cc_aut_G_l1l3_4(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1l3_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3; 
         S[:tprime] = 0.0)
    edge4 = Edge([nothing], cc_aut_G_l1l3_4, us_aut_G_l1l3_4!)
    
    map_edges[:l1][:l3] = [edge1, edge2, edge3, edge4]

    tuple = (:l1, :l4)
    cc_aut_G_l1l4_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2])
    us_aut_G_l1l4_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l4; 
         S[:d] += S[:d] * (getfield(S, :time) - constants[:t1]))
    edge1 = Edge([nothing], cc_aut_G_l1l4_1, us_aut_G_l1l4_1!)
    cc_aut_G_l1l4_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2])
    us_aut_G_l1l4_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l4)
    edge2 = Edge([nothing], cc_aut_G_l1l4_2, us_aut_G_l1l4_2!)
    
    map_edges[:l1][:l4] = [edge1, edge2]

    tuple = (:l1, :l2)
    cc_aut_G_l1l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l1l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2)
    edge1 = Edge([nothing], cc_aut_G_l1l2_1, us_aut_G_l1l2_1!)
    cc_aut_G_l1l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l1l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2; 
         S[:d] = S[:d] * (constants[:t2] - constants[:t1]))
    edge2 = Edge([nothing], cc_aut_G_l1l2_2, us_aut_G_l1l2_2!)
    cc_aut_G_l1l2_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:isabs]) && 
        getfield(S, :time) <= constants[:t1]
    us_aut_G_l1l2_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2; 
         S[:d] = (constants[:t2] - constants[:t1]) *
                   min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))
    edge3 = Edge([nothing], cc_aut_G_l1l2_3, us_aut_G_l1l2_3!)
    cc_aut_G_l1l2_4(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:isabs]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2])
    us_aut_G_l1l2_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2; 
         S[:d] += (constants[:t2] - getfield(S, :time)) * 
                    min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))
    edge4 = Edge([nothing], cc_aut_G_l1l2_4, us_aut_G_l1l2_4!)
    
    map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]

    # l3 loc
    tuple = (:l3, :l1)
    cc_aut_G_l3l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l3l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1; 
         S[:n] = get_value(S, x, sym_obs); 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_G_l3l1_1, us_aut_G_l3l1_1!)
    
    map_edges[:l3][:l1] = [edge1]

    tuple = (:l3, :l2)
    cc_aut_G_l3l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l3l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2;
         S[:d] = S[:d] * (constants[:t2] - constants[:t1]))
    edge2 = Edge([nothing], cc_aut_G_l3l2_2, us_aut_G_l3l2_2!)
    cc_aut_G_l3l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l3l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2)
    edge1 = Edge([nothing], cc_aut_G_l3l2_1, us_aut_G_l3l2_1!)
 
    map_edges[:l3][:l2] = [edge1, edge2]

    # l4 loc
    tuple = (:l4, :l1)
    cc_aut_G_l4l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l4l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1; 
         S[:d] += S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:tprime] = 0.0; 
         S[:n] = get_value(S, x, sym_obs); 
         S[:in] = true; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_G_l4l1_1, us_aut_G_l4l1_1!)
    
    map_edges[:l4][:l1] = [edge1]

    tuple = (:l4, :l2)
    cc_aut_G_l4l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t2])
    us_aut_G_l4l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2; 
         S[:d] +=  S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:tprime] = 0.0)
    edge1 = Edge([nothing], cc_aut_G_l4l2_1, us_aut_G_l4l2_1!)
    
    map_edges[:l4][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA(m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
   
end

export create_automaton_G

