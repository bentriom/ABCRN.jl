
function create_automaton_G_and_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs_G::VariableModel,
                                  x3::Float64, x4::Float64, t3::Float64, t4::Float64, sym_obs_F::VariableModel)

    @assert sym_obs_G in m.g
    @assert sym_obs_F in m.g
    # Locations
    locations = [:l0G, :l1G, :l2G, :l3G, :l4G,
                 :l1F, :l2F, :l3F]

    # Invariant predicates
    true_inv_predicate = (x::Vector{Int}) -> return true 
    Λ_F = Dict(:l0G => true_inv_predicate, :l1G => true_inv_predicate,
               :l2G => true_inv_predicate, :l3G => true_inv_predicate, 
               :l4G => true_inv_predicate,
               :l1F => true_inv_predicate,
               :l2F => true_inv_predicate, :l3F => true_inv_predicate)
    
    ## Init and final loc
    locations_init = [:l0G]
    locations_final = [:l2F]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:tprime => 1, :in => 2,
                                                        :n => 3,  :d => 4, 
                                                        :dprime => 5, :isabs => 6)

    ## Flow of variables
    flow = Dict{Location,Vector{Float64}}(:l0G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l1G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l2G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l3G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l4G => [1.0,0.0,0.0,0.0,0.0,0.0],
                                          :l1F => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l2F => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l3F => [0.0,0.0,0.0,0.0,0.0,0.0])

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end

    istrue(val::Float64) = convert(Bool, val)

    # l0G loc
    tuple = (:l0G, :l1G)
    cc_aut_G_l0Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l0Gl1G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1G; 
         S[:d] = 0; 
         S[:n] = get_value(S, x, sym_obs_G); 
         S[:in] = true; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([nothing], cc_aut_G_l0Gl1G_1, us_aut_G_l0Gl1G_1!)
    map_edges[:l0G][:l1G] = [edge1]

    # l1G loc
    tuple = (:l1G, :l3G)
    cc_aut_G_l1Gl3G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        getfield(S, :time) <= constants[:t1] && 
        S[:n] < constants[:x1] || S[:n] > constants[:x2]
    us_aut_G_l1Gl3G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3G; 
         S[:d] = min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:in] = false)
    edge1 = Edge([nothing], cc_aut_G_l1Gl3G_1, us_aut_G_l1Gl3G_1!)

    cc_aut_G_l1Gl3G_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
         !istrue(S[:in]) && 
         (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
         (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1Gl3G_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3G; 
         S[:d] = S[:d] * (getfield(S, :time) - constants[:t1]); 
         S[:tprime] = 0.0)
    edge3 = Edge([nothing], cc_aut_G_l1Gl3G_3, us_aut_G_l1Gl3G_3!)
   
    cc_aut_G_l1Gl3G_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) <= constants[:t1]) && 
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1Gl3G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3G; 
         S[:d] = 0; 
         S[:in] = false)
    edge2 = Edge([nothing], cc_aut_G_l1Gl3G_2, us_aut_G_l1Gl3G_2!)

    cc_aut_G_l1Gl3G_4(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_G_l1Gl3G_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3G; 
         S[:tprime] = 0.0)
    edge4 = Edge([nothing], cc_aut_G_l1Gl3G_4, us_aut_G_l1Gl3G_4!)
    
    map_edges[:l1G][:l3G] = [edge1, edge2, edge3, edge4]

    tuple = (:l1G, :l4G)
    cc_aut_G_l1Gl4G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2])
    us_aut_G_l1Gl4G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l4G; 
         S[:d] += S[:d] * (getfield(S, :time) - constants[:t1]))
    edge1 = Edge([nothing], cc_aut_G_l1Gl4G_1, us_aut_G_l1Gl4G_1!)
    
    cc_aut_G_l1Gl4G_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2])
    us_aut_G_l1Gl4G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l4G)
    edge2 = Edge([nothing], cc_aut_G_l1Gl4G_2, us_aut_G_l1Gl4G_2!)
    
    map_edges[:l1G][:l4G] = [edge1, edge2]

    tuple = (:l1G, :l2G)
    cc_aut_G_l1Gl2G_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:isabs]) && 
        getfield(S, :time) <= constants[:t1]
    us_aut_G_l1Gl2G_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G; 
        S[:d] = (constants[:t2] - constants[:t1]) *
        min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))
    edge3 = Edge([nothing], cc_aut_G_l1Gl2G_3, us_aut_G_l1Gl2G_3!)
        cc_aut_G_l1Gl2G_4(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:isabs]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2])
    us_aut_G_l1Gl2G_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G; 
        S[:d] += (constants[:t2] - getfield(S, :time)) * 
        min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))
    edge4 = Edge([nothing], cc_aut_G_l1Gl2G_4, us_aut_G_l1Gl2G_4!)
    cc_aut_G_l1Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l1Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G)
    edge1 = Edge([nothing], cc_aut_G_l1Gl2G_1, us_aut_G_l1Gl2G_1!)
    cc_aut_G_l1Gl2G_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l1Gl2G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G; 
         S[:d] = S[:d] * (constants[:t2] - constants[:t1]))
    edge2 = Edge([nothing], cc_aut_G_l1Gl2G_2, us_aut_G_l1Gl2G_2!)
   
    map_edges[:l1G][:l2G] = [edge3, edge4, edge1, edge2]

    # l3G loc
    tuple = (:l3G, :l1G)
    cc_aut_G_l3Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l3Gl1G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1G; 
         S[:n] = get_value(S, x, sym_obs_G); 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_G_l3Gl1G_1, us_aut_G_l3Gl1G_1!)
    map_edges[:l3G][:l1G] = [edge1]

    tuple = (:l3G, :l2G)
    cc_aut_G_l3Gl2G_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l3Gl2G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G;
         S[:d] = S[:d] * (constants[:t2] - constants[:t1]))
    edge2 = Edge([nothing], cc_aut_G_l3Gl2G_2, us_aut_G_l3Gl2G_2!)
    
    cc_aut_G_l3Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        !istrue(S[:in]) && 
        getfield(S, :time) >= constants[:t2]
    us_aut_G_l3Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G)
    edge1 = Edge([nothing], cc_aut_G_l3Gl2G_1, us_aut_G_l3Gl2G_1!)
 
    map_edges[:l3G][:l2G] = [edge1, edge2]

    # l4G loc
    tuple = (:l4G, :l1G)
    cc_aut_G_l4Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_G_l4Gl1G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1G; 
         S[:d] += S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:tprime] = 0.0; 
         S[:n] = get_value(S, x, sym_obs_G); 
         S[:in] = true; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_G_l4Gl1G_1, us_aut_G_l4Gl1G_1!)
    
    map_edges[:l4G][:l1G] = [edge1]

    tuple = (:l4G, :l2G)
    cc_aut_G_l4Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t2])
    us_aut_G_l4Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2G; 
         S[:d] +=  S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
         S[:tprime] = 0.0)
    edge1 = Edge([nothing], cc_aut_G_l4Gl2G_1, us_aut_G_l4Gl2G_1!)

    map_edges[:l4G][:l2G] = [edge1]
    
    # Connection between the two automata: l2G => l1F
    tuple = (:l2G, :l1F)
    cc_aut_F_l2Gl1F_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_F_l2Gl1F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1F; 
         S[:n] = get_value(S, x, sym_obs_F);
         S[:dprime] = Inf; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([nothing], cc_aut_F_l2Gl1F_1, us_aut_F_l2Gl1F_1!)
    map_edges[:l2G][:l1F] = [edge1]

    # l1F loc : we construct  the edges of the form l1F => (..)
    tuple = (:l1F, :l2F)
    cc_aut_F_l1Fl2F_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        getfield(S, :time) >= constants[:t3] &&
        (constants[:x3] <= S[:n] <= constants[:x4])
    us_aut_F_l1Fl2F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2F;
         S[:dprime] = 0)
    edge1 = Edge([nothing], cc_aut_F_l1Fl2F_1, us_aut_F_l1Fl2F_1!)
    
    cc_aut_F_l1Fl2F_4(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        getfield(S, :time) >= constants[:t3] &&
        S[:dprime] == 0 
    us_aut_F_l1Fl2F_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2F)
    edge4 = Edge([nothing], cc_aut_F_l1Fl2F_4, us_aut_F_l1Fl2F_4!)

    cc_aut_F_l1Fl2F_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t4]) && 
        (S[:n] < constants[:x3] || S[:n] > constants[:x4])
    us_aut_F_l1Fl2F_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2F;
         S[:dprime] = min(abs(S[:n] - constants[:x3]), abs(S[:n] - constants[:x4]));
         S[:d] += S[:dprime])
    edge2 = Edge([nothing], cc_aut_F_l1Fl2F_2, us_aut_F_l1Fl2F_2!)

    cc_aut_F_l1Fl2F_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        istrue(S[:isabs]) && getfield(S, :time) <= constants[:t4]
    us_aut_F_l1Fl2F_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2F;
         S[:d] += S[:dprime])
    edge3 = Edge([nothing], cc_aut_F_l1Fl2F_3, us_aut_F_l1Fl2F_3!)

    map_edges[:l1F][:l2F] = [edge1, edge4, edge3, edge2]

    tuple = (:l1F, :l3F)
    cc_aut_F_l1Fl3F_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (constants[:x3] <= S[:n] <= constants[:x4])
    us_aut_F_l1Fl3F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3F;
         S[:dprime] = 0;)
    edge1 = Edge([nothing], cc_aut_F_l1Fl3F_1, us_aut_F_l1Fl3F_1!)
    
    cc_aut_F_l1Fl3F_2(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (S[:n] < constants[:x3] || S[:n] > constants[:x4]) && 
        (getfield(S, :time) <= constants[:t3])
    us_aut_F_l1Fl3F_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3F;
         S[:dprime] = min(sqrt((getfield(S, :time) - constants[:t3])^2 + (S[:n] - constants[:x4])^2), 
                           sqrt((getfield(S, :time) - constants[:t3])^2 + (S[:n] - constants[:x3])^2)))
    edge2 = Edge([nothing], cc_aut_F_l1Fl3F_2, us_aut_F_l1Fl3F_2!)

    cc_aut_F_l1Fl3F_3(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (S[:n] < constants[:x3] || S[:n] > constants[:x4]) && 
        (constants[:t3] <= getfield(S, :time) <= constants[:t4])
    us_aut_F_l1Fl3F_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3F;
         S[:dprime] = min(S[:dprime], min(abs(S[:n] - constants[:x3]), abs(S[:n] - constants[:x4]))))
    edge3 = Edge([nothing], cc_aut_F_l1Fl3F_3, us_aut_F_l1Fl3F_3!)
    
    map_edges[:l1F][:l3F] = [edge1, edge2, edge3]

    # l3F loc
    tuple = (:l3F, :l1F)
    cc_aut_F_l3Fl1F_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = true
    us_aut_F_l3Fl1F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1F;
         S[:n] = get_value(S, x, sym_obs_F);
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_F_l3Fl1F_1, us_aut_F_l3Fl1F_1!)
    
    map_edges[:l3F][:l1F] = [edge1]
    
    tuple = (:l3F, :l2F)
    cc_aut_F_l3Fl2F_1(S::StateLHA, constants::Dict{Symbol,Float64}, xn::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t4])
    us_aut_F_l3Fl2F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2F)
    edge1 = Edge([nothing], cc_aut_F_l3Fl2F_1, us_aut_F_l3Fl2F_1!)
    
    map_edges[:l3F][:l2F] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2,
                                     :x3 => x3,  :x4 => x4, :t3 => t3, :t4 => t4)
    
    A = LHA(m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_G_and_F

