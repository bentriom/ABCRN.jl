
function create_automaton_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    # Locations
    locations = [:l0, :l1, :l2, :l3]

    ## Invariant predicates
    true_inv_predicate = (x::Vector{Int}) -> return true 
    Λ_F = Dict(:l0 => true_inv_predicate, :l1 => true_inv_predicate,
               :l2 => true_inv_predicate, :l3 => true_inv_predicate)
    
    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:l2]

    #S.n <=> S.values[A.map_var_automaton_idx[:n]] 
    #P <=> xn[map_var_model_idx[constants[str_O]] with str_O = :P. On stock str_O dans constants
    # P = get_value(S, x, sym_obs) 
    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:n => 1, :d => 2, :isabs => 3)

    ## Flow of variables
    flow = Dict{Location,Vector{Float64}}(:l0 => [0.0,0.0,0.0], 
                                          :l1 => [0.0,0.0,0.0], 
                                          :l2 => [0.0,0.0,0.0], 
                                          :l3 => [0.0,0.0,0.0])

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end
    istrue(val::Float64) = convert(Bool, val)
    
    # l0 loc : we construct  the edges of the form l0 => (..)
    # "cc" as check_constraints
    tuple = (:l0, :l1)
    cc_aut_F_l0l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = true
    us_aut_F_l0l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1; 
         S[:n] = get_value(S, x, sym_obs);
         S[:d] = Inf; 
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([nothing], cc_aut_F_l0l1_1, us_aut_F_l0l1_1!)
    map_edges[:l0][:l1] = [edge1]

    # l1 loc : we construct  the edges of the form l1 => (..)
    tuple = (:l1, :l2)
    cc_aut_F_l1l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        getfield(S, :time) >= constants[:t1] &&
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_F_l1l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2;
         S[:d] = 0)
    edge1 = Edge([nothing], cc_aut_F_l1l2_1, us_aut_F_l1l2_1!)
    
    cc_aut_F_l1l2_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        getfield(S, :time) >= constants[:t1] &&
        S[:d] == 0 
    us_aut_F_l1l2_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2)
    edge4 = Edge([nothing], cc_aut_F_l1l2_4, us_aut_F_l1l2_4!)

    cc_aut_F_l1l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t2]) && 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2])
    us_aut_F_l1l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2;
         S[:d] = min(abs(S[:n] - constants[:x1]), abs(S[:n] - constants[:x2])))
    edge2 = Edge([nothing], cc_aut_F_l1l2_2, us_aut_F_l1l2_2!)

    cc_aut_F_l1l2_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        istrue(S[:isabs]) && getfield(S, :time) <= constants[:t2]
    us_aut_F_l1l2_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2)
    edge3 = Edge([nothing], cc_aut_F_l1l2_3, us_aut_F_l1l2_3!)

    map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]

    tuple = (:l1, :l3)
    cc_aut_F_l1l3_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (constants[:x1] <= S[:n] <= constants[:x2])
    us_aut_F_l1l3_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3;
         S[:d] = 0;)
    edge1 = Edge([nothing], cc_aut_F_l1l3_1, us_aut_F_l1l3_1!)
    
    cc_aut_F_l1l3_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2]) && 
        (getfield(S, :time) <= constants[:t1])
    us_aut_F_l1l3_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3;
         S[:d] = min(sqrt((getfield(S, :time) - constants[:t1])^2 + (S[:n] - constants[:x2])^2), 
                      sqrt((getfield(S, :time) - constants[:t1])^2 + (S[:n] - constants[:x1])^2)))
    edge2 = Edge([nothing], cc_aut_F_l1l3_2, us_aut_F_l1l3_2!)

    cc_aut_F_l1l3_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S[:n] < constants[:x1] || S[:n] > constants[:x2]) && 
        (constants[:t1] <= getfield(S, :time) <= constants[:t2])
    us_aut_F_l1l3_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l3;
         S[:d] = min(S[:d], min(abs(S[:n] - constants[:x1]), abs(S[:n] - constants[:x2]))))
    edge3 = Edge([nothing], cc_aut_F_l1l3_3, us_aut_F_l1l3_3!)
    map_edges[:l1][:l3] = [edge1, edge2, edge3]

    # l3 loc
    tuple = (:l3, :l1)
    cc_aut_F_l3l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = true
    us_aut_F_l3l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l1;
         S[:n] = get_value(S, x, sym_obs);
         S[:isabs] = getfield(m, :isabsorbing)(getfield(m, :p),x))
    edge1 = Edge([:ALL], cc_aut_F_l3l1_1, us_aut_F_l3l1_1!)
    map_edges[:l3][:l1] = [edge1]
    
    
    tuple = (:l3, :l2)
    cc_aut_F_l3l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (getfield(S, :time) >= constants[:t2])
    us_aut_F_l3l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}) = 
        (S.loc = :l2)
    edge1 = Edge([nothing], cc_aut_F_l3l2_1, us_aut_F_l3l2_1!)
    map_edges[:l3][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA(m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_F

