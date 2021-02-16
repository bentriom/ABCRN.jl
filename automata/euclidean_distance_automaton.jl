
function create_euclidean_distance_automaton(m::ContinuousTimeModel, timeline::AbstractVector{Float64}, observations::AbstractVector{Float64}, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert length(timeline) == length(observations) "Timeline and observations vectors don't have the same length"
    nbr_observations = length(observations)

    # Locations
    locations = [:l0, :l1, :l2]

    ## Invariant predicates
    true_inv_predicate(x::Vector{Int}) = true 
    Λ_F = Dict(:l0 => true_inv_predicate, :l1 => true_inv_predicate,
               :l2 => true_inv_predicate)

    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:l2]

    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:t => 1, :n => 2, 
                                                        :d => 3, :idx => 4)
    vector_flow = [1.0, 0.0, 0.0, 0.0]
    flow = Dict{Location,Vector{Float64}}(:l0 => vector_flow, 
                                          :l1 => vector_flow, 
                                          :l2 => vector_flow)

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end

    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    idx_var_t = map_var_automaton_idx[:t] 
    idx_var_n = map_var_automaton_idx[:n] 
    idx_var_d = map_var_automaton_idx[:d] 
    idx_var_idx = map_var_automaton_idx[:idx] 
    
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_eucl_dist_aut_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote
        # l0 loc
        # l0 => l1
        @everywhere $(func_name(:cc, :l0, :l1, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l0, :l1, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        (setfield!(S, :loc, Symbol("l1")); 
         setindex!(getfield(S, :values), x[$(idx_obs_var)], $(idx_var_n));
         setindex!(getfield(S, :values), 0.0, $(idx_var_d));
         setindex!(getfield(S, :values), 1.0, $(idx_var_idx)))

        # l1 loc
        # l1 => l1
        # Defined below 
        @everywhere $(func_name(:cc, :l1, :l1, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (tml = $(Tuple(timeline));
         tml_idx = tml[convert(Int, getfield(S, :values)[$(idx_var_idx)])];
         getfield(S, :values)[$(idx_var_t)] >= tml_idx)
        @everywhere $(func_name(:us, :l1, :l1, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (y_obs = $(Tuple(observations));
         y_obs_idx = y_obs[convert(Int, getfield(S, :values)[$(idx_var_idx)])];
         setindex!(getfield(S, :values), getfield(S, :values)[$(idx_var_d)] + (getfield(S, :values)[$(idx_var_n)]  - y_obs_idx)^2, 
                                         $(idx_var_d));
         setindex!(getfield(S, :values), getfield(S, :values)[$(idx_var_idx)] + 1.0, $(idx_var_idx)))
        
        @everywhere $(func_name(:cc, :l1, :l1, 2))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = true 
        @everywhere $(func_name(:us, :l1, :l1, 2))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        (setindex!(getfield(S, :values), x[$(idx_obs_var)], $(idx_var_n)))
        
        # l1 => l2
        @everywhere $(func_name(:cc, :l1, :l2, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        getfield(S, :values)[$(idx_var_idx)] >= ($nbr_observations + 1)
        @everywhere $(func_name(:us, :l1, :l2, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        (setfield!(S, :loc, Symbol("l2")); 
         setindex!(getfield(S, :values), sqrt(getfield(S, :values)[$(idx_var_d)]), $(idx_var_d)))
    end
    eval(meta_elementary_functions)

    # l0 loc
    # l0 => l1
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :l1, 1)), getfield(Main, func_name(:us, :l0, :l1, 1)))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l1
    #=
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l1, :l1, 1)), getfield(Main, func_name(:us, :l1, :l2, 1)))
    map_edges[:l1][:l1] = [edge1]
    for i = 1:nbr_observations
        meta_edge_i = quote
            @everywhere $(func_name(:cc, :l1, :l1, 1+i))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
            S[:t] >= $(timeline[i])
            @everywhere $(func_name(:us, :l1, :l1, 1+i))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
            (setindex!(getfield(S, :values), S[:d] + (S[:n]  - $(observations[i]))^2, $(idx_var_d));
             setindex!(getfield(S, :values), S[:idx] + 1.0, $(idx_var_idx)))
        end
        eval(meta_edge_i)
        push!(map_edges[:l1][:l1], Edge(nothing, getfield(Main, func_name(:cc, :l1, :l1, 1+i)), getfield(Main, func_name(:us, :l1, :l1, 1+i))))
    end
    =#
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l1, 1)), getfield(Main, func_name(:us, :l1, :l1, 1)))
    edge2 = Edge([:ALL], getfield(Main, func_name(:cc, :l1, :l1, 2)), getfield(Main, func_name(:us, :l1, :l1, 2)))
    map_edges[:l1][:l1] = [edge1, edge2]

    # l1 => l2
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 1)), getfield(Main, func_name(:us, :l1, :l2, 1)))
    map_edges[:l1][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)
    for i = 1:nbr_observations
        constants[Symbol("tml_$(convert(Float64, i))")] = timeline[i]
        constants[Symbol("y_$(convert(Float64, i))")] = observations[i]
    end

    A = LHA("Euclidean distance", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_euclidean_distance_automaton

