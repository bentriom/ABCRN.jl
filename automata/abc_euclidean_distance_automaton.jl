
@everywhere import FunctionWrappers: FunctionWrapper
@everywhere const ABCCheckConstraintsFunction = FunctionWrapper{Bool,Tuple{Float64,Vector{Float64},Vector{Int},Vector{Float64},Float64}}
@everywhere const ABCUpdateStateFunction = FunctionWrapper{Symbol,Tuple{Float64,Vector{Float64},Vector{Int},Vector{Float64},Float64}}

@everywhere struct EdgeABCEuclideanDistanceAutomaton <: Edge
    transitions::TransitionSet
    check_constraints::ABCCheckConstraintsFunction
    update_state!::ABCUpdateStateFunction
end

# Creation of the automaton types
@everywhere mutable struct ABCEuclideanDistanceAutomaton <: LHA
    transitions::Vector{Transition}
    locations::Vector{Location} 
    Λ::Dict{Location,InvariantPredicateFunction}
    locations_init::Vector{Location}
    locations_final::Vector{Location}
    map_var_automaton_idx::Dict{VariableAutomaton,Int} # nvar keys : str_var => idx in values
    flow::Dict{Location,Vector{Float64}} # output of length nvar
    map_edges::Dict{Location, Dict{Location,Vector{EdgeABCEuclideanDistanceAutomaton}}}
    constants::Dict{Symbol,Float64}
    map_var_model_idx::Dict{VariableModel,Int} # of dim d (of a model)
    ϵ::Float64
end

# A push! method implementend by myself because of preallocation of edge_candidates
@everywhere function _push_edge!(edge_candidates::Vector{<:EdgeABCEuclideanDistanceAutomaton}, 
                                 edge::EdgeABCEuclideanDistanceAutomaton, nbr_candidates::Int)
    if nbr_candidates < length(edge_candidates)
        edge_candidates[nbr_candidates+1] = edge
    else
        push!(edge_candidates, edge)
    end
end

@everywhere function _find_edge_candidates!(edge_candidates::Vector{EdgeABCEuclideanDistanceAutomaton},
                                            edges_from_current_loc::Dict{Location,Vector{EdgeABCEuclideanDistanceAutomaton}},
                                            Λ::Dict{Location,InvariantPredicateFunction},
                                            S_time::Float64, S_values::Vector{Float64},
                                            x::Vector{Int}, p::Vector{Float64}, ϵ::Float64,
                                            only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        if !Λ[target_loc](x) continue end
        for edge in edges_from_current_loc[target_loc]
            if edge.check_constraints(S_time, S_values, x, p, ϵ)
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

@everywhere function _get_edge_index(edge_candidates::Vector{EdgeABCEuclideanDistanceAutomaton}, nbr_candidates::Int,
                                     detected_event::Bool, tr_nplus1::Transition)
    ind_edge = 0
    bool_event = detected_event
    for i = 1:nbr_candidates
        edge = edge_candidates[i]
        # Asynchronous edge detection: we fire it
        if edge.transitions == nothing
            return (i, detected_event)
        end
        # Synchronous detection
        if !detected_event && tr_nplus1 != nothing
            if (edge.transitions[1] == :ALL) || (tr_nplus1 in edge.transitions)
                ind_edge = i
                bool_event = true
            end
        end
    end
    return (ind_edge, bool_event)
end

@everywhere function next_state!(A::ABCEuclideanDistanceAutomaton,
                                 ptr_loc_state::Vector{Symbol}, values_state::Vector{Float64}, ptr_time_state::Vector{Float64},
                                 xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                                 xn::Vector{Int}, p::Vector{Float64},
                                 edge_candidates::Vector{EdgeABCEuclideanDistanceAutomaton}; verbose::Bool = false)
    # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
    detected_event::Bool = false
    turns = 0
    Λ = getfield(A, :Λ)
    flow = getfield(A, :flow)
    map_edges = getfield(A, :map_edges)
    ϵ = A.ϵ
    if verbose 
        println("##### Begin next_state!")
        @show xnplus1, tnplus1, tr_nplus1
    end
    # First, we check the asynchronous transitions
    while true
        turns += 1
        #edge_candidates = empty!(edge_candidates) 
        edges_from_current_loc = map_edges[ptr_loc_state[1]]
        # Save all edges that satisfies transition predicate (asynchronous ones)
        nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                ptr_time_state[1], values_state, xn, p, ϵ, true)
        # Search the one we must chose, here the event is nothing because 
        # we're not processing yet the next event
        ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, nothing)
        # Update the state with the chosen one (if it exists)
        # Should be xn here
        #first_round = false
        if ind_edge > 0
            firing_edge = edge_candidates[ind_edge]
            ptr_loc_state[1] = firing_edge.update_state!(ptr_time_state[1], values_state, xn, p, ϵ)
        else
            if verbose println("No edge fired") end
            break 
        end
        if verbose
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event, nbr_candidates
            @show ptr_loc_state[1]
            @show ptr_time_state[1]
            @show values_state
            if turns == 500
                @warn "We've reached 500 turns"
            end
        end
        # For debug
        #=
        if turns > 100
        println("Number of turns in next_state! is suspicious")
        @show first_round, detected_event
        @show length(edge_candidates)
        )@show tnplus1, tr_nplus1, xnplus1
        @show edge_candidates
        error("Unpredicted behavior automaton")
        end
        =#
    end
    if verbose 
        println("Time flies with the flow...")
    end
    # Now time flies according to the flow
    for i in eachindex(values_state)
        coeff_deriv = flow[ptr_loc_state[1]][i]
        if coeff_deriv > 0
            values_state[i] += coeff_deriv*(tnplus1 - ptr_time_state[1])
        end
    end
    ptr_time_state[1] = tnplus1
    if verbose 
        @show ptr_loc_state[1]
        @show ptr_time_state[1]
        @show values_state
    end
    # Now firing an edge according to the event 
    while true
        turns += 1
        edges_from_current_loc = map_edges[ptr_loc_state[1]]
        # Save all edges that satisfies transition predicate (synchronous ones)
        nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                ptr_time_state[1], values_state, xnplus1, p, ϵ, false)
        # Search the one we must chose
        ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, tr_nplus1)
        # Update the state with the chosen one (if it exists)
        if ind_edge > 0
            firing_edge = edge_candidates[ind_edge]
            ptr_loc_state[1] = firing_edge.update_state!(ptr_time_state[1], values_state, xnplus1, p, ϵ)
        end
        if ind_edge == 0 || detected_event
            if verbose 
                if detected_event    
                    println("Synchronized with $(tr_nplus1)") 
                    @show turns
                    @show edge_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show detected_event
                    @show ptr_loc_state[1]
                    @show ptr_time_state[1]
                    @show values_state
                else
                    println("No edge fired")
                end
            end
            break 
        end
        if verbose
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event, nbr_candidates
            @show detected_event
            @show ptr_loc_state[1]
            @show ptr_time_state[1]
            @show values_state
            if turns == 500
                @warn "We've reached 500 turns"
            end
        end
        # For debug
        #=
        if turns > 100
        println("Number of turns in next_state! is suspicious")
        @show detected_event
        @show length(edge_candidates)
        @show tnplus1, tr_nplus1, xnplus1
        @show edge_candidates
        error("Unpredicted behavior automaton")
        end
        =#
    end
    if verbose 
        println("##### End next_state!") 
    end
end

function create_abc_euclidean_distance_automaton(m::ContinuousTimeModel, timeline::AbstractVector{Float64}, observations::AbstractVector{Float64}, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert length(timeline) == length(observations) "Timeline and observations vectors don't have the same length"
    nbr_observations = length(observations)

    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :ABCEuclideanDistanceAutomaton
    edge_type = :EdgeABCEuclideanDistanceAutomaton
    check_constraints = Symbol("check_constraints_$(lha_name)")
    update_state! = Symbol("update_state_$(lha_name)!")

    # Locations
    locations = [:l0, :l1, :l2]

    ## Invariant predicates
    @everywhere true_inv_predicate(x::Vector{Int}) = true 
    Λ_F = Dict{Location,InvariantPredicateFunction}(:l0 => getfield(Main, :true_inv_predicate), :l1 => getfield(Main, :true_inv_predicate),
                                                    :l2 => getfield(Main, :true_inv_predicate))

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
    idx_obs_var = m.map_var_idx[sym_obs]
    to_idx(var::Symbol) = map_var_automaton_idx[var]

    id = MarkovProcesses.newid()
    function check_constraints(from_loc::Location, to_loc::Location, edge_number::Int)
        return Symbol("check_constraints_$(edge_type)_$(from_loc)$(to_loc)_$(edge_number)_$(model_name)_$(id)")
    end
    function update_state!(from_loc::Location, to_loc::Location, edge_number::Int)
        return Symbol("update_state_$(edge_type)_$(from_loc)$(to_loc)_$(edge_number)_$(model_name)_$(id)!")
    end

    ## check_constraints & update_state!
    meta_funcs = quote
        # l0 loc
        # l0 => l1
        #struct $(edge_name(:l0, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = true
        @everywhere $(update_state!(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:d))] = 0.0;
         S_values[$(to_idx(:idx))] = 1.0;
         :l1)

        # l1 loc
        # l1 => l1
        # Defined below 
        #struct $(edge_name(:l1, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) =
        (tml = $(Tuple(timeline));
         tml_idx = tml[convert(Int, S_values[$(to_idx(:idx))])];
         S_values[$(to_idx(:t))] >= tml_idx)
        @everywhere $(update_state!(:l1, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) =
        (y_obs = $(Tuple(observations));
         y_obs_idx = y_obs[convert(Int, S_values[$(to_idx(:idx))])];
         S_values[$(to_idx(:d))] += (S_values[$(to_idx(:n))]-y_obs_idx)^2;
         S_values[$(to_idx(:idx))] += 1.0;
         :l1)

        #struct $(edge_name(:l1, :l1, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l1, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = true 
        @everywhere $(update_state!(:l1, :l1, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         :l1)

        # l1 => l2
        #struct $(edge_name(:l1, :l2, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        S_values[$(to_idx(:idx))] >= ($nbr_observations + 1)
        @everywhere $(update_state!(:l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        (S_values[$(to_idx(:d))] = sqrt(S_values[$(to_idx(:d))]);
         :l2)

        @everywhere $(check_constraints(:l1, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        (S_values[$(to_idx(:d))] > ϵ^2)
        @everywhere $(update_state!(:l1, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}, ϵ::Float64) = 
        (S_values[$(to_idx(:d))] = Inf;
         :l2)
    end
    eval(meta_funcs)

    @eval begin
        map_edges = Dict{Location,Dict{Location,Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location,Vector{$(edge_type)}}()
        end

        ## Edges
        # l0 loc
        # l0 => l1
        edge1 = EdgeABCEuclideanDistanceAutomaton(nothing, $(check_constraints(:l0, :l1, 1)), $(update_state!(:l0, :l1, 1)))
        map_edges[:l0][:l1] = [edge1]

        # l1 loc
        # l1 => l1
        edge1 = EdgeABCEuclideanDistanceAutomaton(nothing, $(check_constraints(:l1, :l1, 1)), $(update_state!(:l1, :l1, 1)))
        edge2 = EdgeABCEuclideanDistanceAutomaton([:ALL], $(check_constraints(:l1, :l1, 2)), $(update_state!(:l1, :l1, 2)))
        map_edges[:l1][:l1] = [edge1, edge2]

        # l1 => l2
        edge1 = EdgeABCEuclideanDistanceAutomaton(nothing, $(check_constraints(:l1, :l2, 1)), $(update_state!(:l1, :l2, 1)))
        edge2 = EdgeABCEuclideanDistanceAutomaton(nothing, $(check_constraints(:l1, :l2, 2)), $(update_state!(:l1, :l2, 2)))
        map_edges[:l1][:l2] = [edge1,edge2]
    end

    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)
    for i = 1:nbr_observations
        constants[Symbol("tml_$(convert(Float64, i))")] = timeline[i]
        constants[Symbol("y_$(convert(Float64, i))")] = observations[i]
    end

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = ABCEuclideanDistanceAutomaton(m.transitions, locations, Λ_F, locations_init, locations_final, 
                                      map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx, Inf)

    return A
end

