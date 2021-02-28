
# Creation of the automaton types
#@everywhere @eval abstract type EdgeEuclideanDistanceAutomaton <: Edge end
@everywhere struct EdgeEuclideanDistanceAutomaton <: Edge
    transitions::TransitionSet
    check_constraints::CheckConstraintsFunction
    update_state!::UpdateStateFunction
end
@everywhere @eval $(MarkovProcesses.generate_code_lha_type_def(:EuclideanDistanceAutomaton, :EdgeEuclideanDistanceAutomaton))

function create_euclidean_distance_automaton(m::ContinuousTimeModel, timeline::AbstractVector{Float64}, observations::AbstractVector{Float64}, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert length(timeline) == length(observations) "Timeline and observations vectors don't have the same length"
    nbr_observations = length(observations)

    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :EuclideanDistanceAutomaton
    edge_type = :EdgeEuclideanDistanceAutomaton
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
    basename_func = "$(model_name)_$(id)"
    edge_name(from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("Edge_$(lha_name)_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)")
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
        @everywhere $(check_constraints(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(update_state!(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:d))] = 0.0;
         S_values[$(to_idx(:idx))] = 1.0;
         :l1)

        # l1 loc
        # l1 => l1
        # Defined below 
        #struct $(edge_name(:l1, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (tml = $(SVector{length(timeline)}(timeline));
         tml_idx = tml[convert(Int, S_values[$(to_idx(:idx))])];
         S_values[$(to_idx(:t))] >= tml_idx)
        @everywhere $(update_state!(:l1, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (y_obs = $(SVector{length(observations)}(observations));
         y_obs_idx = y_obs[convert(Int, S_values[$(to_idx(:idx))])];
         S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))]+(S_values[$(to_idx(:n))]-y_obs_idx)^2;
         S_values[$(to_idx(:idx))] = S_values[$(to_idx(:idx))]+1.0;
         :l1)

        #struct $(edge_name(:l1, :l1, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l1, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
        @everywhere $(update_state!(:l1, :l1, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         :l1)

        # l1 => l2
        #struct $(edge_name(:l1, :l2, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:idx))] >= ($nbr_observations + 1)
        @everywhere $(update_state!(:l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = sqrt(S_values[$(to_idx(:d))]);
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
        edge1 = EdgeEuclideanDistanceAutomaton(nothing, $(check_constraints(:l0, :l1, 1)), $(update_state!(:l0, :l1, 1)))
        map_edges[:l0][:l1] = [edge1]

        # l1 loc
        # l1 => l1
        edge1 = EdgeEuclideanDistanceAutomaton(nothing, $(check_constraints(:l1, :l1, 1)), $(update_state!(:l1, :l1, 1)))
        edge2 = EdgeEuclideanDistanceAutomaton([:ALL], $(check_constraints(:l1, :l1, 2)), $(update_state!(:l1, :l1, 2)))
        map_edges[:l1][:l1] = [edge1, edge2]

        # l1 => l2
        edge1 = EdgeEuclideanDistanceAutomaton(nothing, $(check_constraints(:l1, :l2, 1)), $(update_state!(:l1, :l2, 1)))
        map_edges[:l1][:l2] = [edge1]
    end

    map_edges_transitions = Dict{Symbol, Dict{Symbol,Vector{TransitionSet}}}()
    map_edges_check_constraints = Dict{Symbol, Dict{Symbol,Vector{CheckConstraintsFunction}}}()
    map_edges_update_state = Dict{Symbol, Dict{Symbol,Vector{UpdateStateFunction}}}()
    
    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)
    for i = 1:nbr_observations
        constants[Symbol("tml_$(convert(Float64, i))")] = timeline[i]
        constants[Symbol("y_$(convert(Float64, i))")] = observations[i]
    end

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = EuclideanDistanceAutomaton(m.transitions, locations, Λ_F, locations_init, locations_final, 
                                   map_var_automaton_idx, flow, 
                                   map_edges, map_edges_transitions, map_edges_check_constraints, map_edges_update_state,
                                   constants, m.map_var_idx)

    return A
end

