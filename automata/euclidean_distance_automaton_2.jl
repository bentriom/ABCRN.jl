
# Creation of the automaton types
#@everywhere @eval abstract type EdgeEuclideanDistanceAutomaton2 <: Edge end
@everywhere struct EdgeEuclideanDistanceAutomaton2 <: Edge
    transitions::TransitionSet
    check_constraints::CheckConstraintsFunction
    update_state!::UpdateStateFunction
end
@everywhere @eval $(MarkovProcesses.generate_code_lha_type_def(:EuclideanDistanceAutomaton2,:EdgeEuclideanDistanceAutomaton2))

function create_euclidean_distance_automaton_2(m::ContinuousTimeModel, timeline::AbstractVector{Float64}, observations::AbstractVector{Float64}, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert length(timeline) == length(observations) "Timeline and observations vectors don't have the same length"
    nbr_observations = length(observations)

    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :EuclideanDistanceAutomaton2
    edge_type = :EdgeEuclideanDistanceAutomaton2
    check_constraints = Symbol("check_constraints_$(lha_name)")
    update_state! = Symbol("update_state_$(lha_name)!")

    # Locations
    locations = [:l0, :lfinal]
    for i = 1:nbr_observations
        push!(locations, Symbol("l$(i)"))
    end

    ## Invariant predicates
    @everywhere true_inv_predicate(x::Vector{Int}) = true
    Λ_F = Dict{Location, Function}()
    for loc in locations
        Λ_F[loc] = getfield(Main, :true_inv_predicate)
    end

    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:lfinal]

    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:t => 1, :n => 2, :d => 3)
    vector_flow = [1.0, 0.0, 0.0]
    flow = Dict{Location, Vector{Float64}}()
    for loc in locations
        flow[loc] = vector_flow
    end

    ## Edges
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
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

    loc_nbr_obs = Symbol("l$(nbr_observations)")
    meta_funcs = quote
        # l0 loc
        # l0 => l1
        #struct $(edge_name(:l0, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(update_state!(:l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:d))] = 0.0;
         :l1)

        # lnbr_obs => lfinal
        #struct $(edge_name(loc_nbr_obs, :lfinal, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(loc_nbr_obs, :lfinal, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $(timeline[nbr_observations])
        @everywhere $(update_state!(loc_nbr_obs, :lfinal, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))]+(S_values[$(to_idx(:n))]-$(observations[nbr_observations]))^2;
         S_values[$(to_idx(:d))] = sqrt(S_values[$(to_idx(:d))]);
         :lfinal)

        # lnbr_obs => lnbr_obs
        #struct $(edge_name(loc_nbr_obs, loc_nbr_obs, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(loc_nbr_obs, loc_nbr_obs, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
        @everywhere $(update_state!(loc_nbr_obs, loc_nbr_obs, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         $(Meta.quot(loc_nbr_obs)))
    end
    eval(meta_funcs)

    @eval begin
        map_edges = Dict{Location,Dict{Location,Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location,Vector{$(edge_type)}}()
        end

        # l0 loc
        # l0 => l1
        edge1 = EdgeEuclideanDistanceAutomaton2(nothing, $(check_constraints(:l0, :l1, 1)), $(update_state!(:l0, :l1, 1)))
        map_edges[:l0][:l1] = [edge1]

        # lnbr_obs => lfinal
        edge1 = EdgeEuclideanDistanceAutomaton2(nothing, $(check_constraints(loc_nbr_obs, :lfinal, 1)), $(update_state!(loc_nbr_obs, :lfinal, 1)))
        map_edges[$(Meta.quot(loc_nbr_obs))][:lfinal] = [edge1]
        # lnbr_obs => lnbr_obs
        edge1 = EdgeEuclideanDistanceAutomaton2([:ALL], $(check_constraints(loc_nbr_obs, loc_nbr_obs, 1)), $(update_state!(loc_nbr_obs, loc_nbr_obs, 1)))
        map_edges[$(Meta.quot(loc_nbr_obs))][$(Meta.quot(loc_nbr_obs))] = [edge1]
    end

    for i = 1:(nbr_observations-1)
        loci = Symbol("l$(i)")
        locip1 = Symbol("l$(i+1)")
        meta_funcs = quote
            #struct $(edge_name(loci, locip1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
            @everywhere $(check_constraints(loci, locip1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            S_values[$(to_idx(:t))] >= $(timeline[i])
            @everywhere $(update_state!(loci, locip1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))]+(S_values[$(to_idx(:n))]-$(observations[i]))^2;
             $(Meta.quot(locip1)))

            #struct $(edge_name(loci, loci, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
            @everywhere $(check_constraints(loci, loci, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
            @everywhere $(update_state!(loci, loci, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
            (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
             $(Meta.quot(loci)))
        end
        eval(meta_funcs)

        # loci => loci+1
        edge1 = EdgeEuclideanDistanceAutomaton2(nothing, getfield(Main, check_constraints(loci, locip1, 1)), getfield(Main, update_state!(loci, locip1, 1)))
        map_edges[loci][locip1] = [edge1]
        # loci => loci
        edge1 = EdgeEuclideanDistanceAutomaton2([:ALL], getfield(Main, check_constraints(loci, loci, 1)), getfield(Main, update_state!(loci, loci, 1)))
        map_edges[loci][loci] = [edge1]
    end

    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)

    map_edges_transitions = Dict{Symbol, Dict{Symbol,Vector{TransitionSet}}}()
    map_edges_check_constraints = Dict{Symbol, Dict{Symbol,Vector{CheckConstraintsFunction}}}()
    map_edges_update_state = Dict{Symbol, Dict{Symbol,Vector{UpdateStateFunction}}}()

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = EuclideanDistanceAutomaton2(m.transitions, locations, Λ_F, locations_init, locations_final, 
                                    map_var_automaton_idx, flow, 
                                    map_edges, map_edges_transitions, map_edges_check_constraints, map_edges_update_state,
                                    constants, m.map_var_idx)
    return A
end

