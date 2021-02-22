
# Creation of the automaton types
@everywhere @eval abstract type EdgeEuclideanDistanceAutomaton2 <: Edge end
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

    loc_nbr_obs = Symbol("l$(nbr_observations)")
    @everywhere @eval begin
        # l0 loc
        # l0 => l1
        @everywhere struct $(edge_name(:l0, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l0, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(update_state!)(edge::$(edge_name(:l0, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:d))] = 0.0;
         :l1)

        # lnbr_obs => lfinal
        @everywhere struct $(edge_name(loc_nbr_obs, :lfinal, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(loc_nbr_obs, :lfinal, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $(timeline[nbr_observations])
        @everywhere $(update_state!)(edge::$(edge_name(loc_nbr_obs, :lfinal, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))]+(S_values[$(to_idx(:n))]-$(observations[nbr_observations]))^2;
         S_values[$(to_idx(:d))] = sqrt(S_values[$(to_idx(:d))]);
         :lfinal)

        # lnbr_obs => lnbr_obs
        @everywhere struct $(edge_name(loc_nbr_obs, loc_nbr_obs, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(loc_nbr_obs, loc_nbr_obs, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
        @everywhere $(update_state!)(edge::$(edge_name(loc_nbr_obs, loc_nbr_obs, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         $(Meta.quot(loc_nbr_obs)))
    end

    @eval begin
        map_edges = Dict{Location,Dict{Location,Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location,Vector{$(edge_type)}}()
        end
        
        # l0 loc
        # l0 => l1
        edge1 = $(edge_name(:l0, :l1, 1))(nothing)
        map_edges[:l0][:l1] = [edge1]

        # lnbr_obs => lfinal
        edge1 = $(edge_name(loc_nbr_obs, :lfinal, 1))(nothing)
        map_edges[$(Meta.quot(loc_nbr_obs))][:lfinal] = [edge1]
        # lnbr_obs => lnbr_obs
        edge1 = $(edge_name(loc_nbr_obs, loc_nbr_obs, 1))([:ALL])
        map_edges[$(Meta.quot(loc_nbr_obs))][$(Meta.quot(loc_nbr_obs))] = [edge1]
    end

    function generate_code_loci_functions(i::Int, loci::Symbol, locip1::Symbol)
        return quote
            struct $(edge_name(loci, locip1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
            $(check_constraints)(edge::$(edge_name(loci, locip1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            S_values[$(to_idx(:t))] >= $(timeline[i])
            $(update_state!)(edge::$(edge_name(loci, locip1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))]+(S_values[$(to_idx(:n))]-$(observations[i]))^2;
             $(Meta.quot(locip1)))

            struct $(edge_name(loci, loci, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
            $(check_constraints)(edge::$(edge_name(loci, loci, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
            $(update_state!)(edge::$(edge_name(loci, loci, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
            (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
             $(Meta.quot(loci)))
        end
    end

    function generate_code_loci_edges(loci::Symbol, locip1::Symbol)
        return quote
            # loci => loci+1
            edge1 = $(edge_name(loci, locip1, 1))(nothing)
            map_edges[$(Meta.quot(loci))][$(Meta.quot(locip1))] = [edge1]
            # loci => loci
            edge1 = $(edge_name(loci, loci, 1))([:ALL])
            map_edges[$(Meta.quot(loci))][$(Meta.quot(loci))] = [edge1]
        end
    end

    for i = 1:(nbr_observations-1)
        loci = Symbol("l$(i)")
        locip1 = Symbol("l$(i+1)")
        @everywhere @eval $(generate_code_loci_functions(i, loci, locip1))
        @everywhere @eval $(generate_code_loci_edges(loci, locip1))
    end
    
    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type, check_constraints, update_state!))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = EuclideanDistanceAutomaton2(m.transitions, locations, Λ_F, locations_init, locations_final, 
                                    map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

