
function create_euclidean_distance_automaton_2(m::ContinuousTimeModel, timeline::AbstractVector{Float64}, observations::AbstractVector{Float64}, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert length(timeline) == length(observations) "Timeline and observations vectors don't have the same length"
    nbr_observations = length(observations)

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
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end
   
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    to_idx(var::Symbol) = map_var_automaton_idx[var] 
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_eucl_dist_aut_2_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    loc_nbr_obs = Symbol("l$(nbr_observations)")
    meta_elementary_functions = quote
        # l0 loc
        # l0 => l1
        @everywhere $(func_name(:cc, :l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l0, :l1, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1"); 
         setindex!(S_values, x[$(idx_obs_var)], $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:d))))
       
        # lnbr_obs => lfinal
        @everywhere $(func_name(:cc, loc_nbr_obs, :lfinal, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $(timeline[nbr_observations])
        @everywhere $(func_name(:us, loc_nbr_obs, :lfinal, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("lfinal"); 
         setindex!(S_values, S_values[$(to_idx(:d))] + (S_values[$(to_idx(:n))]-$(observations[nbr_observations]))^2, 
                                         $(to_idx(:d)));
         setindex!(S_values, sqrt(S_values[$(to_idx(:d))]), $(to_idx(:d))))

        # lnbr_obs => lnbr_obs
        @everywhere $(func_name(:cc, loc_nbr_obs, loc_nbr_obs, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
        @everywhere $(func_name(:us, loc_nbr_obs, loc_nbr_obs, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (setindex!(S_values, x[$(idx_obs_var)], $(to_idx(:n))))
    end
    eval(meta_elementary_functions)
    # l0 loc
    # l0 => l1
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :l1, 1)), getfield(Main, func_name(:us, :l0, :l1, 1)))
    map_edges[:l0][:l1] = [edge1]

    # lnbr_obs => lfinal
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, loc_nbr_obs, :lfinal, 1)), getfield(Main, func_name(:us, loc_nbr_obs, :lfinal, 1)))
    map_edges[loc_nbr_obs][:lfinal] = [edge1]
    # lnbr_obs => lnbr_obs
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, loc_nbr_obs, loc_nbr_obs, 1)), getfield(Main, func_name(:us, loc_nbr_obs, loc_nbr_obs, 1)))
    map_edges[loc_nbr_obs][loc_nbr_obs] = [edge1]

    for i = 1:(nbr_observations-1)
        loci = Symbol("l$(i)")
        locip1 = Symbol("l$(i+1)")
        meta_elementary_functions_loci = quote
            # l1 loc
            # l1 => l1
            # Defined below 
            @everywhere $(func_name(:cc, loci, locip1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            S_values[$(to_idx(:t))] >= $(timeline[i])
            @everywhere $(func_name(:us, loci, locip1, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
            (ptr_loc[1] = $(Meta.quot(locip1)); 
             setindex!(S_values, S_values[$(to_idx(:d))] + (S_values[$(to_idx(:n))]-$(observations[i]))^2, 
                                             $(to_idx(:d))))

            @everywhere $(func_name(:cc, loci, loci, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true 
            @everywhere $(func_name(:us, loci, loci, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
            (setindex!(S_values, x[$(idx_obs_var)], $(to_idx(:n))))
        end
        eval(meta_elementary_functions_loci)

        # loci => loci+1
        edge1 = Edge(nothing, getfield(Main, func_name(:cc, loci, locip1, 1)), getfield(Main, func_name(:us, loci, locip1, 1)))
        map_edges[loci][locip1] = [edge1]
        # loci => loci
        edge1 = Edge([:ALL], getfield(Main, func_name(:cc, loci, loci, 1)), getfield(Main, func_name(:us, loci, loci, 1)))
        map_edges[loci][loci] = [edge1]
    end
    
    ## Constants
    constants = Dict{Symbol,Float64}(:nbr_obs => nbr_observations)

    A = LHA("Euclidean distance", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_euclidean_distance_automaton_2

