
# l0 loc
# l0 => l1
@everywhere cc_eucl_dist_aut_l0l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l1 loc
# l1 => l1
@everywhere cc_eucl_dist_aut_l1l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:t] >= constants[Symbol("tml_$(S[:idx])")]
@everywhere us_eucl_dist_aut_l1l1_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S[:d] += (S[:n]  - constants[Symbol("y_$(S[:idx])")])^2;
 S[:idx] += 1.0)

# l1 => l2
@everywhere cc_eucl_dist_aut_l1l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:idx] >= constants[:nbr_obs] + 1
@everywhere us_eucl_dist_aut_l1l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2; 
 S[:d] = sqrt(S[:d]))

@everywhere cc_eucl_dist_aut_l1l1_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true 

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
    nbr_first_vars = 4
    for i = 1:nbr_observations
        map_var_automaton_idx[Symbol("tml_$(convert(Float64, i))")] = nbr_first_vars + i
    end
    for i = 1:nbr_observations
        map_var_automaton_idx[Symbol("y_$(convert(Float64, i))")] = nbr_first_vars + nbr_observations + i
    end
    
    vector_flow = zeros(nbr_first_vars + 2*nbr_observations)
    vector_flow[1] = 1.0
    flow = Dict{Location,Vector{Float64}}(:l0 => vector_flow, 
                                          :l1 => vector_flow, 
                                          :l2 => vector_flow)

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end

    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    
    # l0 loc
    # l0 => l1
    sym_func_us_l0l1_1 = Symbol("us_eucl_dist_$(basename_func)_l0l1_1!")
    str_us_l0l1_1 = "
    @everywhere $(sym_func_us_l0l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
        (S.loc = :l1; \n
         S[:n] = x[$(idx_obs_var)];\n
         S[:d] = 0.0;\n
         S[:idx] = 1.0)"
    eval(Meta.parse(str_us_l0l1_1))
    edge1 = Edge([nothing], getfield(Main, :cc_eucl_dist_aut_l0l1_1), getfield(Main, sym_func_us_l0l1_1))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l1
    sym_func_us_l1l2_2 = Symbol("us_eucl_dist_$(basename_func)_l1l2_2!")
    str_us_l1l2_2 = "
    @everywhere $(sym_func_us_l1l2_2)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =\n 
    (S[:n] = x[$(idx_obs_var)])"
    eval(Meta.parse(str_us_l1l2_2))
    edge1 = Edge([nothing], getfield(Main, :cc_eucl_dist_aut_l1l1_1), getfield(Main, :us_eucl_dist_aut_l1l1_1!))
    edge2 = Edge([:ALL], getfield(Main, :cc_eucl_dist_aut_l1l1_2), getfield(Main, sym_func_us_l1l2_2))
    map_edges[:l1][:l1] = [edge1, edge2]

    # l1 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_eucl_dist_aut_l1l2_1), getfield(Main, :us_eucl_dist_aut_l1l2_1!))
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

