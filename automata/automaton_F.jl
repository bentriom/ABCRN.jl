
@everywhere istrue(val::Float64) = convert(Bool, val)

# Invariant predicates functions
@everywhere true_inv_predicate(x::Vector{Int}) = true 

# Check constraints and update state functions

# l0 loc : we construct  the edges of the form l0 => (..)
# "cc" as check_constraints and "us" as update_state
# l0 => l1
@everywhere cc_aut_F_l0l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true
# us_aut_F_l0l1_1! inside create_automaton_F

# l1 loc
# l1 => l2
@everywhere cc_aut_F_l1l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) >= constants[:t1] &&
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_F_l1l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2;
 S[:d] = 0)

@everywhere cc_aut_F_l1l2_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) >= constants[:t1] &&
S[:d] == 0 
@everywhere us_aut_F_l1l2_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2)

@everywhere cc_aut_F_l1l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t2]) && 
(S[:n] < constants[:x1] || S[:n] > constants[:x2])
@everywhere us_aut_F_l1l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2;
 S[:d] = min(abs(S[:n] - constants[:x1]), abs(S[:n] - constants[:x2])))

@everywhere cc_aut_F_l1l2_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && getfield(S, :time) <= constants[:t2]
@everywhere us_aut_F_l1l2_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2)

# l1 => l3
@everywhere cc_aut_F_l1l3_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_F_l1l3_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3;
 S[:d] = 0;)

@everywhere cc_aut_F_l1l3_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S[:n] < constants[:x1] || S[:n] > constants[:x2]) && 
(getfield(S, :time) <= constants[:t1])
@everywhere us_aut_F_l1l3_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3;
 S[:d] = min(sqrt((getfield(S, :time) - constants[:t1])^2 + (S[:n] - constants[:x2])^2), 
             sqrt((getfield(S, :time) - constants[:t1])^2 + (S[:n] - constants[:x1])^2)))

@everywhere cc_aut_F_l1l3_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S[:n] < constants[:x1] || S[:n] > constants[:x2]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2])
@everywhere us_aut_F_l1l3_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3;
 S[:d] = min(S[:d], min(abs(S[:n] - constants[:x1]), abs(S[:n] - constants[:x2]))))

# l3 loc
# l3 => l1
@everywhere cc_aut_F_l3l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l3 => l2
@everywhere cc_aut_F_l3l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t2])
@everywhere us_aut_F_l3l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2)

function create_automaton_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    # Locations
    locations = [:l0, :l1, :l2, :l3]

    ## Invariant predicates
    Λ_F = Dict(:l0 => getfield(Main, :true_inv_predicate), :l1 => getfield(Main, :true_inv_predicate),
               :l2 => getfield(Main, :true_inv_predicate), :l3 => getfield(Main, :true_inv_predicate))
    
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
    
    sym_isabs_func = Symbol(m.isabsorbing)
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')

    # l0 loc
    # l0 => l1
    sym_func_us_l0l1_1 = Symbol("us_aut_F_$(basename_func)_l0l1_1!")
    str_us_l0l1_1 = "
    @everywhere $(sym_func_us_l0l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1; \n
     S[:n] = x[$(idx_obs_var)];\n
     S[:d] = Inf; \n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l0l1_1))
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l0l1_1), getfield(Main, sym_func_us_l0l1_1))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l1l2_1), getfield(Main, :us_aut_F_l1l2_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_F_l1l2_2), getfield(Main, :us_aut_F_l1l2_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_F_l1l2_3), getfield(Main, :us_aut_F_l1l2_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_F_l1l2_4), getfield(Main, :us_aut_F_l1l2_4!))
    map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]

    # l1 => l3
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l1l3_1), getfield(Main, :us_aut_F_l1l3_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_F_l1l3_2), getfield(Main, :us_aut_F_l1l3_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_F_l1l3_3), getfield(Main, :us_aut_F_l1l3_3!))
    map_edges[:l1][:l3] = [edge1, edge2, edge3]
  
    # l3 loc
    # l3 => l1
    sym_func_us_l3l1_1 = Symbol("us_aut_F_$(basename_func)_l0l1_1!")
    str_us_l3l1_1 = 
    "@everywhere $(sym_func_us_l3l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1;\n
    S[:n] = x[$(idx_obs_var)];\n
    S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l3l1_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_F_l3l1_1), getfield(Main, sym_func_us_l3l1_1))
    map_edges[:l3][:l1] = [edge1]
    
    # l3 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l3l2_1), getfield(Main, :us_aut_F_l3l2_1!))
    map_edges[:l3][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA("F property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_F

