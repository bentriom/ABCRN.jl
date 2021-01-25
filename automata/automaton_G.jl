
@everywhere istrue(val::Float64) = convert(Bool, val)

# Invariant predicate functions
@everywhere true_inv_predicate(x::Vector{Int}) = true 

# l0 loc
# l0 => l1
@everywhere cc_aut_G_l0l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l1 => l3
@everywhere cc_aut_G_l1l3_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) <= constants[:t1] && 
S[:n] < constants[:x1] || S[:n] > constants[:x2]
@everywhere us_aut_G_l1l3_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3; 
 S[:d] = min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
 S[:in] = false)

@everywhere cc_aut_G_l1l3_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) <= constants[:t1]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1l3_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3; 
 S[:d] = 0; 
 S[:in] = false)

@everywhere cc_aut_G_l1l3_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1l3_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3; 
 S[:d] = S[:d] * (getfield(S, :time) - constants[:t1]); 
 S[:tprime] = 0.0)

@everywhere cc_aut_G_l1l3_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1l3_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3; 
 S[:tprime] = 0.0)

# l1 => l4
@everywhere cc_aut_G_l1l4_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(S[:n] < constants[:x1] || S[:n] > constants[:x2])
@everywhere us_aut_G_l1l4_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l4; 
 S[:d] += S[:d] * (getfield(S, :time) - constants[:t1]))

@everywhere cc_aut_G_l1l4_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(S[:n] < constants[:x1] || S[:n] > constants[:x2])
@everywhere us_aut_G_l1l4_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l4)


# l1 => l2
@everywhere cc_aut_G_l1l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l1l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2)

@everywhere cc_aut_G_l1l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l1l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2; 
 S[:d] = S[:d] * (constants[:t2] - constants[:t1]))

@everywhere cc_aut_G_l1l2_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && 
getfield(S, :time) <= constants[:t1]
@everywhere us_aut_G_l1l2_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2; 
 S[:d] = (constants[:t2] - constants[:t1]) *
min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))

@everywhere cc_aut_G_l1l2_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2])
@everywhere us_aut_G_l1l2_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2; 
 S[:d] += (constants[:t2] - getfield(S, :time)) * 
min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))

# l3 => l1
@everywhere cc_aut_G_l3l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l4 => l1
@everywhere cc_aut_G_l4l1_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l2 => l1
@everywhere cc_aut_G_l3l2_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l3l2_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2;
 S[:d] = S[:d] * (constants[:t2] - constants[:t1]))

@everywhere cc_aut_G_l3l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l3l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2)

# l4 => l2
@everywhere cc_aut_G_l4l2_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t2])
@everywhere us_aut_G_l4l2_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2; 
 S[:d] +=  S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
 S[:tprime] = 0.0)

function create_automaton_G(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    @assert sym_obs in m.g
    # Locations
    locations = [:l0, :l1, :l2, :l3, :l4]

    # Invariant predicates
    Λ_F = Dict(:l0 => getfield(Main, :true_inv_predicate), :l1 => getfield(Main, :true_inv_predicate),
               :l2 => getfield(Main, :true_inv_predicate), :l3 => getfield(Main, :true_inv_predicate), 
               :l4 => getfield(Main, :true_inv_predicate))
    
    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:l2]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:tprime => 1, :in => 2,
                                                        :n => 3,  :d => 4, :isabs => 5)

    ## Flow of variables
    flow = Dict{Location,Vector{Float64}}(:l0 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l1 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l2 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l3 => [0.0,0.0,0.0,0.0,0.0], 
                                          :l4 => [1.0,0.0,0.0,0.0,0.0])

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
    sym_func_us_l0l1_1 = Symbol("us_aut_G_$(basename_func)_l0l1_1!")
    str_us_l0l1_1 = "
    @everywhere $(sym_func_us_l0l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1; \n
    S[:d] = 0; \n
    S[:n] = x[$(idx_obs_var)]; \n
    S[:in] = true; \n
    S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l0l1_1))
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l0l1_1), getfield(Main, sym_func_us_l0l1_1))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l3
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1l3_1), getfield(Main, :us_aut_G_l1l3_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1l3_2), getfield(Main, :us_aut_G_l1l3_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_G_l1l3_3), getfield(Main, :us_aut_G_l1l3_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_G_l1l3_4), getfield(Main, :us_aut_G_l1l3_4!))
    map_edges[:l1][:l3] = [edge1, edge2, edge3, edge4]

    # l1 => l4
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1l4_1), getfield(Main, :us_aut_G_l1l4_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1l4_2), getfield(Main, :us_aut_G_l1l4_2!))
    map_edges[:l1][:l4] = [edge1, edge2]
   
    # l1 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1l2_1), getfield(Main, :us_aut_G_l1l2_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1l2_2), getfield(Main, :us_aut_G_l1l2_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_G_l1l2_3), getfield(Main, :us_aut_G_l1l2_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_G_l1l2_4), getfield(Main, :us_aut_G_l1l2_4!))
    map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]

    # l3 loc
    # l3 => l1
    sym_func_us_l3l1_1 = Symbol("us_aut_G_$(basename_func)_l3l1_1!")
    str_us_l3l1_1 = "
    @everywhere $(sym_func_us_l3l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
    (S.loc = :l1; 
    S[:n] = x[$(idx_obs_var)]; 
    S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l3l1_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_G_l3l1_1), getfield(Main, sym_func_us_l3l1_1))
    map_edges[:l3][:l1] = [edge1]

    # l3 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l3l2_1), getfield(Main, :us_aut_G_l3l2_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l3l2_2), getfield(Main, :us_aut_G_l3l2_2!))
    map_edges[:l3][:l2] = [edge1, edge2]

    # l4 loc
    # l4 => l1
    sym_func_us_l4l1_1 = Symbol("us_aut_G_$(basename_func)_l4l1_1!")
    str_us_l4l1_1 = "
    @everywhere $(sym_func_us_l4l1_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1; \n
    S[:d] += S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); \n
    S[:tprime] = 0.0; \n
    S[:n] = x[$(idx_obs_var)]; \n
    S[:in] = true; \n
    S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l4l1_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_G_l4l1_1), getfield(Main, sym_func_us_l4l1_1))
    map_edges[:l4][:l1] = [edge1]

    # l4 => l2
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l4l2_1), getfield(Main, :us_aut_G_l4l2_1!))
    map_edges[:l4][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA("G property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
   
end

export create_automaton_G

