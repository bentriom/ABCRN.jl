
@everywhere istrue(val::Float64) = convert(Bool, val)

## Invariant predicate functions

@everywhere true_inv_predicate(x::Vector{Int}) = true 

## Edges check constraint and update state functions

# l0G loc
# l0G => l1G
@everywhere cc_aut_G_l0Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l1G loc
# l1G => l3G
@everywhere cc_aut_G_l1Gl3G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) <= constants[:t1] && 
S[:n] < constants[:x1] || S[:n] > constants[:x2]
@everywhere us_aut_G_l1Gl3G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3G; 
 S[:d] = min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
 S[:in] = false)

@everywhere cc_aut_G_l1Gl3G_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) <= constants[:t1]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1Gl3G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3G; 
 S[:d] = 0; 
 S[:in] = false)

@everywhere cc_aut_G_l1Gl3G_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1Gl3G_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3G; 
 S[:d] = S[:d] * (getfield(S, :time) - constants[:t1]); 
 S[:tprime] = 0.0)

@everywhere cc_aut_G_l1Gl3G_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(constants[:x1] <= S[:n] <= constants[:x2])
@everywhere us_aut_G_l1Gl3G_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3G; 
 S[:tprime] = 0.0)

# l1G => l4G
@everywhere cc_aut_G_l1Gl4G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(S[:n] < constants[:x1] || S[:n] > constants[:x2])
@everywhere us_aut_G_l1Gl4G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l4G; 
 S[:d] += S[:d] * (getfield(S, :time) - constants[:t1]))

@everywhere cc_aut_G_l1Gl4G_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2]) && 
(S[:n] < constants[:x1] || S[:n] > constants[:x2])
@everywhere us_aut_G_l1Gl4G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l4G)

# l1G => l2G
@everywhere cc_aut_G_l1Gl2G_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && 
getfield(S, :time) <= constants[:t1]
@everywhere us_aut_G_l1Gl2G_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G; 
 S[:d] = (constants[:t2] - constants[:t1]) *
 min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))

@everywhere cc_aut_G_l1Gl2G_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && 
(constants[:t1] <= getfield(S, :time) <= constants[:t2])
@everywhere us_aut_G_l1Gl2G_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G; 
 S[:d] += (constants[:t2] - getfield(S, :time)) * 
 min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])))

@everywhere cc_aut_G_l1Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l1Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G)

@everywhere cc_aut_G_l1Gl2G_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l1Gl2G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G; 
 S[:d] = S[:d] * (constants[:t2] - constants[:t1]))

# l3G loc
# l3G => l1G
@everywhere cc_aut_G_l3Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l3G => l2G
@everywhere cc_aut_G_l3Gl2G_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l3Gl2G_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G;
 S[:d] = S[:d] * (constants[:t2] - constants[:t1]))

@everywhere cc_aut_G_l3Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
!istrue(S[:in]) && 
getfield(S, :time) >= constants[:t2]
@everywhere us_aut_G_l3Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G)

# l4G loc
# l4G => l1G
@everywhere cc_aut_G_l4Gl1G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l4G => l2G
@everywhere cc_aut_G_l4Gl2G_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t2])
@everywhere us_aut_G_l4Gl2G_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2G; 
 S[:d] +=  S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); 
 S[:tprime] = 0.0)

# Connection between the two automata: l2G => l1F
@everywhere cc_aut_F_l2Gl1F_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l1F loc : we construct  the edges of the form l1F => (..)
# l1F => l2F
@everywhere cc_aut_F_l1Fl2F_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) >= constants[:t3] &&
(constants[:x3] <= S[:n] <= constants[:x4])
@everywhere us_aut_F_l1Fl2F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2F;
 S[:dprime] = 0)

@everywhere cc_aut_F_l1Fl2F_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
getfield(S, :time) >= constants[:t3] &&
S[:dprime] == 0 
@everywhere us_aut_F_l1Fl2F_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2F)

@everywhere cc_aut_F_l1Fl2F_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t4]) && 
(S[:n] < constants[:x3] || S[:n] > constants[:x4])
@everywhere us_aut_F_l1Fl2F_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2F;
 S[:dprime] = min(abs(S[:n] - constants[:x3]), abs(S[:n] - constants[:x4]));
 S[:d] += S[:dprime])

@everywhere cc_aut_F_l1Fl2F_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
istrue(S[:isabs]) && getfield(S, :time) <= constants[:t4]
@everywhere us_aut_F_l1Fl2F_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2F;
 S[:d] += S[:dprime])

# l1F => l3F
@everywhere cc_aut_F_l1Fl3F_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(constants[:x3] <= S[:n] <= constants[:x4])
@everywhere us_aut_F_l1Fl3F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3F;
 S[:dprime] = 0;)

@everywhere cc_aut_F_l1Fl3F_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S[:n] < constants[:x3] || S[:n] > constants[:x4]) && 
(getfield(S, :time) <= constants[:t3])
@everywhere us_aut_F_l1Fl3F_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3F;
 S[:dprime] = min(sqrt((getfield(S, :time) - constants[:t3])^2 + (S[:n] - constants[:x4])^2), 
                  sqrt((getfield(S, :time) - constants[:t3])^2 + (S[:n] - constants[:x3])^2)))

@everywhere cc_aut_F_l1Fl3F_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S[:n] < constants[:x3] || S[:n] > constants[:x4]) && 
(constants[:t3] <= getfield(S, :time) <= constants[:t4])
@everywhere us_aut_F_l1Fl3F_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l3F;
 S[:dprime] = min(S[:dprime], min(abs(S[:n] - constants[:x3]), abs(S[:n] - constants[:x4]))))

# l3F loc
# l3F => l1F
@everywhere cc_aut_F_l3Fl1F_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = true

# l3F => l2F
@everywhere cc_aut_F_l3Fl2F_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(getfield(S, :time) >= constants[:t4])
@everywhere us_aut_F_l3Fl2F_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(S.loc = :l2F)

function create_automaton_G_and_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs_G::VariableModel,
                                  x3::Float64, x4::Float64, t3::Float64, t4::Float64, sym_obs_F::VariableModel)
    # Requirements for the automaton
    @assert sym_obs_G in m.g && sym_obs_F in m.g "$(sym_obs_G) or $(sym_obs_F) are not observed."
    @assert (x1 <= x2) "x1 > x2 impossible for G and F automaton."
    @assert (t1 <= t2) "t1 > t2 impossible for G and F automaton."
    @assert (x3 <= x4) "x3 > x3 impossible for G and F automaton."
    @assert (t3 <= t4) "t3 > t4 impossible for G and F automaton."
    @assert (t2 <= t3) "t3 > t2 impossible for G and F automaton."
    
    # Locations
    locations = [:l0G, :l1G, :l2G, :l3G, :l4G,
                 :l1F, :l2F, :l3F]

    # Invariant predicates
    Λ_F = Dict(:l0G => getfield(Main, :true_inv_predicate), :l1G => getfield(Main, :true_inv_predicate),
               :l2G => getfield(Main, :true_inv_predicate), :l3G => getfield(Main, :true_inv_predicate), 
               :l4G => getfield(Main, :true_inv_predicate),
               :l1F => getfield(Main, :true_inv_predicate),
               :l2F => getfield(Main, :true_inv_predicate), :l3F => getfield(Main, :true_inv_predicate))

    ## Init and final loc
    locations_init = [:l0G]
    locations_final = [:l2F]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:tprime => 1, :in => 2,
                                                        :n => 3,  :d => 4, 
                                                        :dprime => 5, :isabs => 6)

    ## Flow of variables
    flow = Dict{Location,Vector{Float64}}(:l0G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l1G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l2G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l3G => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l4G => [1.0,0.0,0.0,0.0,0.0,0.0],
                                          :l1F => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l2F => [0.0,0.0,0.0,0.0,0.0,0.0], 
                                          :l3F => [0.0,0.0,0.0,0.0,0.0,0.0])

    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end

    sym_isabs_func = Symbol(m.isabsorbing)
    idx_obs_var_F = getfield(m, :map_var_idx)[sym_obs_F]
    idx_obs_var_G = getfield(m, :map_var_idx)[sym_obs_G]
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    
    # l0G loc
    # l0G => l1G
    sym_func_us_l0Gl1G_1 = Symbol("us_aut_G_$(basename_func)_l0Gl1G_1!")
    str_us_l0Gl1G_1 = "
    @everywhere $(sym_func_us_l0Gl1G_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1G; \n
     S[:d] = 0; \n
     S[:n] = x[$(idx_obs_var_G)]; \n
     S[:in] = true; \n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l0Gl1G_1))
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l0Gl1G_1), getfield(Main, sym_func_us_l0Gl1G_1))
    map_edges[:l0G][:l1G] = [edge1]

    # l1G => l3G
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl3G_1), getfield(Main, :us_aut_G_l1Gl3G_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl3G_2), getfield(Main, :us_aut_G_l1Gl3G_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl3G_3), getfield(Main, :us_aut_G_l1Gl3G_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl3G_4), getfield(Main, :us_aut_G_l1Gl3G_4!))
    map_edges[:l1G][:l3G] = [edge1, edge2, edge3, edge4]
    
    # l1G => l4G
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl4G_1), getfield(Main, :us_aut_G_l1Gl4G_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl4G_2), getfield(Main, :us_aut_G_l1Gl4G_2!))
    map_edges[:l1G][:l4G] = [edge1, edge2]
    
    # l1G => l2G
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl2G_1), getfield(Main, :us_aut_G_l1Gl2G_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl2G_2), getfield(Main, :us_aut_G_l1Gl2G_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl2G_3), getfield(Main, :us_aut_G_l1Gl2G_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_G_l1Gl2G_4), getfield(Main, :us_aut_G_l1Gl2G_4!))
    map_edges[:l1G][:l2G] = [edge3, edge4, edge1, edge2]

    # l3G loc
    # l3G => l1G
    sym_func_us_l3Gl1G_1 = Symbol("us_aut_G_$(basename_func)_l3Gl1G_1!")
    str_us_l3Gl1G_1 = "
    @everywhere $(sym_func_us_l3Gl1G_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1G; \n
     S[:n] = x[$(idx_obs_var_G)]; \n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l3Gl1G_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_G_l3Gl1G_1), getfield(Main, sym_func_us_l3Gl1G_1))
    map_edges[:l3G][:l1G] = [edge1]
    
    # l3G => l2G
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l3Gl2G_1), getfield(Main, :us_aut_G_l3Gl2G_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_G_l3Gl2G_2), getfield(Main, :us_aut_G_l3Gl2G_2!))
    map_edges[:l3G][:l2G] = [edge1, edge2]

    # l4 loc
    # l4G => l1G
    sym_func_us_l4Gl1G_1 = Symbol("us_aut_G_$(basename_func)_l4Gl1G_1!")
    str_us_l4Gl1G_1 = "
    @everywhere $(sym_func_us_l4Gl1G_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1G; \n
     S[:d] += S[:tprime] * min(abs(constants[:x1] - S[:n]), abs(constants[:x2] - S[:n])); \n
     S[:tprime] = 0.0; \n
     S[:n] = x[$(idx_obs_var_G)]; \n
     S[:in] = true; \n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l4Gl1G_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_G_l4Gl1G_1), getfield(Main, sym_func_us_l4Gl1G_1))
    map_edges[:l4G][:l1G] = [edge1]

    # l4G => l2G
    edge1 = Edge([nothing], getfield(Main, :cc_aut_G_l4Gl2G_1), getfield(Main, :us_aut_G_l4Gl2G_1!))
    map_edges[:l4G][:l2G] = [edge1]

    # l2G loc
    # l2G => l1F : Transition from autF to autG
    sym_func_us_l2Gl1F_1 = Symbol("us_aut_G_$(basename_func)_l2Gl1F_1!")
    str_us_l2Gl1F_1 = "
    @everywhere $(sym_func_us_l2Gl1F_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1F; \n
     S[:n] = x[$(idx_obs_var_F)];\n
     S[:dprime] = Inf; \n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l2Gl1F_1))
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l2Gl1F_1), getfield(Main, sym_func_us_l2Gl1F_1))
    map_edges[:l2G][:l1F] = [edge1]

    # l1F loc
    # l1F => l3F
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl2F_1), getfield(Main, :us_aut_F_l1Fl2F_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl2F_2), getfield(Main, :us_aut_F_l1Fl2F_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl2F_3), getfield(Main, :us_aut_F_l1Fl2F_3!))
    edge4 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl2F_4), getfield(Main, :us_aut_F_l1Fl2F_4!))
    map_edges[:l1F][:l2F] = [edge1, edge4, edge3, edge2]

    # l1F => l3F
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl3F_1), getfield(Main, :us_aut_F_l1Fl3F_1!))
    edge2 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl3F_2), getfield(Main, :us_aut_F_l1Fl3F_2!))
    edge3 = Edge([nothing], getfield(Main, :cc_aut_F_l1Fl3F_3), getfield(Main, :us_aut_F_l1Fl3F_3!))
    map_edges[:l1F][:l3F] = [edge1, edge2, edge3]

    # l3F loc
    # l3F => l1F
    sym_func_us_l3Fl1F_1 = Symbol("us_aut_G_$(basename_func)_l3Fl1F_1!")
    str_us_l3Fl1F_1 = "
    @everywhere $(sym_func_us_l3Fl1F_1)(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = \n
    (S.loc = :l1F;\n
     S[:n] = x[$(idx_obs_var_F)];\n
     S[:isabs] = getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x))"
    eval(Meta.parse(str_us_l3Fl1F_1))
    edge1 = Edge([:ALL], getfield(Main, :cc_aut_F_l3Fl1F_1), getfield(Main, sym_func_us_l3Fl1F_1))
    map_edges[:l3F][:l1F] = [edge1]

    # l3F => l2F
    edge1 = Edge([nothing], getfield(Main, :cc_aut_F_l3Fl2F_1), getfield(Main, :us_aut_F_l3Fl2F_1!))    
    map_edges[:l3F][:l2F] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2,
                                     :x3 => x3,  :x4 => x4, :t3 => t3, :t4 => t4)

    A = LHA("G and F property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_G_and_F

