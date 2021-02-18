
function create_automaton_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert (x1 <= x2) "x1 > x2 impossible for F automaton."
    @assert (t1 <= t2) "t1 > t2 impossible for F automaton."

    # Locations
    locations = [:l0, :l1, :l2, :l3]

    ## Invariant predicates
    @everywhere true_inv_predicate(x::Vector{Int}) = true 
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

    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    idx_var_n = map_var_automaton_idx[:n] 
    idx_var_d = map_var_automaton_idx[:d] 
    idx_var_isabs = map_var_automaton_idx[:isabs] 

    nbr_rand = rand(1:100000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    sym_isabs_func = Symbol(m.isabsorbing)
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_aut_F_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote 
        @everywhere istrue(val::Float64) = convert(Bool, val)

        ## Check constraints and update state functions
        # l0 loc : we construct  the edges of the form l0 => (..)
        # "cc" as check_constraints and "us" as update_state
        # l0 => l1
        @everywhere $(func_name(:cc, :l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(idx_var_n)] = x[$(idx_obs_var)];
         setindex!(S_values, Inf, $(idx_var_d)); 
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs));
         :l1)
        # l1 loc
        # l1 => l2
        #=
        @everywhere $(func_name(:cc, :l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t1 &&
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (setindex!(S_values, 0, $(idx_var_d));
        :l2)        =#
        #=
        @everywhere $(func_name(:cc, :l1, :l2, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && S_time <= $t2
        @everywhere $(func_name(:us, :l1, :l2, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2)
        =#
        @everywhere $(func_name(:cc, :l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t1 &&
        S_values[$(idx_var_d)] == 0 
        @everywhere $(func_name(:us, :l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2)

        @everywhere $(func_name(:cc, :l1, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2) && 
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2;)
        #setindex!(S_values, min(abs(S_values[$(idx_var_n)] - $x1), abs(S_values[$(idx_var_n)] - $x2)), $(idx_var_d)))

        # l1 => l3
        @everywhere $(func_name(:cc, :l1, :l3, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t1) &&
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1, :l3, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (setindex!(S_values, min(sqrt((S_time - $t1)^2 + (S_values[$(idx_var_n)] - $x2)^2), 
                                 sqrt((S_time - $t1)^2 + (S_values[$(idx_var_n)] - $x1)^2)), $(idx_var_d));
        :l3)        
        @everywhere $(func_name(:cc, :l1, :l3, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1, :l3, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (setindex!(S_values, 0, $(idx_var_d));
         :l3)
        @everywhere $(func_name(:cc, :l1, :l3, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t1) &&
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1, :l3, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (val_min = min(S_values[$(idx_var_d)], 
                       min(abs(S_values[$(idx_var_n)] - $x1), abs(S_values[$(idx_var_n)] - $x2)));
         setindex!(S_values, val_min, $(idx_var_d));
        :l3)
        # l3 loc
        # l3 => l1
        @everywhere $(func_name(:cc, :l3, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l3, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(idx_var_n)] = x[$(idx_obs_var)];
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs));
         :l1)
        # l3 => l2
        @everywhere $(func_name(:cc, :l3, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2)
    end
    eval(meta_elementary_functions)

    # l0 loc
    # l0 => l1
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :l1, 1)), getfield(Main, func_name(:us, :l0, :l1, 1)))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l2
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 1)), getfield(Main, func_name(:us, :l1, :l2, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 2)), getfield(Main, func_name(:us, :l1, :l2, 2)))
    map_edges[:l1][:l2] = [edge1, edge2]
    #edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 3)), getfield(Main, func_name(:us, :l1, :l2, 3)))
    #edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 4)), getfield(Main, func_name(:us, :l1, :l2, 4)))
    #map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]

    # l1 => l3
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 1)), getfield(Main, func_name(:us, :l1, :l3, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 2)), getfield(Main, func_name(:us, :l1, :l3, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 3)), getfield(Main, func_name(:us, :l1, :l3, 3)))
    map_edges[:l1][:l3] = [edge1, edge2, edge3]

    # l3 loc
    # l3 => l1
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l3, :l1, 1)), getfield(Main, func_name(:us, :l3, :l1, 1)))
    map_edges[:l3][:l1] = [edge1]

    # l3 => l2
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l3, :l2, 1)), getfield(Main, func_name(:us, :l3, :l2, 1)))
    map_edges[:l3][:l2] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA("F property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_F

