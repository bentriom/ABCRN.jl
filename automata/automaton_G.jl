

function create_automaton_G(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert (x1 <= x2) "x1 > x2 impossible for G automaton."
    @assert (t1 <= t2) "t1 > t2 impossible for G automaton."
    
    # Locations
    locations = [:l0, :l1, :l2, :l3, :l4]

    # Invariant predicates
    @everywhere true_inv_predicate(x::Vector{Int}) = true 
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

    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    idx_var_n = map_var_automaton_idx[:n] 
    idx_var_d = map_var_automaton_idx[:d] 
    idx_var_isabs = map_var_automaton_idx[:isabs] 
    idx_var_in = map_var_automaton_idx[:in] 
    idx_var_tprime = map_var_automaton_idx[:tprime]

    nbr_rand = rand(1:100000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    sym_isabs_func = Symbol(m.isabsorbing)
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_aut_G_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote
        @everywhere istrue(val::Float64) = convert(Bool, val)
        # l0 loc
        # l0 => l1
        @everywhere $(func_name(:cc, :l0, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l0, :l1, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1"); 
         setindex!(S_values, 0, $(idx_var_d));
         setindex!(S_values, x[$(idx_obs_var)], $(idx_var_n));
         setindex!(S_values, true, $(idx_var_in));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l1 => l3
        @everywhere $(func_name(:cc, :l1, :l3, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time <= $t1 && 
        S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2
        @everywhere $(func_name(:us, :l1, :l3, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3"); 
         setindex!(S_values, min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));
         setindex!(S_values, false, $(idx_var_in)))

        @everywhere $(func_name(:cc, :l1, :l3, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t1) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1, :l3, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3"); 
         setindex!(S_values, 0, $(idx_var_d));
         setindex!(S_values, false, $(idx_var_in)))

        @everywhere $(func_name(:cc, :l1, :l3, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1, :l3, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3"); 
         setindex!(S_values, S_values[$(idx_var_d)] * (S_time - $t1), $(idx_var_d));
         setindex!(S_values, 0.0, $(idx_var_tprime)))

        @everywhere $(func_name(:cc, :l1, :l3, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1, :l3, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3"); 
         setindex!(S_values, 0.0, $(idx_var_tprime)))

        # l1 => l4
        @everywhere $(func_name(:cc, :l1, :l4, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1, :l4, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l4"); 
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_d)] * (S_time - $t1), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1, :l4, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1, :l4, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l4"))


        # l1 => l2
        #=
        @everywhere $(func_name(:cc, :l1, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        S_time >= $t2
        @everywhere $(func_name(:us, :l1, :l2, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"))

        @everywhere $(func_name(:cc, :l1, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        S_time >= $t2
        @everywhere $(func_name(:us, :l1, :l2, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"); 
         setindex!(S_values, S_values[$(idx_var_d)] * ($t2 - $t1), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1, :l2, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && 
        S_time <= $t1
        @everywhere $(func_name(:us, :l1, :l2, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"); 
         setindex!(S_values, ($t2 - $t1) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1, :l2, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && 
        ($t1 <= S_time <= $t2)
        @everywhere $(func_name(:us, :l1, :l2, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"); 
         setindex!(S_values, S_values[$(idx_var_d)] + ($t2 - S_time) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)))
        =#

        # l3 => l1
        @everywhere $(func_name(:cc, :l3, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l3, :l1, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1"); 
         setindex!(S_values, x[$(idx_obs_var)], $(idx_var_n));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l4 => l1
        @everywhere $(func_name(:cc, :l4, :l1, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l4, :l1, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1"); 
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_tprime)] * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));
         setindex!(S_values, 0.0, $(idx_var_tprime));
         setindex!(S_values, x[$(idx_obs_var)], $(idx_var_n));
         setindex!(S_values, true, $(idx_var_in));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l3 => l2
        @everywhere $(func_name(:cc, :l3, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        (S_time >= $t2 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3, :l2, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2");
         setindex!(S_values, S_values[$(idx_var_d)] * ($t2 - $t1), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l3, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        (S_time >= $t2 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3, :l2, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"))

        # l4 => l2
        @everywhere $(func_name(:cc, :l4, :l2, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l4, :l2, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"); 
         setindex!(S_values, S_values[$(idx_var_d)] + ($t2 - S_time) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));
         setindex!(S_values, 0.0, $(idx_var_tprime)))

        @everywhere $(func_name(:cc, :l4, :l2, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2)
        @everywhere $(func_name(:us, :l4, :l2, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2"); 
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_tprime)] * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));
         setindex!(S_values, 0.0, $(idx_var_tprime)))
 
    end
    eval(meta_elementary_functions)

    # l0 loc
    # l0 => l1
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :l1, 1)), getfield(Main, func_name(:us, :l0, :l1, 1)))
    map_edges[:l0][:l1] = [edge1]

    # l1 loc
    # l1 => l3
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 1)), getfield(Main, func_name(:us, :l1, :l3, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 2)), getfield(Main, func_name(:us, :l1, :l3, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 3)), getfield(Main, func_name(:us, :l1, :l3, 3)))
    edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l3, 4)), getfield(Main, func_name(:us, :l1, :l3, 4)))
    map_edges[:l1][:l3] = [edge1, edge2, edge3, edge4]

    # l1 => l4
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l4, 1)), getfield(Main, func_name(:us, :l1, :l4, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l4, 2)), getfield(Main, func_name(:us, :l1, :l4, 2)))
    map_edges[:l1][:l4] = [edge1, edge2]
   
    # l1 => l2
    #=
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 1)), getfield(Main, func_name(:us, :l1, :l2, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 2)), getfield(Main, func_name(:us, :l1, :l2, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 3)), getfield(Main, func_name(:us, :l1, :l2, 3)))
    edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1, :l2, 4)), getfield(Main, func_name(:us, :l1, :l2, 4)))
    map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]
    =#

    # l3 loc
    # l3 => l1
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l3, :l1, 1)), getfield(Main, func_name(:us, :l3, :l1, 1)))
    map_edges[:l3][:l1] = [edge1]

    # l3 => l2
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l3, :l2, 1)), getfield(Main, func_name(:us, :l3, :l2, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l3, :l2, 2)), getfield(Main, func_name(:us, :l3, :l2, 2)))
    map_edges[:l3][:l2] = [edge1, edge2]

    # l4 loc
    # l4 => l1
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l4, :l1, 1)), getfield(Main, func_name(:us, :l4, :l1, 1)))
    map_edges[:l4][:l1] = [edge1]

    # l4 => l2
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l4, :l2, 1)), getfield(Main, func_name(:us, :l4, :l2, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l4, :l2, 2)), getfield(Main, func_name(:us, :l4, :l2, 2)))
    map_edges[:l4][:l2] = [edge1,edge2]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    A = LHA("G property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
   
end

export create_automaton_G

