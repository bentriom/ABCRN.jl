
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
    @everywhere true_inv_predicate(x::Vector{Int}) = true 
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

    idx_obs_var_F = getfield(m, :map_var_idx)[sym_obs_F]
    idx_obs_var_G = getfield(m, :map_var_idx)[sym_obs_G]
    idx_var_n = map_var_automaton_idx[:n] 
    idx_var_d = map_var_automaton_idx[:d] 
    idx_var_dprime = map_var_automaton_idx[:dprime] 
    idx_var_isabs = map_var_automaton_idx[:isabs] 
    idx_var_in = map_var_automaton_idx[:in] 
    idx_var_tprime = map_var_automaton_idx[:tprime]
    
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    sym_isabs_func = Symbol(m.isabsorbing)
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_aut_F_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote
        @everywhere istrue(val::Float64) = convert(Bool, val)
        ## Edges check constraint and update state functions

        # l0G loc
        # l0G => l1G
        @everywhere $(func_name(:cc, :l0G, :l1G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l0G, :l1G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1G"); 
         setindex!(S_values, 0, $(idx_var_d));
         setindex!(S_values, x[$(idx_obs_var_G)], $(idx_var_n)); 
         setindex!(S_values, true, $(idx_var_in));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l1G loc
        # l1G => l3G
        @everywhere $(func_name(:cc, :l1G, :l3G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time <= $t1 && 
        S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2
        @everywhere $(func_name(:us, :l1G, :l3G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3G"); 
         setindex!(S_values, min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)); 
         setindex!(S_values, false, $(idx_var_in)))

        @everywhere $(func_name(:cc, :l1G, :l3G, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t1) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1G, :l3G, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3G"); 
         setindex!(S_values, 0, $(idx_var_d));
         setindex!(S_values, false, $(idx_var_in)))

        @everywhere $(func_name(:cc, :l1G, :l3G, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1G, :l3G, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3G"); 
         setindex!(S_values, S_values[$(idx_var_d)] * (S_time - $t1), $(idx_var_d)); 
         setindex!(S_values, 0.0, $(idx_var_tprime)))

        @everywhere $(func_name(:cc, :l1G, :l3G, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(idx_var_n)] <= $x2)
        @everywhere $(func_name(:us, :l1G, :l3G, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3G"); 
         setindex!(S_values, 0.0, $(idx_var_tprime)))

        # l1G => l4G
        @everywhere $(func_name(:cc, :l1G, :l4G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1G, :l4G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l4G"); 
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_d)] * (S_time - $t1), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1G, :l4G, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(idx_var_n)] < $x1 || S_values[$(idx_var_n)] > $x2)
        @everywhere $(func_name(:us, :l1G, :l4G, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l4G"))

        # l1G => l2G
        #=
        @everywhere $(func_name(:cc, :l1G, :l2G, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && 
        S_time <= $t1
        @everywhere $(func_name(:us, :l1G, :l2G, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"); 
         setindex!(S_values, ($t2 - $t1) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)))
        
         @everywhere $(func_name(:cc, :l1G, :l2G, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && 
        ($t1 <= S_time <= $t2)
        @everywhere $(func_name(:us, :l1G, :l2G, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"); 
         setindex!(S_values, S_values[$(idx_var_d)] + ($t2 - S_time) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1G, :l2G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        S_time >= $t2
        @everywhere $(func_name(:us, :l1G, :l2G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"))

        @everywhere $(func_name(:cc, :l1G, :l2G, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        S_time >= $t2
        @everywhere $(func_name(:us, :l1G, :l2G, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"); 
         setindex!(S_values, S_values[$(idx_var_d)] * ($t2 - $t1), $(idx_var_d)))
        =#

        # l3G loc
        # l3G => l1G
        @everywhere $(func_name(:cc, :l3G, :l1G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l3G, :l1G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1G"); 
         setindex!(S_values, x[$(idx_obs_var_G)], $(idx_var_n)); 
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l3G => l2G
        @everywhere $(func_name(:cc, :l3G, :l2G, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_in)]) && 
        (S_time >= $t2 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3G, :l2G, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G");
         setindex!(S_values, S_values[$(idx_var_d)] * ($t2 - $t1), $(idx_var_d)))

        @everywhere $(func_name(:cc, :l3G, :l2G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(idx_var_in)]) && 
        (S_time >= $t2 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3G, :l2G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"))

        # l4G loc
        # l4G => l1G
        @everywhere $(func_name(:cc, :l4G, :l1G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l4G, :l1G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1G"); 
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_tprime)] * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d)); 
         setindex!(S_values, 0.0, $(idx_var_tprime));
         setindex!(S_values, x[$(idx_obs_var_G)], $(idx_var_n)); 
         setindex!(S_values, true, $(idx_var_in));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l4G => l2G
        @everywhere $(func_name(:cc, :l4G, :l2G, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l4G, :l2G, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"); 
         setindex!(S_values, S_values[$(idx_var_d)] +  ($t2 - S_time) * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));)
        
        @everywhere $(func_name(:cc, :l4G, :l2G, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2)
        @everywhere $(func_name(:us, :l4G, :l2G, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2G"); 
         setindex!(S_values, S_values[$(idx_var_d)] +  S_values[$(idx_var_tprime)] * min(abs($x1 - S_values[$(idx_var_n)]), abs($x2 - S_values[$(idx_var_n)])), $(idx_var_d));)


        # Connection between the two automata: l2G => l1F
        @everywhere $(func_name(:cc, :l2G, :l1F, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l2G, :l1F, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1F"); 
         setindex!(S_values, x[$(idx_obs_var_F)], $(idx_var_n));
         setindex!(S_values, Inf, $(idx_var_dprime));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l1F loc : we construct  the edges of the form l1F => (..)
        # l1F => l2F
        @everywhere $(func_name(:cc, :l1F, :l2F, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t3 &&
        S_values[$(idx_var_dprime)] == 0 
        @everywhere $(func_name(:us, :l1F, :l2F, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2F"));
         #setindex!(S_values, 0, $(idx_var_dprime)))

        @everywhere $(func_name(:cc, :l1F, :l2F, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t4) && 
        (S_values[$(idx_var_n)] < $x3 || S_values[$(idx_var_n)] > $x4)
        @everywhere $(func_name(:us, :l1F, :l2F, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2F");
         #setindex!(S_values, min(abs(S_values[$(idx_var_n)] - $x3), abs(S_values[$(idx_var_n)] - $x4)), $(idx_var_dprime));
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_dprime)], $(idx_var_d)))
        #=
        @everywhere $(func_name(:cc, :l1F, :l2F, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(idx_var_isabs)]) && S_time <= $t4
        @everywhere $(func_name(:us, :l1F, :l2F, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2F");
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_dprime)], $(idx_var_d)))

        @everywhere $(func_name(:cc, :l1F, :l2F, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t3 &&
        S_values[$(idx_var_dprime)] == 0 
        @everywhere $(func_name(:us, :l1F, :l2F, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2F"))
        =#

        # l1F => l3F
        @everywhere $(func_name(:cc, :l1F, :l3F, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t3) &&
        (S_values[$(idx_var_n)] < $x3 || S_values[$(idx_var_n)] > $x4)
        @everywhere $(func_name(:us, :l1F, :l3F, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3F");
         setindex!(S_values, min(sqrt((S_time - $t3)^2 + (S_values[$(idx_var_n)] - $x4)^2), 
                                             sqrt((S_time - $t3)^2 + (S_values[$(idx_var_n)] - $x3)^2)), $(idx_var_dprime)))

        @everywhere $(func_name(:cc, :l1F, :l3F, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        ($x3 <= S_values[$(idx_var_n)] <= $x4)
        @everywhere $(func_name(:us, :l1F, :l3F, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3F");
         setindex!(S_values, 0, $(idx_var_dprime));)

        @everywhere $(func_name(:cc, :l1F, :l3F, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t3) &&
        (S_values[$(idx_var_n)] < $x3 || S_values[$(idx_var_n)] > $x4)
        @everywhere $(func_name(:us, :l1F, :l3F, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l3F");
         setindex!(S_values, min(S_values[$(idx_var_dprime)], min(abs(S_values[$(idx_var_n)] - $x3), abs(S_values[$(idx_var_n)] - $x4))), $(idx_var_dprime)))

        # l3F loc
        # l3F => l1F
        @everywhere $(func_name(:cc, :l3F, :l1F, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        @everywhere $(func_name(:us, :l3F, :l1F, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l1F");
         setindex!(S_values, x[$(idx_obs_var_F)], $(idx_var_n));
         setindex!(S_values, getfield(Main, $(Meta.quot(sym_isabs_func)))(p, x), $(idx_var_isabs)))

        # l3F => l2F
        @everywhere $(func_name(:cc, :l3F, :l2F, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t4 || istrue(S_values[$(idx_var_isabs)]))
        @everywhere $(func_name(:us, :l3F, :l2F, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (ptr_loc[1] = Symbol("l2F");
         setindex!(S_values, S_values[$(idx_var_d)] + S_values[$(idx_var_dprime)], $(idx_var_d)))
    end
    eval(meta_elementary_functions)

    # l0G loc
    # l0G => l1G
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l0G, :l1G, 1)), getfield(Main, func_name(:us, :l0G, :l1G, 1)))
    map_edges[:l0G][:l1G] = [edge1]

    # l1G => l3G
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l3G, 1)), getfield(Main, func_name(:us, :l1G, :l3G, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l3G, 2)), getfield(Main, func_name(:us, :l1G, :l3G, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l3G, 3)), getfield(Main, func_name(:us, :l1G, :l3G, 3)))
    edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l3G, 4)), getfield(Main, func_name(:us, :l1G, :l3G, 4)))
    map_edges[:l1G][:l3G] = [edge1, edge2, edge3, edge4]

    # l1G => l4G
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l4G, 1)), getfield(Main, func_name(:us, :l1G, :l4G, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l4G, 2)), getfield(Main, func_name(:us, :l1G, :l4G, 2)))
    map_edges[:l1G][:l4G] = [edge1, edge2]

    # l1G => l2G
    #=
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l2G, 1)), getfield(Main, func_name(:us, :l1G, :l2G, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l2G, 2)), getfield(Main, func_name(:us, :l1G, :l2G, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l2G, 3)), getfield(Main, func_name(:us, :l1G, :l2G, 3)))
    edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1G, :l2G, 4)), getfield(Main, func_name(:us, :l1G, :l2G, 4)))
    map_edges[:l1G][:l2G] = [edge3, edge4, edge1, edge2]
    =#

    # l3G loc
    # l3G => l1G
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l3G, :l1G, 1)), getfield(Main, func_name(:us, :l3G, :l1G, 1)))
    map_edges[:l3G][:l1G] = [edge1]

    # l3G => l2G
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l3G, :l2G, 1)), getfield(Main, func_name(:us, :l3G, :l2G, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l3G, :l2G, 2)), getfield(Main, func_name(:us, :l3G, :l2G, 2)))
    map_edges[:l3G][:l2G] = [edge1, edge2]

    # l4 loc
    # l4G => l1G
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l4G, :l1G, 1)), getfield(Main, func_name(:us, :l4G, :l1G, 1)))
    map_edges[:l4G][:l1G] = [edge1]

    # l4G => l2G
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l4G, :l2G, 1)), getfield(Main, func_name(:us, :l4G, :l2G, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l4G, :l2G, 2)), getfield(Main, func_name(:us, :l4G, :l2G, 2)))
    map_edges[:l4G][:l2G] = [edge1,edge2]

    # l2G loc
    # l2G => l1F : Transition from autF to autG
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l2G, :l1F, 1)), getfield(Main, func_name(:us, :l2G, :l1F, 1)))
    map_edges[:l2G][:l1F] = [edge1]

    # l1F loc
    # l1F => l3F
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l2F, 1)), getfield(Main, func_name(:us, :l1F, :l2F, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l2F, 2)), getfield(Main, func_name(:us, :l1F, :l2F, 2)))
    map_edges[:l1F][:l2F] = [edge1, edge2]
    #edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l2F, 3)), getfield(Main, func_name(:us, :l1F, :l2F, 3)))
    #edge4 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l2F, 4)), getfield(Main, func_name(:us, :l1F, :l2F, 4)))
    #map_edges[:l1F][:l2F] = [edge1, edge4, edge3, edge2]

    # l1F => l3F
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l3F, 1)), getfield(Main, func_name(:us, :l1F, :l3F, 1)))
    edge2 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l3F, 2)), getfield(Main, func_name(:us, :l1F, :l3F, 2)))
    edge3 = Edge(nothing, getfield(Main, func_name(:cc, :l1F, :l3F, 3)), getfield(Main, func_name(:us, :l1F, :l3F, 3)))
    map_edges[:l1F][:l3F] = [edge1, edge2, edge3]

    # l3F loc
    # l3F => l1F
    edge1 = Edge([:ALL], getfield(Main, func_name(:cc, :l3F, :l1F, 1)), getfield(Main, func_name(:us, :l3F, :l1F, 1)))
    map_edges[:l3F][:l1F] = [edge1]

    # l3F => l2F
    edge1 = Edge(nothing, getfield(Main, func_name(:cc, :l3F, :l2F, 1)), getfield(Main, func_name(:us, :l3F, :l2F, 1)))    
    map_edges[:l3F][:l2F] = [edge1]

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2,
                                     :x3 => x3,  :x4 => x4, :t3 => t3, :t4 => t4)

    A = LHA("G and F property", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_automaton_G_and_F

