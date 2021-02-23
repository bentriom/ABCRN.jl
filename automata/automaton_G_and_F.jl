
# Creation of the automaton types
#@everywhere @eval abstract type EdgeAutomatonGandF <: Edge end
@everywhere struct EdgeAutomatonGandF{T} <: Edge transitions::Union{Nothing,Vector{Symbol}} end
@everywhere @eval $(MarkovProcesses.generate_code_lha_type_def(:AutomatonGandF, :EdgeAutomatonGandF))

function create_automaton_G_and_F(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs_G::VariableModel,
                                  x3::Float64, x4::Float64, t3::Float64, t4::Float64, sym_obs_F::VariableModel)
    # Requirements for the automaton
    @assert sym_obs_G in m.g && sym_obs_F in m.g "$(sym_obs_G) or $(sym_obs_F) are not observed."
    @assert (x1 <= x2) "x1 > x2 impossible for G and F automaton."
    @assert (t1 <= t2) "t1 > t2 impossible for G and F automaton."
    @assert (x3 <= x4) "x3 > x3 impossible for G and F automaton."
    @assert (t3 <= t4) "t3 > t4 impossible for G and F automaton."
    @assert (t2 <= t3) "t3 > t2 impossible for G and F automaton."

    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :AutomatonGandF
    edge_type = :EdgeAutomatonGandF
    check_constraints = Symbol("check_constraints_$(lha_name)")
    update_state! = Symbol("update_state_$(lha_name)!")

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
    idx_obs_var_F = getfield(m, :map_var_idx)[sym_obs_F]
    idx_obs_var_G = getfield(m, :map_var_idx)[sym_obs_G]
    to_idx(var::Symbol) = map_var_automaton_idx[var]

    id = MarkovProcesses.newid()
    #Symbol("Edge_$(lha_name)_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)")
    function generate_edge_type(from_loc::Location, to_loc::Location, edge_number::Int)
        tag_type = Symbol("$(from_loc)$(to_loc)_$(edge_number)_$(model_name)_$(id)")
        @eval type = EdgeAutomatonGandF{$(Meta.quot(tag_type))}
        return @eval($type)
    end

    ## check_constraints & update_state!
    @everywhere @eval begin
    #return quote
        istrue(val::Float64) = convert(Bool, val)
        ## Edges check constraint and update state functions

        # l0G loc
        # l0G => l1G
        #struct $(generate_edge_type(:l0G, :l1G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l0G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(generate_edge_type(:l0G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = 0;
         S_values[$(to_idx(:n))] = x[$(idx_obs_var_G)];
         S_values[$(to_idx(:in))] = true;
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1G)

        # l1G loc
        # l1G => l3G
        #struct $(generate_edge_type(:l1G, :l3G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l3G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time <= $t1 && 
        S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l3G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:in))] = false;
         :l3G)

        #struct $(generate_edge_type(:l1G, :l3G, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l3G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t1) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l3G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = 0;
         S_values[$(to_idx(:in))] = false;
         :l3G)

        #struct $(generate_edge_type(:l1G, :l3G, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l3G, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l3G, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * (S_time - $t1);
         S_values[$(to_idx(:tprime))] = 0.0;
         :l3G)

        #struct $(generate_edge_type(:l1G, :l3G, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l3G, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l3G, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:tprime))] = 0.0;
         :l3G)

        # l1G => l4G
        #struct $(generate_edge_type(:l1G, :l4G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l4G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l4G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:d))] * (S_time - $t1);
         :l4G)

        #struct $(generate_edge_type(:l1G, :l4G, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l4G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l4G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l4G)

        # l1G => l2G
        #=
        #struct $(generate_edge_type(:l1G, :l2G, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l2G, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:isabs))]) && 
        S_time <= $t1
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l2G, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = ($t2 - $t1) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
        :l2G)

        #struct $(generate_edge_type(:l1G, :l2G, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l2G, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:isabs))]) && 
        ($t1 <= S_time <= $t2)
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l2G, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + ($t2 - S_time) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
        :l2G)

        #struct $(generate_edge_type(:l1G, :l2G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        S_time >= $t2
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2G)

        #struct $(generate_edge_type(:l1G, :l2G, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        S_time >= $t2
        $(update_state!)(edge::$(generate_edge_type(:l1G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * ($t2 - $t1);
        :l2G)
        =#

        # l3G loc
        # l3G => l1G
        #struct $(generate_edge_type(:l3G, :l1G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l3G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(generate_edge_type(:l3G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var_G)];
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1G)

        # l3G => l2G
        #struct $(generate_edge_type(:l3G, :l2G, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l3G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        (S_time >= $t2 || istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(generate_edge_type(:l3G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * ($t2 - $t1);
         :l2G)

        #struct $(generate_edge_type(:l3G, :l2G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l3G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        (S_time >= $t2 || istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(generate_edge_type(:l3G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2G)

        # l4G loc
        # l4G => l1G
        #struct $(generate_edge_type(:l4G, :l1G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l4G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(generate_edge_type(:l4G, :l1G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:tprime))] * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:tprime))] = 0.0;
         S_values[$(to_idx(:n))] = x[$(idx_obs_var_G)];
         S_values[$(to_idx(:in))] = true;
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1G)

        # l4G => l2G
        #struct $(generate_edge_type(:l4G, :l2G, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l4G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(generate_edge_type(:l4G, :l2G, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] +  ($t2 - S_time) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         :l2G)

        #struct $(generate_edge_type(:l4G, :l2G, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l4G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2)
        $(update_state!)(edge::$(generate_edge_type(:l4G, :l2G, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] +  S_values[$(to_idx(:tprime))] * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         :l2G)


        # Connection between the two automata: l2G => l1F
        #struct $(generate_edge_type(:l2G, :l1F, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l2G, :l1F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(generate_edge_type(:l2G, :l1F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var_F)];
         S_values[$(to_idx(:dprime))] = Inf;
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1F)

        # l1F loc : we con#struct  the edges of the form l1F => (..)
        # l1F => l2F
        #struct $(generate_edge_type(:l1F, :l2F, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l2F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t3 &&
        S_values[$(to_idx(:dprime))] == 0 
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l2F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (#S_values[$(to_idx(:dprime))] = 0;
         :l2F)

        #struct $(generate_edge_type(:l1F, :l2F, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l2F, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t4) && 
        (S_values[$(to_idx(:n))] < $x3 || S_values[$(to_idx(:n))] > $x4)
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l2F, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (#S_values[$(to_idx(:dprime))] = min(abs(S_values[$(to_idx(:n))] - $x3), abs(S_values[$(to_idx(:n))] - $x4));
         S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:dprime))];
         :l2F)
        #=
        #struct $(generate_edge_type(:l1F, :l2F, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l2F, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:isabs))]) && S_time <= $t4
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l2F, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:dprime))];
        :l2F)

        #struct $(generate_edge_type(:l1F, :l2F, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l2F, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time >= $t3 &&
        S_values[$(to_idx(:dprime))] == 0 
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l2F, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2F)
        =#

        # l1F => l3F
        #struct $(generate_edge_type(:l1F, :l3F, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l3F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t3) &&
        (S_values[$(to_idx(:n))] < $x3 || S_values[$(to_idx(:n))] > $x4)
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l3F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:dprime))] = min(sqrt((S_time - $t3)^2 + (S_values[$(to_idx(:n))] - $x4)^2), 
                                            sqrt((S_time - $t3)^2 + (S_values[$(to_idx(:n))] - $x3)^2));
         :l3F)

        #struct $(generate_edge_type(:l1F, :l3F, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l3F, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        ($x3 <= S_values[$(to_idx(:n))] <= $x4)
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l3F, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:dprime))] = 0;
         :l3F)

        #struct $(generate_edge_type(:l1F, :l3F, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l1F, :l3F, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t3) &&
        (S_values[$(to_idx(:n))] < $x3 || S_values[$(to_idx(:n))] > $x4)
        $(update_state!)(edge::$(generate_edge_type(:l1F, :l3F, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:dprime))] = min(S_values[$(to_idx(:dprime))], min(abs(S_values[$(to_idx(:n))] - $x3), abs(S_values[$(to_idx(:n))] - $x4)));
         :l3F)

        # l3F loc
        # l3F => l1F
        #struct $(generate_edge_type(:l3F, :l1F, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l3F, :l1F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(generate_edge_type(:l3F, :l1F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var_F)];
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1F)

        # l3F => l2F
        #struct $(generate_edge_type(:l3F, :l2F, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(generate_edge_type(:l3F, :l2F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t4 || istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(generate_edge_type(:l3F, :l2F, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:dprime))];
         :l2F)
    end

    @eval begin
        map_edges = Dict{Location, Dict{Location, Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location, Vector{$(edge_type)}}()
        end

        # l0G loc
        # l0G => l1G
        edge1 = $(generate_edge_type(:l0G, :l1G, 1))(nothing)
        map_edges[:l0G][:l1G] = [edge1]

        # l1G => l3G
        edge1 = $(generate_edge_type(:l1G, :l3G, 1))(nothing)
        edge2 = $(generate_edge_type(:l1G, :l3G, 2))(nothing)
        edge3 = $(generate_edge_type(:l1G, :l3G, 3))(nothing)
        edge4 = $(generate_edge_type(:l1G, :l3G, 4))(nothing)
        map_edges[:l1G][:l3G] = [edge1, edge2, edge3, edge4]

        # l1G => l4G
        edge1 = $(generate_edge_type(:l1G, :l4G, 1))(nothing)
        edge2 = $(generate_edge_type(:l1G, :l4G, 2))(nothing)
        map_edges[:l1G][:l4G] = [edge1, edge2]

        # l1G => l2G
        #=
        edge1 = $(generate_edge_type(:l1G, :l2G, 1))(nothing)
        edge2 = $(generate_edge_type(:l1G, :l2G, 2))(nothing)
        edge3 = $(generate_edge_type(:l1G, :l2G, 3))(nothing)
        edge4 = $(generate_edge_type(:l1G, :l2G, 4))(nothing)
        map_edges[:l1G][:l2G] = [edge3, edge4, edge1, edge2]
        =#

        # l3G loc
        # l3G => l1G
        edge1 = $(generate_edge_type(:l3G, :l1G, 1))([:ALL])
        map_edges[:l3G][:l1G] = [edge1]

        # l3G => l2G
        edge1 = $(generate_edge_type(:l3G, :l2G, 1))(nothing)
        edge2 = $(generate_edge_type(:l3G, :l2G, 2))(nothing)
        map_edges[:l3G][:l2G] = [edge1, edge2]

        # l4 loc
        # l4G => l1G
        edge1 = $(generate_edge_type(:l4G, :l1G, 1))([:ALL])
        map_edges[:l4G][:l1G] = [edge1]

        # l4G => l2G
        edge1 = $(generate_edge_type(:l4G, :l2G, 1))(nothing)
        edge2 = $(generate_edge_type(:l4G, :l2G, 2))(nothing)
        map_edges[:l4G][:l2G] = [edge1,edge2]

        # l2G loc
        # l2G => l1F : Transition from autF to autG
        edge1 = $(generate_edge_type(:l2G, :l1F, 1))(nothing)
        map_edges[:l2G][:l1F] = [edge1]

        # l1F loc
        # l1F => l3F
        edge1 = $(generate_edge_type(:l1F, :l2F, 1))(nothing)
        edge2 = $(generate_edge_type(:l1F, :l2F, 2))(nothing)
        map_edges[:l1F][:l2F] = [edge1, edge2]
        #edge3 = $(generate_edge_type(:l1F, :l2F, 3))(nothing)
        #edge4 = $(generate_edge_type(:l1F, :l2F, 4))(nothing)
        #map_edges[:l1F][:l2F] = [edge1, edge4, edge3, edge2]

        # l1F => l3F
        edge1 = $(generate_edge_type(:l1F, :l3F, 1))(nothing)
        edge2 = $(generate_edge_type(:l1F, :l3F, 2))(nothing)
        edge3 = $(generate_edge_type(:l1F, :l3F, 3))(nothing)
        map_edges[:l1F][:l3F] = [edge1, edge2, edge3]

        # l3F loc
        # l3F => l1F
        edge1 = $(generate_edge_type(:l3F, :l1F, 1))([:ALL])
        map_edges[:l3F][:l1F] = [edge1]

        # l3F => l2F
        edge1 = $(generate_edge_type(:l3F, :l2F, 1))(nothing)    
        map_edges[:l3F][:l2F] = [edge1]
    end

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2,
                                     :x3 => x3,  :x4 => x4, :t3 => t3, :t4 => t4)

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type, check_constraints, update_state!))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = AutomatonGandF(m.transitions, locations, Λ_F, locations_init, locations_final, 
                       map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

