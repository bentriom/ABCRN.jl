
# Creation of the automaton types
@everywhere @eval abstract type EdgeAutomatonG <: Edge end
@everywhere @eval $(MarkovProcesses.generate_code_lha_type_def(:AutomatonG, :EdgeAutomatonG))

function create_automaton_G(m::ContinuousTimeModel, x1::Float64, x2::Float64, t1::Float64, t2::Float64, sym_obs::VariableModel)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert (x1 <= x2) "x1 > x2 impossible for G automaton."
    @assert (t1 <= t2) "t1 > t2 impossible for G automaton."

    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :AutomatonG
    edge_type = :EdgeAutomatonG
    check_constraints = Symbol("check_constraints_$(lha_name)")
    update_state! = Symbol("update_state_$(lha_name)!")

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
    to_idx(var::Symbol) = map_var_automaton_idx[var]
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    id = MarkovProcesses.newid()
    basename_func = "$(model_name)_$(id)"
    edge_name(from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("Edge_$(lha_name)_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)")

    ## check_constraints & update_state!
    @everywhere @eval begin
        istrue(val::Float64) = convert(Bool, val)
        # l0 loc
        # l0 => l1
        struct $(edge_name(:l0, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l0, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(edge_name(:l0, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = 0;
         S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:in))] = true;
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p,x);
         :l1)

        # l1 => l3
        struct $(edge_name(:l1, :l3, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l3, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_time <= $t1 && 
        S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2
        $(update_state!)(edge::$(edge_name(:l1, :l3, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:in))] = false;
         :l3)

        struct $(edge_name(:l1, :l3, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l3, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time <= $t1) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(edge_name(:l1, :l3, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = 0;
         S_values[$(to_idx(:in))] = false;
         :l3)

        struct $(edge_name(:l1, :l3, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l3, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(edge_name(:l1, :l3, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * (S_time - $t1);
         S_values[$(to_idx(:tprime))] = 0.0;
         :l3)

        struct $(edge_name(:l1, :l3, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l3, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        ($x1 <= S_values[$(to_idx(:n))] <= $x2)
        $(update_state!)(edge::$(edge_name(:l1, :l3, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:tprime))] = 0.0;
         :l3)

        # l1 => l4
        struct $(edge_name(:l1, :l4, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l4, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2)
        $(update_state!)(edge::$(edge_name(:l1, :l4, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:d))] * (S_time - $t1);
         :l4)

        struct $(edge_name(:l1, :l4, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l4, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        ($t1 <= S_time <= $t2) && 
        (S_values[$(to_idx(:n))] < $x1 || S_values[$(to_idx(:n))] > $x2)
        $(update_state!)(edge::$(edge_name(:l1, :l4, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l4)

        # l1 => l2
        #=
        struct $(edge_name(:l1, :l2, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        S_time >= $t2
        $(update_state!)(edge::$(edge_name(:l1, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2)

        struct $(edge_name(:l1, :l2, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        S_time >= $t2
        $(update_state!)(edge::$(edge_name(:l1, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * ($t2 - $t1);
        :l2)

        struct $(edge_name(:l1, :l2, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l2, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:isabs))]) && 
        S_time <= $t1
        $(update_state!)(edge::$(edge_name(:l1, :l2, 3)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = ($t2 - $t1) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
        :l2)

        struct $(edge_name(:l1, :l2, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l1, :l2, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:isabs))]) && 
        ($t1 <= S_time <= $t2)
        $(update_state!)(edge::$(edge_name(:l1, :l2, 4)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + ($t2 - S_time) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
        :l2)
        =#

        # l3 => l1
        struct $(edge_name(:l3, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l3, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(edge_name(:l3, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p, x);
         :l1)

        # l4 => l1
        struct $(edge_name(:l4, :l1, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l4, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = true
        $(update_state!)(edge::$(edge_name(:l4, :l1, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:tprime))] * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:tprime))] = 0.0;
         S_values[$(to_idx(:n))] = x[$(idx_obs_var)];
         S_values[$(to_idx(:in))] = true;
         S_values[$(to_idx(:isabs))] = $(m.isabsorbing)(p,x);
         :l1)

        # l3 => l2
        struct $(edge_name(:l3, :l2, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l3, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        istrue(S_values[$(to_idx(:in))]) && 
        (S_time >= $t2 || istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(edge_name(:l3, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] * ($t2 - $t1);
         :l2)

        struct $(edge_name(:l3, :l2, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l3, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        !istrue(S_values[$(to_idx(:in))]) && 
        (S_time >= $t2 || istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(edge_name(:l3, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (:l2)

        # l4 => l2
        struct $(edge_name(:l4, :l2, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l4, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (istrue(S_values[$(to_idx(:isabs))]))
        $(update_state!)(edge::$(edge_name(:l4, :l2, 1)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + ($t2 - S_time) * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:tprime))] = 0.0;
         :l2)

        struct $(edge_name(:l4, :l2, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        $(check_constraints)(edge::$(edge_name(:l4, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_time >= $t2)
        $(update_state!)(edge::$(edge_name(:l4, :l2, 2)), S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:d))] = S_values[$(to_idx(:d))] + S_values[$(to_idx(:tprime))] * min(abs($x1 - S_values[$(to_idx(:n))]), abs($x2 - S_values[$(to_idx(:n))]));
         S_values[$(to_idx(:tprime))] = 0.0;
         :l2)

    end

    @eval begin
        map_edges = Dict{Location, Dict{Location, Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location, Vector{$(edge_type)}}()
        end

        # l0 loc
        # l0 => l1
        edge1 = $(edge_name(:l0, :l1, 1))(nothing)
        map_edges[:l0][:l1] = [edge1]

        # l1 loc
        # l1 => l3
        edge1 = $(edge_name(:l1, :l3, 1))(nothing)
        edge2 = $(edge_name(:l1, :l3, 2))(nothing)
        edge3 = $(edge_name(:l1, :l3, 3))(nothing)
        edge4 = $(edge_name(:l1, :l3, 4))(nothing)
        map_edges[:l1][:l3] = [edge1, edge2, edge3, edge4]

        # l1 => l4
        edge1 = $(edge_name(:l1, :l4, 1))(nothing)
        edge2 = $(edge_name(:l1, :l4, 2))(nothing)
        map_edges[:l1][:l4] = [edge1, edge2]

        # l1 => l2
        #=
        edge1 = $(edge_name(:l1, :l2, 1))(nothing)
        edge2 = $(edge_name(:l1, :l2, 2))(nothing)
        edge3 = $(edge_name(:l1, :l2, 3))(nothing)
        edge4 = $(edge_name(:l1, :l2, 4))(nothing)
        map_edges[:l1][:l2] = [edge1, edge2, edge3, edge4]
        =#

        # l3 loc
        # l3 => l1
        edge1 = $(edge_name(:l3, :l1, 1))([:ALL])
        map_edges[:l3][:l1] = [edge1]

        # l3 => l2
        edge1 = $(edge_name(:l3, :l2, 1))(nothing)
        edge2 = $(edge_name(:l3, :l2, 2))(nothing)
        map_edges[:l3][:l2] = [edge1, edge2]

        # l4 loc
        # l4 => l1
        edge1 = $(edge_name(:l4, :l1, 1))([:ALL])
        map_edges[:l4][:l1] = [edge1]

        # l4 => l2
        edge1 = $(edge_name(:l4, :l2, 1))(nothing)
        edge2 = $(edge_name(:l4, :l2, 2))(nothing)
        map_edges[:l4][:l2] = [edge1,edge2]
    end

    ## Constants
    constants = Dict{Symbol,Float64}(:x1 => x1,  :x2 => x2, :t1 => t1, :t2 => t2)

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type, check_constraints, update_state!))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = AutomatonG(m.transitions, locations, Λ_F, locations_init, locations_final, 
                   map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A

end

