
@everywhere f_mean_tp(mean_tp::Float64, tp::Float64, n::Float64) =
(mean_tp * (n-1) + tp) / n
@everywhere g_var_tp(var_tp::Float64, mean_tp::Float64, tp::Float64, n::Float64) =
(n-2)/(n-1)*var_tp + (tp-mean_tp)^2/n
#((n-2)*var_tp + (tp-mean_tp)*(tp - f_mean_tp(mean_tp, tp, n))) / (n-1)

@everywhere mean_error(mean_tp::Float64, var_tp::Float64, ref_mean_tp::Float64, ref_var_tp::Float64) =
abs(mean_tp - ref_mean_tp)

@everywhere min_mean_var_relative_error(mean_tp::Float64, var_tp::Float64, ref_mean_tp::Float64, ref_var_tp::Float64) =
min(abs((mean_tp - ref_mean_tp)/ref_mean_tp), sqrt(var_tp)/ref_mean_tp)

# Creation of the automaton types
#@everywhere @eval abstract type EdgePeriodAutomaton <: Edge end
@everywhere struct EdgePeriodAutomaton <: Edge 
    transitions::TransitionSet 
    check_constraints::CheckConstraintsFunction
    update_state!::UpdateStateFunction
end
@everywhere @eval $(MarkovProcesses.generate_code_lha_type_def(:PeriodAutomaton, :EdgePeriodAutomaton))

function create_period_automaton(m::ContinuousTimeModel, L::Float64, H::Float64, N::Int, sym_obs::VariableModel;
                                 initT::Float64 = 0.0, ref_mean_tp::Float64 = 0.0, ref_var_tp::Float64 = 0.0, error_func::Symbol = :mean_error)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert (L < H) "L >= H impossible for period automaton."
    @assert (N >= 1) "N < 1 impossible for period automaton."

    N = convert(Float64, N)
    # Automaton types and functions
    model_name = Symbol(typeof(m))
    lha_name = :PeriodAutomaton
    edge_type = :EdgePeriodAutomaton
    check_constraints = Symbol("check_constraints_$(lha_name)")
    update_state! = Symbol("update_state_$(lha_name)!")

    ## Locations
    locations = [:l0, :l0prime, :low, :mid, :high, :final]

    ## Invariant predicates
    idx_sym_obs = getfield(m, :map_var_idx)[sym_obs]
    id = MarkovProcesses.newid()
    basename_func = "$(model_name)_$(id)"
    sym_name_L = Symbol("val_L_aut_per_$(basename_func)")
    sym_name_H = Symbol("val_H_aut_per_$(basename_func)")

    @everywhere true_predicate(x::Vector{Int}) = true
    @everywhere low_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] <= $L
    @everywhere not_low_predicate(x::Vector{Int}) = !low_predicate(x)
    @everywhere mid_predicate(x::Vector{Int}) = $L < x[$(Meta.quot(idx_sym_obs))] < $H
    @everywhere high_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] >= $H

    Λ_F = Dict{Location,InvariantPredicateFunction}(:l0 => getfield(Main, :true_predicate), :l0prime => getfield(Main, :not_low_predicate),
                                                    :low => getfield(Main, :low_predicate), :mid => getfield(Main, :mid_predicate), 
                                                    :high => getfield(Main, :high_predicate), :final => getfield(Main, :true_predicate))

    ## Init and final loc
    locations_init = [:l0]
    locations_final = [:final]

    ## Map of automaton variables
    map_var_automaton_idx = Dict{VariableAutomaton,Int}(:t => 1, :n => 2, :top => 3, :tp => 4,
                                                        :mean_tp => 5, :var_tp => 6, :d => 7)

    flow = Dict{Location,Vector{Float64}}(:l0 => [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
                                          :l0prime => [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
                                          :low => [1.0,0.0,0.0,1.0,0.0,0.0,0.0],
                                          :mid => [1.0,0.0,0.0,1.0,0.0,0.0,0.0],
                                          :high => [1.0,0.0,0.0,1.0,0.0,0.0,0.0],
                                          :final => [1.0,0.0,0.0,0.0,0.0,0.0,0.0]) 
    ## Edges
    to_idx(var::Symbol) = map_var_automaton_idx[var]
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    edge_name(from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("Edge_$(lha_name)_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)")
    function check_constraints(from_loc::Location, to_loc::Location, edge_number::Int)
        return Symbol("check_constraints_$(edge_type)_$(from_loc)$(to_loc)_$(edge_number)_$(model_name)_$(id)")
    end
    function update_state!(from_loc::Location, to_loc::Location, edge_number::Int)
        return Symbol("update_state_$(edge_type)_$(from_loc)$(to_loc)_$(edge_number)_$(model_name)_$(id)!")
    end

    ## check_constraints & update_state!
    meta_funcs = quote
        # l0 loc 
        # * l0 => l0
        #struct $(edge_name(:l0, :l0, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0, :l0, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(update_state!(:l0, :l0, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:l0)

        # * l0 => l0prime
        #struct $(edge_name(:l0, :l0prime, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $initT
        @everywhere $(update_state!(:l0, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:d))] = Inf;
         :l0prime)

        # * l0 => low
        #struct $(edge_name(:l0, :low, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $initT
        @everywhere $(update_state!(:l0, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:t))] = 0.0;
         S_values[$(to_idx(:top))] = 0.0;
         S_values[$(to_idx(:n))] = -1;
         S_values[$(to_idx(:tp))] = 0.0;
         S_values[$(to_idx(:d))] = Inf;
         :low)

        # l0prime
        # * l0prime => l0prime
        #struct $(edge_name(:l0prime, :l0prime, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0prime, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(update_state!(:l0prime, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:l0prime)

        # * l0prime => low
        #struct $(edge_name(:l0prime, :low, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:l0prime, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(update_state!(:l0prime, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:t))] = 0.0;
         S_values[$(to_idx(:top))] = 0.0;
         S_values[$(to_idx(:n))] = -1;
         S_values[$(to_idx(:tp))] = 0.0;
         :low)

        # low 
        # * low => low
        #struct $(edge_name(:low, :low, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:low, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:low, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:low)

        # * low => mid 
        #struct $(edge_name(:low, :mid, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:low, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:low, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:mid)

        # * low => final
        #struct $(edge_name(:low, :final, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:low, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(update_state!(:low, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (val_d = $(error_func)(S_values[$(to_idx(:mean_tp))], 
                               S_values[$(to_idx(:var_tp))], 
                               $(ref_mean_tp), $(ref_var_tp));
         S_values[$(to_idx(:d))] = val_d;
         :final)

        # mid
        # * mid => mid
        #struct $(edge_name(:mid, :mid, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:mid, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:mid)

        # * mid => low 
        #struct $(edge_name(:mid, :low, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N &&
        S_values[$(to_idx(:top))] == 0.0
        @everywhere $(update_state!(:mid, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:low)

        #struct $(edge_name(:mid, :low, 2)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :low, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == -1.0 &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(update_state!(:mid, :low, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:n))] = S_values[$(to_idx(:n))] + 1;
         S_values[$(to_idx(:top))] = 0.0;
         S_values[$(to_idx(:tp))] = 0.0;
         :low)

        #struct $(edge_name(:mid, :low, 3)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :low, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] == 0.0) &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(update_state!(:mid, :low, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:n))] = S_values[$(to_idx(:n))] + 1;
         S_values[$(to_idx(:top))] = 0.0;
         S_values[$(to_idx(:mean_tp))] = f_mean_tp(S_values[$(to_idx(:mean_tp))], 
                                                   S_values[$(to_idx(:tp))],
                                                   S_values[$(to_idx(:n))]);
         S_values[$(to_idx(:tp))] = 0.0;
         :low)

        #struct $(edge_name(:mid, :low, 4)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :low, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (1 <= S_values[$(to_idx(:n))] < $N) &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(update_state!(:mid, :low, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:n))] = S_values[$(to_idx(:n))] + 1;
         S_values[$(to_idx(:top))] = 0.0;
         S_values[$(to_idx(:var_tp))] = g_var_tp(S_values[$(to_idx(:var_tp))], 
                                                 S_values[$(to_idx(:mean_tp))],
                                                 S_values[$(to_idx(:tp))],
                                                 S_values[$(to_idx(:n))]);
         S_values[$(to_idx(:mean_tp))] = f_mean_tp(S_values[$(to_idx(:mean_tp))], 
                                                   S_values[$(to_idx(:tp))],
                                                   S_values[$(to_idx(:n))]);
         S_values[$(to_idx(:tp))] = 0.0;
         :low)

        # * mid => high
        #struct $(edge_name(:mid, :high, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:mid, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (S_values[$(to_idx(:top))] = 1.0;
         :high)

        # * mid => final
        #struct $(edge_name(:mid, :final, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:mid, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(update_state!(:mid, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (val_d = $(error_func)(S_values[$(to_idx(:mean_tp))],
                               S_values[$(to_idx(:var_tp))],
                               $(ref_mean_tp), $(ref_var_tp));
         S_values[$(to_idx(:d))] = val_d;
         :final)

        # high 
        # * high => high
        #struct $(edge_name(:high, :high, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:high, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:high, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:high)

        # * high => mid
        #struct $(edge_name(:high, :mid, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:high, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(update_state!(:high, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (:mid)

        # * high => final
        #struct $(edge_name(:high, :final, 1)) <: $(edge_type) transitions::Union{Nothing,Vector{Symbol}} end
        @everywhere $(check_constraints(:high, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(update_state!(:high, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (val_d = $(error_func)(S_values[$(to_idx(:mean_tp))],
                               S_values[$(to_idx(:var_tp))],
                               $(ref_mean_tp), $(ref_var_tp));
         S_values[$(to_idx(:d))] = val_d;
         :final)
    end
    eval(meta_funcs)

    @eval begin
        map_edges = Dict{Location, Dict{Location, Vector{$(edge_type)}}}()
        for loc in $(locations)
            map_edges[loc] = Dict{Location, Vector{$(edge_type)}}()
        end

        # l0 loc 
        # * l0 => l0
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:l0, :l0, 1)), $(update_state!(:l0, :l0, 1)))
        map_edges[:l0][:l0] = [edge1]
        # * l0 => l0prime
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:l0, :l0prime, 1)), $(update_state!(:l0, :l0prime, 1)))
        map_edges[:l0][:l0prime] = [edge1]
        # * l0 => low
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:l0, :low, 1)), $(update_state!(:l0, :low, 1)))
        map_edges[:l0][:low] = [edge1]

        # l0prime
        # * l0prime => l0prime
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:l0prime, :l0prime, 1)), $(update_state!(:l0prime, :l0prime, 1)))
        map_edges[:l0prime][:l0prime] = [edge1]
        # * l0prime => low
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:l0prime, :low, 1)), $(update_state!(:l0prime, :low, 1)))
        map_edges[:l0prime][:low] = [edge1]

        # low 
        # * low => low
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:low, :low, 1)), $(update_state!(:low, :low, 1)))
        map_edges[:low][:low] = [edge1]
        # * low => mid 
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:low, :mid, 1)), $(update_state!(:low, :mid, 1)))
        map_edges[:low][:mid] = [edge1]
        # * low => final
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:low, :final, 1)), $(update_state!(:low, :final, 1)))
        map_edges[:low][:final] = [edge1]

        # mid
        # * mid => mid
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :mid, 1)), $(update_state!(:mid, :mid, 1)))
        map_edges[:mid][:mid] = [edge1]
        # * mid => low 
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :low, 1)), $(update_state!(:mid, :low, 1)))
        edge2 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :low, 2)), $(update_state!(:mid, :low, 2)))
        edge3 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :low, 3)), $(update_state!(:mid, :low, 3)))
        edge4 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :low, 4)), $(update_state!(:mid, :low, 4)))
        map_edges[:mid][:low] = [edge1, edge2, edge3, edge4]
        # * mid => high
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:mid, :high, 1)), $(update_state!(:mid, :high, 1)))
        map_edges[:mid][:high] = [edge1]
        # * mid => final
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:mid, :final, 1)), $(update_state!(:mid, :final, 1)))
        map_edges[:mid][:final] = [edge1]

        # high 
        # * high => high
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:high, :high, 1)), $(update_state!(:high, :high, 1)))
        map_edges[:high][:high] = [edge1]
        # * high => mid
        edge1 = EdgePeriodAutomaton([:ALL], $(check_constraints(:high, :mid, 1)), $(update_state!(:high, :mid, 1)))
        map_edges[:high][:mid] = [edge1]
        # * high => final
        edge1 = EdgePeriodAutomaton(nothing, $(check_constraints(:high, :final, 1)), $(update_state!(:high, :final, 1)))
        map_edges[:high][:final] = [edge1] 
    end

    ## Create data separately
    map_edges_transitions = Dict{Symbol, Dict{Symbol,Vector{TransitionSet}}}()
    map_edges_check_constraints = Dict{Symbol, Dict{Symbol,Vector{CheckConstraintsFunction}}}()
    map_edges_update_state = Dict{Symbol, Dict{Symbol,Vector{UpdateStateFunction}}}()

    ## Constants
    constants = Dict{Symbol,Float64}(:N => N, :L => L, :H => H, :initT => initT)

    # Updating types and simulation methods
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_model_type_def(model_name, lha_name))
    @everywhere @eval $(MarkovProcesses.generate_code_next_state(lha_name, edge_type))
    @everywhere @eval $(MarkovProcesses.generate_code_synchronized_simulation(model_name, lha_name, edge_type, m.f!, m.isabsorbing))

    A = PeriodAutomaton(m.transitions, locations, Λ_F, locations_init, locations_final, 
                        map_var_automaton_idx, flow, map_edges, 
                        map_edges_transitions, map_edges_check_constraints, map_edges_update_state,
                        constants, m.map_var_idx)
    return A
end

export create_period_automaton

