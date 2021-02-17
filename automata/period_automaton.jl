#(S[:mean_tp] * (S[:n]-1) + S[:tp]) / S[:n]
#(S[:var_tp] * (S[:n]-1) + (S[:mean_tp]-S[:tp])^2) / S[:n]

@everywhere f_mean_tp(mean_tp::Float64, tp::Float64, n::Float64) =
(mean_tp * (n-1) + tp) / n
@everywhere g_var_tp(var_tp::Float64, mean_tp::Float64, tp::Float64, n::Float64) =
(n-2)/(n-1)*var_tp + (tp-mean_tp)^2/n
#((n-2)*var_tp + (tp-mean_tp)*(tp - f_mean_tp(mean_tp, tp, n))) / (n-1)

@everywhere mean_error(mean_tp::Float64, var_tp::Float64, ref_mean_tp::Float64, ref_var_tp::Float64) =
abs(mean_tp - ref_mean_tp)

function create_period_automaton(m::ContinuousTimeModel, L::Float64, H::Float64, N::Int, sym_obs::VariableModel;
                                 initT::Float64 = 0.0, ref_mean_tp::Float64 = 0.0, ref_var_tp::Float64 = 0.0, error_func::Symbol = :mean_error)
    # Requirements for the automaton
    @assert sym_obs in m.g "$(sym_obs) is not observed."
    @assert (L < H) "L >= H impossible for period automaton."
    @assert (N >= 1) "N < 1 impossible for period automaton."
    
    N = convert(Float64, N)
    nbr_rand = rand(1:1000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')

    ## Locations
    locations = [:l0, :l0prime, :low, :mid, :high, :final]

    ## Invariant predicates
    idx_sym_obs = getfield(m, :map_var_idx)[sym_obs]
    sym_name_L = Symbol("val_L_aut_per_$(basename_func)")
    sym_name_H = Symbol("val_H_aut_per_$(basename_func)")
    
    @everywhere true_predicate(x::Vector{Int}) = true
    @everywhere low_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] <= $L
    @everywhere not_low_predicate(x::Vector{Int}) = !low_predicate(x)
    @everywhere mid_predicate(x::Vector{Int}) = $L < x[$(Meta.quot(idx_sym_obs))] < $H
    @everywhere high_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] >= $H

    Λ_F = Dict(:l0 => getfield(Main, :true_predicate), :l0prime => getfield(Main, :not_low_predicate),
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
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end
 
    to_idx(var::Symbol) = map_var_automaton_idx[var] 
    nbr_rand = rand(1:100000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_aut_per_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote
        ## Edge functions

        # l0 loc 
        # * l0 => l0
        @everywhere $(func_name(:cc, :l0, :l0, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0, :l0, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * l0 => l0prime
        @everywhere $(func_name(:cc, :l0, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $initT
        @everywhere $(func_name(:us, :l0, :l0prime, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("l0prime");
         setindex!(S_values, Inf, $(to_idx(:d))))

        # * l0 => low
        @everywhere $(func_name(:cc, :l0, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:t))] >= $initT
        @everywhere $(func_name(:us, :l0, :low, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low");
         setindex!(S_values, 0.0, $(to_idx(:t)));
         setindex!(S_values, 0.0, $(to_idx(:top)));
         setindex!(S_values, -1, $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:tp)));
         setindex!(S_values, Inf, $(to_idx(:d))))

        # l0prime
        # * l0prime => l0prime
        @everywhere $(func_name(:cc, :l0prime, :l0prime, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0prime, :l0prime, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * l0prime => low
        @everywhere $(func_name(:cc, :l0prime, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0prime, :low, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low");
         setindex!(S_values, 0.0, $(to_idx(:t)));
         setindex!(S_values, 0.0, $(to_idx(:top)));
         setindex!(S_values, -1, $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:tp))))

        # low 
        # * low => low
        @everywhere $(func_name(:cc, :low, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :low, :low, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * low => mid 
        @everywhere $(func_name(:cc, :low, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :low, :mid, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("mid"))

        # * low => final
        @everywhere $(func_name(:cc, :low, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(func_name(:us, :low, :final, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("final");
         val_d = getfield(Main, $(Meta.quot(error_func)))(S_values[$(to_idx(:mean_tp))], 
                                                          S_values[$(to_idx(:var_tp))], 
                                                          $(ref_mean_tp), $(ref_var_tp));
         setindex!(S_values, val_d, $(to_idx(:d))))

        # mid
        # * mid => mid
        @everywhere $(func_name(:cc, :mid, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :mid, :mid, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * mid => low 
        @everywhere $(func_name(:cc, :mid, :low, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N &&
        S_values[$(to_idx(:top))] == 0.0
        @everywhere $(func_name(:us, :mid, :low, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low"))

        @everywhere $(func_name(:cc, :mid, :low, 2))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == -1.0 &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 2))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low");
         setindex!(S_values, S_values[$(to_idx(:n))] + 1, $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:top)));
         setindex!(S_values, 0.0, $(to_idx(:tp))))

        @everywhere $(func_name(:cc, :mid, :low, 3))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (S_values[$(to_idx(:n))] == 0.0) &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 3))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low");
         setindex!(S_values, S_values[$(to_idx(:n))] + 1, $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:top)));
         setindex!(S_values, f_mean_tp(S_values[$(to_idx(:mean_tp))], 
                                                   S_values[$(to_idx(:tp))],
                                                   S_values[$(to_idx(:n))]), $(to_idx(:mean_tp)));
         setindex!(S_values, 0.0, $(to_idx(:tp))))

        @everywhere $(func_name(:cc, :mid, :low, 4))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        (1 <= S_values[$(to_idx(:n))] < $N) &&
        S_values[$(to_idx(:top))] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 4))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("low");
         setindex!(S_values, S_values[$(to_idx(:n))] + 1, $(to_idx(:n)));
         setindex!(S_values, 0.0, $(to_idx(:top)));
         setindex!(S_values, g_var_tp(S_values[$(to_idx(:var_tp))], 
                                                   S_values[$(to_idx(:mean_tp))],
                                                   S_values[$(to_idx(:tp))],
                                                   S_values[$(to_idx(:n))]), $(to_idx(:var_tp)));
         setindex!(S_values, f_mean_tp(S_values[$(to_idx(:mean_tp))], 
                                                   S_values[$(to_idx(:tp))],
                                                   S_values[$(to_idx(:n))]), $(to_idx(:mean_tp)));
         setindex!(S_values, 0.0, $(to_idx(:tp))))

        # * mid => high
        @everywhere $(func_name(:cc, :mid, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :mid, :high, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("high");
         setindex!(S_values, 1.0, $(to_idx(:top))))

        # * mid => final
        @everywhere $(func_name(:cc, :mid, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(func_name(:us, :mid, :final, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("final");
         val_d = getfield(Main, Meta.quot($error_func))(S_values[$(to_idx(:mean_tp))],
                                                        S_values[$(to_idx(:var_tp))],
                                                        $(ref_mean_tp), $(ref_var_tp));
         setindex!(S_values, val_d, $(to_idx(:d))))

        # high 
        # * high => high
        @everywhere $(func_name(:cc, :high, :high, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :high, :high, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * high => mid
        @everywhere $(func_name(:cc, :high, :mid, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] < $N
        @everywhere $(func_name(:us, :high, :mid, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("mid"))

        # * high => final
        @everywhere $(func_name(:cc, :high, :final, 1))(S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) = 
        S_values[$(to_idx(:n))] == $N
        @everywhere $(func_name(:us, :high, :final, 1))(ptr_loc::Vector{Symbol}, S_time::Float64, S_values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64}) =
        (ptr_loc[1] = Symbol("final");
         val_d = getfield(Main, Meta.quot($error_func))(S_values[$(to_idx(:mean_tp))],
                                                        S_values[$(to_idx(:var_tp))],
                                                        $(ref_mean_tp), $(ref_var_tp));
         setindex!(S_values, val_d, $(to_idx(:d))))
    end
    eval(meta_elementary_functions)

    # l0 loc 
    # * l0 => l0
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :l0, :l0, 1)), getfield(Main, func_name(:us, :l0, :l0, 1)))
    map_edges[:l0][:l0] = [edge_1]
    # * l0 => l0prime
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :l0prime, 1)), getfield(Main, func_name(:us, :l0, :l0prime, 1)))
    map_edges[:l0][:l0prime] = [edge_1]
    # * l0 => low
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :l0, :low, 1)), getfield(Main, func_name(:us, :l0, :low, 1)))
    map_edges[:l0][:low] = [edge_1]

    # l0prime
    # * l0prime => l0prime
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :l0prime, :l0prime, 1)), getfield(Main, func_name(:us, :l0prime, :l0prime, 1)))
    map_edges[:l0prime][:l0prime] = [edge_1]
    # * l0prime => low
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :l0prime, :low, 1)), getfield(Main, func_name(:us, :l0prime, :low, 1)))
    map_edges[:l0prime][:low] = [edge_1]

    # low 
    # * low => low
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :low, :low, 1)), getfield(Main, func_name(:us, :low, :low, 1)))
    map_edges[:low][:low] = [edge_1]
    # * low => mid 
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :low, :mid, 1)), getfield(Main, func_name(:us, :low, :mid, 1)))
    map_edges[:low][:mid] = [edge_1]
    # * low => final
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :low, :final, 1)), getfield(Main, func_name(:us, :low, :final, 1)))
    map_edges[:low][:final] = [edge_1]

    # mid
    # * mid => mid
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :mid, 1)), getfield(Main, func_name(:us, :mid, :mid, 1)))
    map_edges[:mid][:mid] = [edge_1]
    # * mid => low 
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :low, 1)), getfield(Main, func_name(:us, :mid, :low, 1)))
    edge_2 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :low, 2)), getfield(Main, func_name(:us, :mid, :low, 2)))
    edge_3 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :low, 3)), getfield(Main, func_name(:us, :mid, :low, 3)))
    edge_4 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :low, 4)), getfield(Main, func_name(:us, :mid, :low, 4)))
    map_edges[:mid][:low] = [edge_1, edge_2, edge_3, edge_4]
    # * mid => high
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :mid, :high, 1)), getfield(Main, func_name(:us, :mid, :high, 1)))
    map_edges[:mid][:high] = [edge_1]
    # * mid => final
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :mid, :final, 1)), getfield(Main, func_name(:us, :mid, :final, 1)))
    map_edges[:mid][:final] = [edge_1]

    # high 
    # * high => high
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :high, :high, 1)), getfield(Main, func_name(:us, :high, :high, 1)))
    map_edges[:high][:high] = [edge_1]
    # * high => mid
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :high, :mid, 1)), getfield(Main, func_name(:us, :high, :mid, 1)))
    map_edges[:high][:mid] = [edge_1]
    # * high => final
    edge_1 = Edge(nothing, getfield(Main, func_name(:cc, :high, :final, 1)), getfield(Main, func_name(:us, :high, :final, 1)))
    map_edges[:high][:final] = [edge_1] 

    ## Constants
    constants = Dict{Symbol,Float64}(:N => N, :L => L, :H => H, :initT => initT)

    A = LHA("Period", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_period_automaton

