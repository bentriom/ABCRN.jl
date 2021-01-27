
function create_period_automaton(m::ContinuousTimeModel, L::Float64, H::Float64, N::Int, sym_obs::VariableModel;
                                 initT::Float64 = 0.0)
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
                                                        :mean_tp => 5, :var_tp => 6)

    flow = Dict{Location,Vector{Float64}}(:l0 => [1.0,0.0,0.0,0.0,0.0,0.0],
                                          :l0prime => [1.0,0.0,0.0,0.0,0.0,0.0],
                                          :low => [1.0,0.0,0.0,1.0,0.0,0.0],
                                          :mid => [1.0,0.0,0.0,1.0,0.0,0.0],
                                          :high => [1.0,0.0,0.0,1.0,0.0,0.0],
                                          :final => [1.0,0.0,0.0,0.0,0.0,0.0]) 
    ## Edges
    map_edges = Dict{Location, Dict{Location, Vector{Edge}}}()
    for loc in locations 
        map_edges[loc] = Dict{Location, Vector{Edge}}()
    end
    
    idx_obs_var = getfield(m, :map_var_idx)[sym_obs]
    idx_var_t = map_var_automaton_idx[:t] 
    idx_var_tp = map_var_automaton_idx[:tp] 
    idx_var_mean_tp = map_var_automaton_idx[:mean_tp] 
    idx_var_var_tp = map_var_automaton_idx[:var_tp] 
    idx_var_n = map_var_automaton_idx[:n] 
    idx_var_top = map_var_automaton_idx[:top]
    
    nbr_rand = rand(1:100000)
    basename_func = "$(replace(m.name, ' '=>'_'))_$(nbr_rand)"
    basename_func = replace(basename_func, '-'=>'_')
    func_name(type_func::Symbol, from_loc::Location, to_loc::Location, edge_number::Int) = 
    Symbol("$(type_func)_aut_per_$(basename_func)_$(from_loc)$(to_loc)_$(edge_number)$(type_func == :us ? "!" : "")")
    meta_elementary_functions = quote
        ## Edge functions

        # l0 loc 
        # * l0 => l0
        @everywhere $(func_name(:cc, :l0, :l0, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0, :l0, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * l0 => l0prime
        @everywhere $(func_name(:cc, :l0, :l0prime, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:t] >= $initT
        @everywhere $(func_name(:us, :l0, :l0prime, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("l0prime")))

        # * l0 => low
        @everywhere $(func_name(:cc, :l0, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:t] >= $initT
        @everywhere $(func_name(:us, :l0, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low"));
         setindex!(getfield(S, :values), 0.0, $(idx_var_t));
         setindex!(getfield(S, :values), 0.0, $(idx_var_top));
         setindex!(getfield(S, :values), -1, $(idx_var_n)))

        # l0prime
        # * l0prime => l0prime
        @everywhere $(func_name(:cc, :l0prime, :l0prime, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0prime, :l0prime, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * l0prime => low
        @everywhere $(func_name(:cc, :l0prime, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        true
        @everywhere $(func_name(:us, :l0prime, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low"));
         setindex!(getfield(S, :values), 0.0, $(idx_var_t));
         setindex!(getfield(S, :values), 0.0, $(idx_var_top));
         setindex!(getfield(S, :values), -1, $(idx_var_n)))

        # low 
        # * low => low
        @everywhere $(func_name(:cc, :low, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :low, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * low => mid 
        @everywhere $(func_name(:cc, :low, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :low, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("mid")))

        # * low => final
        @everywhere $(func_name(:cc, :low, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] == $N
        @everywhere $(func_name(:us, :low, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("final")))

        # mid
        # * mid => mid
        @everywhere $(func_name(:cc, :mid, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :mid, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * mid => low 
        @everywhere $(func_name(:cc, :mid, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N &&
        S[:top] == 0.0
        @everywhere $(func_name(:us, :mid, :low, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low")))

        @everywhere $(func_name(:cc, :mid, :low, 2))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] == -1.0 &&
        S[:top] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 2))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low"));
         setindex!(getfield(S, :values), S[:n] + 1, $(idx_var_n));
         setindex!(getfield(S, :values), 0.0, $(idx_var_t));
         setindex!(getfield(S, :values), 0.0, $(idx_var_top)))

        @everywhere $(func_name(:cc, :mid, :low, 3))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        (0 <= S[:n] <= 1) &&
        S[:top] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 3))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low"));
         setindex!(getfield(S, :values), S[:n] + 1, $(idx_var_n));
         setindex!(getfield(S, :values), 0.0, $(idx_var_top));
         setindex!(getfield(S, :values), (S[:mean_tp] * (S[:n]-1) + S[:tp]) / S[:n], $(idx_var_mean_tp));
         setindex!(getfield(S, :values), 0.0, $(idx_var_tp)))

        @everywhere $(func_name(:cc, :mid, :low, 4))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        (2 <= S[:n] < $N) &&
        S[:top] == 1.0
        @everywhere $(func_name(:us, :mid, :low, 4))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("low"));
         setindex!(getfield(S, :values), S[:n] + 1, $(idx_var_n));
         setindex!(getfield(S, :values), 0.0, $(idx_var_top));
         setindex!(getfield(S, :values), (S[:mean_tp] * (S[:n]-1) + S[:tp]) / S[:n], $(idx_var_mean_tp));
         setindex!(getfield(S, :values), (S[:var_tp] * (S[:n]-1) + (S[:mean_tp]-S[:tp])^2) / S[:n], $(idx_var_var_tp)))

        # * mid => high
        @everywhere $(func_name(:cc, :mid, :high, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :mid, :high, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("high"));
         setindex!(getfield(S, :values), 1.0, $(idx_var_top)))

        # * mid => final
        @everywhere $(func_name(:cc, :mid, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] == $N
        @everywhere $(func_name(:us, :mid, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("final")))

        # high 
        # * high => high
        @everywhere $(func_name(:cc, :high, :high, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :high, :high, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (nothing)

        # * high => mid
        @everywhere $(func_name(:cc, :high, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] < $N
        @everywhere $(func_name(:us, :high, :mid, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("mid")))

        # * high => final
        @everywhere $(func_name(:cc, :high, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) = 
        S[:n] == $N
        @everywhere $(func_name(:us, :high, :final, 1))(S::StateLHA, x::Vector{Int}, p::Vector{Float64}) =
        (setfield!(S, :loc, Symbol("final")))
    end
    eval(meta_elementary_functions)

    # l0 loc 
    # * l0 => l0
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :l0, :l0, 1)), getfield(Main, func_name(:us, :l0, :l0, 1)))
    map_edges[:l0][:l0] = [edge_1]
    # * l0 => l0prime
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :l0, :l0prime, 1)), getfield(Main, func_name(:us, :l0, :l0prime, 1)))
    map_edges[:l0][:l0prime] = [edge_1]
    # * l0 => low
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :l0, :low, 1)), getfield(Main, func_name(:us, :l0, :low, 1)))
    map_edges[:l0][:low] = [edge_1]

    # l0prime
    # * l0prime => l0prime
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :l0prime, :l0prime, 1)), getfield(Main, func_name(:us, :l0prime, :l0prime, 1)))
    map_edges[:l0prime][:l0prime] = [edge_1]
    # * l0prime => low
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :l0prime, :low, 1)), getfield(Main, func_name(:us, :l0prime, :low, 1)))
    map_edges[:l0prime][:low] = [edge_1]

    # low 
    # * low => low
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :low, :low, 1)), getfield(Main, func_name(:us, :low, :low, 1)))
    map_edges[:low][:low] = [edge_1]
    # * low => mid 
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :low, :mid, 1)), getfield(Main, func_name(:us, :low, :mid, 1)))
    map_edges[:low][:mid] = [edge_1]
    # * low => final
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :low, :final, 1)), getfield(Main, func_name(:us, :low, :final, 1)))
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
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :mid, :final, 1)), getfield(Main, func_name(:us, :mid, :final, 1)))
    map_edges[:mid][:final] = [edge_1]

    # high 
    # * high => high
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :high, :high, 1)), getfield(Main, func_name(:us, :high, :high, 1)))
    map_edges[:high][:high] = [edge_1]
    # * high => mid
    edge_1 = Edge([:ALL], getfield(Main, func_name(:cc, :high, :mid, 1)), getfield(Main, func_name(:us, :high, :mid, 1)))
    map_edges[:high][:mid] = [edge_1]
    # * high => final
    edge_1 = Edge([nothing], getfield(Main, func_name(:cc, :high, :final, 1)), getfield(Main, func_name(:us, :high, :final, 1)))
    map_edges[:high][:final] = [edge_1] 

    ## Constants
    constants = Dict{Symbol,Float64}(:N => N, :L => L, :H => H, :initT => initT)

    A = LHA("Period", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_period_automaton

