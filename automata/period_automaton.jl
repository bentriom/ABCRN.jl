
## Edge functions

# l0 loc 
# * l0 => l0
@everywhere cc_aut_per_l0l0_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
true
@everywhere us_aut_per_l0l0_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(nothing)

# * l0 => l0prime
@everywhere cc_aut_per_l0l0prime_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:t] >= constants[:initT]
@everywhere us_aut_per_l0l0prime_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :l0prime)

# * l0 => low
@everywhere cc_aut_per_l0low_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:t] >= constants[:initT]
@everywhere us_aut_per_l0low_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low;
 S[:t] = 0.0;
 S[:top] = 0.0;
 S[:n] = -1)

# l0prime
# * l0prime => l0prime
@everywhere cc_aut_per_l0primel0prime_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
true
@everywhere us_aut_per_l0primel0prime_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(nothing)

# * l0prime => low
@everywhere cc_aut_per_l0primelow_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
true
@everywhere us_aut_per_l0primelow_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low;
 S[:t] = 0.0;
 S[:top] = 0.0;
 S[:n] = -1)

# low 
# * low => low
@everywhere cc_aut_per_lowlow_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_lowlow_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(nothing)

# * low => mid 
@everywhere cc_aut_per_lowmid_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_lowmid_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :mid)

# * low => final
@everywhere cc_aut_per_lowfinal_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] == constants[:N]
@everywhere us_aut_per_lowfinal_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :final)

# mid
# * mid => mid
@everywhere cc_aut_per_midmid_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_midmid_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(nothing)

# * mid => low 
@everywhere cc_aut_per_midlow_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N] &&
S[:top] == 0.0
@everywhere us_aut_per_midlow_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low)

@everywhere cc_aut_per_midlow_2(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] == -1.0 &&
S[:top] == 1.0
@everywhere us_aut_per_midlow_2!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low;
 S[:n] += 1;
 S[:t] = 0.0;
 S[:top] = 0.0)

@everywhere cc_aut_per_midlow_3(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(0 <= S[:n] <= 1) &&
S[:top] == 1.0
@everywhere us_aut_per_midlow_3!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low;
 S[:n] += 1;
 S[:top] = 0.0;
 S[:mean_tp] = (S[:mean_tp] * (S[:n]-1) + S[:tp]) / S[:n];
 S[:tp] = 0.0)

@everywhere cc_aut_per_midlow_4(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
(2 <= S[:n] < constants[:N]) &&
S[:top] == 1.0
@everywhere us_aut_per_midlow_4!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :low;
 S[:n] += 1;
 S[:top] = 0.0;
 S[:mean_tp] = (S[:mean_tp] * (S[:n]-1) + S[:tp]) / S[:n];
 S[:var_tp] = (S[:var_tp] * (S[:n]-1) + (S[:mean_tp]-S[:tp])^2) / S[:n];
 S[:tp] = 0.0)

# * mid => high
@everywhere cc_aut_per_midhigh_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_midhigh_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :high;
 S[:top] = 1.0)

# * mid => final
@everywhere cc_aut_per_midfinal_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] == constants[:N]
@everywhere us_aut_per_midfinal_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :final)

# high 
# * high => high
@everywhere cc_aut_per_highhigh_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_highhigh_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(nothing)

# * high => mid
@everywhere cc_aut_per_highmid_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] < constants[:N]
@everywhere us_aut_per_highmid_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :mid)

# * high => final
@everywhere cc_aut_per_highfinal_1(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) = 
S[:n] == constants[:N]
@everywhere us_aut_per_highfinal_1!(S::StateLHA, constants::Dict{Symbol,Float64}, x::Vector{Int}, p::Vector{Float64}) =
(S.loc = :final)

function create_period_automaton(m::ContinuousTimeModel, L::Float64, H::Float64, N::Int, sym_obs::VariableModel;
                                 initT::Float64 = 0.0)
    @assert sym_obs in m.g "$(sym_obs) is not observed."
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
    #=
    eval(Meta.parse("@everywhere true_predicate(x::Vector{Int}) = true"))
    eval(Meta.parse("@everywhere low_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] <= $L"))
    eval(Meta.parse("@everywhere not_low_predicate(x::Vector{Int}) = !low_predicate(x)"))
    eval(Meta.parse("@everywhere mid_predicate(x::Vector{Int}) = L < x[$(Meta.quot(idx_sym_obs))] < $H"))
    eval(Meta.parse("@everywhere high_predicate(x::Vector{Int}) = x[$(Meta.quot(idx_sym_obs))] >= $H"))
    =#

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

    # l0 loc 
    # * l0 => l0
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_l0l0_1), getfield(Main, :us_aut_per_l0l0_1!))
    map_edges[:l0][:l0] = [edge_1]
    # * l0 => l0prime
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_l0l0prime_1), getfield(Main, :us_aut_per_l0l0prime_1!))
    map_edges[:l0][:l0prime] = [edge_1]
    # * l0 => low
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_l0low_1), getfield(Main, :us_aut_per_l0low_1!))
    map_edges[:l0][:low] = [edge_1]

    # l0prime
    # * l0prime => l0prime
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_l0primel0prime_1), getfield(Main, :us_aut_per_l0primel0prime_1!))
    map_edges[:l0prime][:l0prime] = [edge_1]
    # * l0prime => low
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_l0primelow_1), getfield(Main, :us_aut_per_l0primelow_1!))
    map_edges[:l0prime][:low] = [edge_1]

    # low 
    # * low => low
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_lowlow_1), getfield(Main, :us_aut_per_lowlow_1!))
    map_edges[:low][:low] = [edge_1]
    # * low => mid 
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_lowmid_1), getfield(Main, :us_aut_per_lowmid_1!))
    map_edges[:low][:mid] = [edge_1]
    # * low => final
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_lowfinal_1), getfield(Main, :us_aut_per_lowfinal_1!))
    map_edges[:low][:final] = [edge_1]

    # mid
    # * mid => mid
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_midmid_1), getfield(Main, :us_aut_per_midmid_1!))
    map_edges[:mid][:mid] = [edge_1]
    # * mid => low 
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_midlow_1), getfield(Main, :us_aut_per_midlow_1!))
    edge_2 = Edge([:ALL], getfield(Main, :cc_aut_per_midlow_2), getfield(Main, :us_aut_per_midlow_2!))
    edge_3 = Edge([:ALL], getfield(Main, :cc_aut_per_midlow_3), getfield(Main, :us_aut_per_midlow_3!))
    edge_4 = Edge([:ALL], getfield(Main, :cc_aut_per_midlow_4), getfield(Main, :us_aut_per_midlow_4!))
    map_edges[:mid][:low] = [edge_1, edge_2, edge_3, edge_4]
    # * mid => high
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_midhigh_1), getfield(Main, :us_aut_per_midhigh_1!))
    map_edges[:mid][:high] = [edge_1]
    # * mid => final
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_midfinal_1), getfield(Main, :us_aut_per_midfinal_1!))
    map_edges[:mid][:final] = [edge_1]

    # high 
    # * high => high
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_highhigh_1), getfield(Main, :us_aut_per_highhigh_1!))
    map_edges[:high][:high] = [edge_1]
    # * high => mid
    edge_1 = Edge([:ALL], getfield(Main, :cc_aut_per_highmid_1), getfield(Main, :us_aut_per_highmid_1!))
    map_edges[:high][:mid] = [edge_1]
    # * high => final
    edge_1 = Edge([nothing], getfield(Main, :cc_aut_per_highfinal_1), getfield(Main, :us_aut_per_highfinal_1!))
    map_edges[:high][:final] = [edge_1] 

    ## Constants
    constants = Dict{Symbol,Float64}(:N => N, :L => L, :H => H, :initT => initT)

    A = LHA("Period", m.transitions, locations, Λ_F, locations_init, locations_final, 
            map_var_automaton_idx, flow, map_edges, constants, m.map_var_idx)
    return A
end

export create_period_automaton

