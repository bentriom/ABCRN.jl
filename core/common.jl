
abstract type Model end 
abstract type AbstractTrajectory end

const Transition = Union{String,Nothing}
const Location = String
const VariableAutomaton = String

mutable struct ContinuousTimeModel <: Model
    d::Int # state space dim
    k::Int # parameter space dim
    map_var_idx::Dict{String,Int} # maps str to full state space
    _map_obs_var_idx::Dict{String,Int} # maps str to observed state space
    map_param_idx::Dict{String,Int} # maps str in parameter space
    l_transitions::Vector{Transition}
    p::Vector{Float64}
    x0::Vector{Int}
    t0::Float64
    f!::Function
    g::Vector{String} # of dimension dobs
    _g_idx::Vector{Int} # of dimension dobs
    isabsorbing::Function
    time_bound::Float64
    buffer_size::Int
end

struct Trajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::Vector{Vector{Int}}
    times::Vector{Float64}
    transitions::Vector{Transition}
end

struct Edge
    transitions::Vector{Transition}
    check_constraints::Function
    update_state!::Function
end

struct LHA
    l_transitions::Vector{Transition}
    l_loc::Vector{Location} 
    Λ::Dict{Location,Function}
    l_loc_init::Vector{Location}
    l_loc_final::Vector{Location}
    map_var_automaton_idx::Dict{VariableAutomaton,Int} # nvar keys : str_var => idx in l_var
    l_flow::Dict{Location,Vector{Float64}} # output of length nvar
    map_edges::Dict{Tuple{Location,Location},Vector{Edge}}
    l_ctes::Dict{String,Float64}
    map_var_model_idx::Dict{String,Int} # of dim d (of a model)
end

mutable struct StateLHA
    A::LHA
    loc::Location
    l_var::Vector{Float64}
    time::Float64
end

mutable struct SynchronizedModel <: Model
    m::ContinuousTimeModel
    automaton::LHA
end

struct SynchronizedTrajectory <: AbstractTrajectory
    S::StateLHA
    m::SynchronizedModel
    values::Vector{Vector{Int}}
    times::Vector{Float64}
    transitions::Vector{Transition}
end

# Constructors
function ContinuousTimeModel(d::Int, k::Int, map_var_idx::Dict, map_param_idx::Dict, l_transitions::Vector{String}, 
              p::Vector{Float64}, x0::Vector{Int}, t0::Float64, 
              f!::Function, isabsorbing::Function; 
              g::Vector{String} = keys(map_var_idx), time_bound::Float64 = Inf, buffer_size::Int = 10)
    dobs = length(g)
    _map_obs_var_idx = Dict()
    _g_idx = Vector{Int}(undef, dobs)
    for i = 1:dobs
        _g_idx[i] = map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::String => idx in state space )
        _map_obs_var_idx[g[i]] = i
    end
  
    if length(methods(f!)) >= 2
        @warn "You have possibly redefined a function Model.f! used in a previously instantiated model."
    end
    if length(methods(isabsorbing)) >= 2
        @warn "You have possibly redefined a function Model.isabsorbing used in a previously instantiated model."
    end

    return ContinuousTimeModel(d, k, map_var_idx, _map_obs_var_idx, map_param_idx, l_transitions, p, x0, t0, f!, g, _g_idx, isabsorbing, time_bound, buffer_size)
end

LHA(A::LHA, map_var::Dict{String,Int}) = LHA(A.l_transitions, A.l_loc, A.Λ, 
                                             A.l_loc_init, A.l_loc_final, A.map_var_automaton_idx, A.l_flow,
                                             A.map_edges, A.l_ctes, map_var)
Base.:*(m::ContinuousTimeModel, A::LHA) = SynchronizedModel(m, A)
Base.:*(A::LHA, m::ContinuousTimeModel) = SynchronizedModel(m, A)

