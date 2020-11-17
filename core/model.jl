
import StaticArrays: SVector

abstract type Model end 
abstract type ContinuousTimeModel <: Model end 
abstract type DiscreteTimeModel <: Model end 

mutable struct CTMC <: ContinuousTimeModel
    d::Int # state space dim
    k::Int # parameter space dim
    map_var_idx::Dict # maps str to full state space
    _map_obs_var_idx::Dict # maps str to observed state space
    map_param_idx::Dict # maps str in parameter space
    l_name_transitions::Vector{String}
    p::Vector{Float64}
    x0::Vector{Int}
    t0::Float64
    f!::Function
    g::Vector{String} # of dimension dobs
    _g_idx::Vector{Int} # of dimension dobs
    is_absorbing::Function
    time_bound::Float64
    buffer_size::Int
end

function CTMC(d::Int, k::Int, map_var_idx::Dict, map_param_idx::Dict, l_name_transitions::Vector{String}, 
              p::Vector{Float64}, x0::Vector{Int}, t0::Float64, 
              f!::Function, is_absorbing::Function; 
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
    if length(methods(is_absorbing)) >= 2
        @warn "You have possibly redefined a function Model.is_absorbing used in a previously instantiated model."
    end

    return CTMC(d, k, map_var_idx, _map_obs_var_idx, map_param_idx, l_name_transitions, p, x0, t0, f!, g, _g_idx, is_absorbing, time_bound, buffer_size)
end

function simulate(m::ContinuousTimeModel)
    # trajectory fields
    full_values = zeros(1, m.d)
    full_values[1,:] = m.x0
    times = Float64[m.t0]
    transitions = Union{String,Nothing}[nothing]
    # values at time n
    n = 0
    xn = m.x0
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{String}(undef, m.buffer_size)
    is_absorbing = m.is_absorbing(m.p,xn)::Bool
    while !is_absorbing && (tn <= m.time_bound)
        i = 0
        while i < m.buffer_size && !is_absorbing && (tn <= m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            xn = @view mat_x[i,:]
            tn = l_t[i]
            is_absorbing = m.is_absorbing(m.p,xn)::Bool
        end
        full_values = vcat(full_values, @view mat_x[1:i,:])
        append!(times, @view l_t[1:i])
        append!(transitions,  @view l_tr[1:i])
        n += i
        is_absorbing = m.is_absorbing(m.p,xn)::Bool
    end
    values = @view full_values[:,m._g_idx] 
    if is_bounded(m)
        if times[end] > m.time_bound
            values[end,:] = values[end-1,:]
            times[end] = m.time_bound
            transitions[end] = nothing
        else
            vcat(values, values[end,:])
            push!(times, m.time_bound)
            push!(transitions, nothing)
        end
    end
    return Trajectory(m, values, times, transitions)
end

function simulate(m::ContinuousTimeModel, n::Int)
    obs = ContinuousObservations(undef, n)
    for i = 1:n
        obs[i] = simulate(m)
    end
    return obs
end

function set_observed_var!(m::Model,g::Vector{String})
    dobs = length(g)
    _map_obs_var_idx = Dict()
    _g_idx = Vector{Int}(undef, dobs)
    for i = 1:dobs
        _g_idx[i] = m.map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::String => idx in state space )
        _map_obs_var_idx[g[i]] = i
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = _map_obs_var_idx
end

is_bounded(m::Model) = m.time_bound < Inf
function check_consistency(m::Model) end
function simulate(m::Model, n::Int; bound::Float64 = Inf)::AbstractObservations end
function set_param!(m::Model, p::Vector{Float64})::Nothing end
function get_param(m::Model)::Vector{Float64} end

get_module_path() = dirname(dirname(pathof(@__MODULE__)))
function load_model(name_model::String)
    include(get_module_path() * "/models/" * name_model * ".jl")
end

