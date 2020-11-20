
# Simulation
function simulate(m::ContinuousTimeModel)
    # trajectory fields
    full_values = Matrix{Int}(undef, 1, m.d)
    full_values[1,:] = m.x0
    times = Float64[m.t0]
    transitions = Union{String,Nothing}[nothing]
    # values at time n
    n = 0
    xn = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{String}(undef, m.buffer_size)
    is_absorbing::Bool = m.is_absorbing(m.p,xn)
    while !is_absorbing && (tn <= m.time_bound)
        i = 0
        while i < m.buffer_size && !is_absorbing && (tn <= m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            xn = view(mat_x, i, :)
            tn = l_t[i]
            is_absorbing = m.is_absorbing(m.p,xn)
        end
        full_values = vcat(full_values, view(mat_x, 1:i, :))
        append!(times, view(l_t, 1:i))
        append!(transitions,  view(l_tr, 1:i))
        n += i
        is_absorbing = m.is_absorbing(m.p,xn)
    end
    if is_bounded(m)
        if times[end] > m.time_bound
            full_values[end,:] = full_values[end-1,:]
            times[end] = m.time_bound
            transitions[end] = nothing
        else
            full_values = vcat(full_values, reshape(full_values[end,:], 1, m.d))
            push!(times, m.time_bound)
            push!(transitions, nothing)
        end
    end
    values = view(full_values, :, m._g_idx)
    return Trajectory(m, values, times, transitions)
end

function simulate(product::SynchronizedModel)
    # trajectory fields
    m = product.m
    A = product.automaton
    full_values = Matrix{Int}(undef, 1, m.d)
    full_values[1,:] = m.x0
    times = Float64[m.t0]
    transitions = Union{String,Nothing}[nothing]
    reshaped_x0 = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
    S0 = init_state(A, reshaped_x0, t0) 
    # values at time n
    n = 0
    xn = reshaped_x0 
    tn = m.t0
    Sn = copy(S0)
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{String}(undef, m.buffer_size)
    is_absorbing::Bool = m.is_absorbing(m.p,xn)
    Snplus1 = copy(Sn)
    while !is_absorbing && (tn < m.time_bound)
        i = 0
        while i < m.buffer_size && !is_absorbing && (tn < m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            xn = view(mat_x, i, :)
            tn = l_t[i]
            tr_n = l_tr[n]
            A.next_state!(Snplus1, A, xn, tn, tr_n, Sn)
            Sn = Snplus1
            is_absorbing = m.is_absorbing(m.p,xn)
        end
        full_values = vcat(full_values, view(mat_x, 1:i, :))
        append!(times, view(l_t, 1:i))
        append!(transitions,  view(l_tr, 1:i))
        n += i
        is_absorbing = m.is_absorbing(m.p,xn)
    end
    if is_bounded(m)
        if times[end] > m.time_bound
            full_values[end,:] = full_values[end-1,:]
            times[end] = m.time_bound
            transitions[end] = nothing
        else
            full_values = vcat(full_values, reshape(full_values[end,:], 1, m.d))
            push!(times, m.time_bound)
            push!(transitions, nothing)
        end
    end
    values = view(full_values, :, m._g_idx)
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

