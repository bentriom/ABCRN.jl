
load_model(name_model::String) = include(get_module_path() * "/models/" * name_model * ".jl")

# Simulation
function simulate(m::ContinuousTimeModel)
    # trajectory fields
    full_values = Vector{Vector{Int}}(undef, m.d)
    for i = 1:m.d full_values[i] = Int[m.x0[i]] end
    for i = 1:m.d sizehint!(full_values[i], m.estim_min_states) end
    times = Float64[m.t0]
    transitions = Transition[nothing]
    # values at time n
    n = 0
    xn = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{String}(undef, m.buffer_size)
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    end_idx = -1
    # use sizehint! ?
    while !isabsorbing && tn <= m.time_bound
        for i = 1:m.buffer_size
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            tn = l_t[i]
            if tn > m.time_bound
                i -= 1 # 0 is an ok value, 1:0 is allowed
                end_idx = i
                break
            end
            xn = view(mat_x, i, :)
            isabsorbing = m.isabsorbing(m.p,xn)
            if isabsorbing 
                end_idx = i
                break
            end
        end
        if end_idx != -1
            break 
        end
        for k = 1:m.d
            append!(full_values[k], view(mat_x, :, k))
        end
        append!(times, l_t)
        append!(transitions,  l_tr)
        n += m.buffer_size
    end
    for k = 1:m.d
        append!(full_values[k], view(mat_x, 1:end_idx, k))
    end
    append!(times, view(l_t, 1:end_idx))
    append!(transitions,  view(l_tr, 1:end_idx))
    if isbounded(m)
        for k = 1:m.d
            push!(full_values[k], full_values[k][end])
        end
        push!(times, m.time_bound)
        push!(transitions, nothing)
    end
    values = full_values[m._g_idx]
    return Trajectory(m, values, times, transitions)
end

function simulate(product::SynchronizedModel)
    # trajectory fields
    m = product.m
    A = product.automaton
    full_values = Vector{Vector{Int}}(undef, m.d)
    for i = 1:m.d full_values[i] = Int[m.x0[i]] end
    for i = 1:m.d sizehint!(full_values[i], m.estim_min_states) end
    times = Float64[m.t0]
    transitions = Union{String,Nothing}[nothing]
    reshaped_x0 = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
    S0 = init_state(A, m.x0, m.t0)
    # values at time n
    n = 0
    xn = reshaped_x0 
    tn = m.t0
    Sn = copy(S0)
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{String}(undef, m.buffer_size)
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    isacceptedLHA::Bool = isaccepted(Sn)
    Snplus1 = copy(Sn)
    while !isabsorbing && tn <= m.time_bound && !isacceptedLHA
        i = 0
        while i < m.buffer_size && !isabsorbing && tn <= m.time_bound && !isacceptedLHA
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            tn = l_t[i]
            if tn > m.time_bound
                i -= 1 # 0 is an ok value, 1:0 is allowed
                break
            end
            xn = view(mat_x, i, :)
            tr_n = l_tr[i]
            next_state!(Snplus1, A, xn, tn, tr_n, Sn)
            Sn = Snplus1
            isabsorbing = m.isabsorbing(m.p,xn)
            isacceptedLHA = isaccepted(Snplus1)
        end
        for k = 1:m.d
            append!(full_values[k], view(mat_x, 1:i, k))
        end
        append!(times, view(l_t, 1:i))
        append!(transitions,  view(l_tr, 1:i))
        n += i
    end
    # When the trajectory is accepted, we should not add an end value
    if isbounded(m) && !isaccepted(Sn)
        @assert times[end] < m.time_bound
        for k = 1:m.d
            push!(full_values[k], full_values[k][end])
        end
        push!(times, m.time_bound)
        push!(transitions, nothing)
    end
    values = full_values[m._g_idx]
    return SynchronizedTrajectory(Snplus1, product, values, times, transitions)
end

function simulate(m::ContinuousTimeModel, n::Int)
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

isbounded(m::ContinuousTimeModel) = m.time_bound < Inf
function check_consistency(m::ContinuousTimeModel) 
    @assert m.d == length(m.map_var_idx) 
    @assert m.k == length(m.map_param_idx)
    @assert m.k == length(m.p)
    @assert length(m.g) <= m.d
    @assert length(m._g_idx) == length(m.g)
    @assert m.buffer_size >= 0
    @assert typeof(m.isabsorbing(m.p, view(reshape(m.x0, 1, m.d), 1, :))) == Bool
    return true
end

set_param!(m::ContinuousTimeModel, p::Vector{Float64}) = (m.p = p)
set_param!(m::ContinuousTimeModel, name_p::String, p_i::Float64) = (m.p[m.map_param_idx[name_p]] = p_i)
function set_param!(m::ContinuousTimeModel, l_name_p::Vector{String}, p::Vector{Float64}) 
    nb_param = length(l_name_p)
    for i = 1:nb_param
        set_param!(m, l_name_p[i], p[i])
    end
end

get_param(m::ContinuousTimeModel) = m.p
getindex(m::ContinuousTimeModel, name_p::String) = m.p[m.map_param_idx[name_p]]
set_time_bound!(m::ContinuousTimeModel, b::Float64) = (m.time_bound = b)
set_time_bound!(sm::SynchronizedModel, b::Float64) = set_time_bound!(sm.m, b)

function getproperty(m::ContinuousTimeModel, sym::Symbol)
    if sym == :dobs
        return length(m.g)
    else
        return getfield(m, sym)
    end
end
get_proba_model(m::ContinuousTimeModel) = m
get_proba_model(sm::SynchronizedModel) = sm.m

# Prior methods
function draw!(Π::ModelPrior)
    dict_dist = Π.map_l_param_dist
    for l_name in keys(dict_dist)
        if length(l_name) == 1
            set_param!(get_proba_model(Π.m), l_name[1], rand(dict_dist[l_name]))
        else
            set_param!(get_proba_model(Π.m), l_name, rand(dict_dist[l_name]))
        end
    end
end

