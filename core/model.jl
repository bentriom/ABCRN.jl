
load_model(name_model::String) = include(get_module_path() * "/models/" * name_model * ".jl")

function _resize_trajectory!(values::Vector{Vector{Int}}, times::Vector{Float64}, 
                             transitions::Vector{Transition}, size::Int)
    for i = eachindex(values) resize!(values[i], size) end
    resize!(times, size)
    resize!(transitions, size)
end


function _finish_bounded_trajectory!(values::Vector{Vector{Int}}, times::Vector{Float64}, 
                                    transitions::Vector{Transition}, time_bound::Float64)
    for i = eachindex(values) push!(values[i], values[i][end]) end
    push!(times, time_bound)
    push!(transitions, nothing)
end

function _update_values!(values::Vector{Vector{Int}}, times::Vector{Float64}, transitions::Vector{Transition},
                         xn::Vector{Int}, tn::Float64, tr_n::Transition, idx::Int)
    for k = eachindex(values) values[k][idx] = xn[k] end
    times[idx] = tn
    transitions[idx] = tr_n
end

"""
    `simulate(m)`

Simulates a model. If `m::SynchronizedModel`, the simulation is synchronized with a 
Linear Hybrid Automaton.
"""
function simulate(m::ContinuousTimeModel)
    # First alloc
    full_values = Vector{Vector{Int}}(undef, m.d)
    for i = eachindex(full_values) full_values[i] = zeros(Int, m.estim_min_states) end
    times = zeros(Float64, m.estim_min_states)
    transitions = Vector{Transition}(undef, m.estim_min_states)
    # Initial values
    for i = eachindex(full_values) full_values[i][1] = m.x0[i] end
    times[1] = m.t0
    transitions[1] = nothing
    # Values at time n
    n = 1
    xn = m.x0 # View for type stability
    tn = m.t0 
    # at time n+1
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    # If x0 is absorbing
    if isabsorbing
        _resize_trajectory!(full_values, times, transitions, 1)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return Trajectory(m, values, times, transitions)
    end
    # Alloc of vectors where we stock n+1 values
    vec_x = zeros(Int, m.d)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    # First we fill the allocated array
    for i = 2:m.estim_min_states
        m.f!(vec_x, l_t, l_tr, xn, tn, m.p)
        tn = l_t[1]
        if tn > m.time_bound
            break
        end
        n += 1
        xn = vec_x
        _update_values!(full_values, times, transitions, xn, tn, l_tr[1], i)
        isabsorbing = m.isabsorbing(m.p,xn)
        if isabsorbing 
            break
        end
    end
    # If simulation ended before the estimation of states
    if n < m.estim_min_states
        _resize_trajectory!(full_values, times, transitions, n)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return Trajectory(m, values, times, transitions)
    end
    # Otherwise, buffering system
    size_tmp = 0
    while !isabsorbing && tn <= m.time_bound
        # Alloc buffer
        _resize_trajectory!(full_values, times, transitions, m.estim_min_states+size_tmp+m.buffer_size)
        i = 0
        while i < m.buffer_size
            i += 1
            m.f!(vec_x, l_t, l_tr, xn, tn, m.p)
            tn = l_t[1]
            if tn > m.time_bound
                i -= 1
                break
            end
            xn = vec_x
            _update_values!(full_values, times, transitions, 
                            xn, tn, l_tr[1], m.estim_min_states+size_tmp+i)
            isabsorbing = m.isabsorbing(m.p,xn)
            if isabsorbing 
                break
            end
        end
        # If simulation ended before the end of buffer
        if i < m.buffer_size
            _resize_trajectory!(full_values, times, transitions, m.estim_min_states+size_tmp+i)
        end
        size_tmp += i
        n += i
    end
    values = full_values[m._g_idx]
    if isbounded(m)
        # Add last value: the convention is the last transition is nothing,
        # the trajectory is bounded
        _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
    end
    return Trajectory(m, values, times, transitions)
end
function simulate(product::SynchronizedModel)
    m = product.m
    A = product.automaton
    # First alloc
    full_values = Vector{Vector{Int}}(undef, m.d)
    for i = eachindex(full_values) full_values[i] = zeros(Int, m.estim_min_states) end
    times = zeros(Float64, m.estim_min_states)
    transitions = Vector{Transition}(undef, m.estim_min_states)
    # Initial values
    for i = eachindex(full_values) full_values[i][1] = m.x0[i] end
    times[1] = m.t0
    transitions[1] = nothing
    S0 = init_state(A, m.x0, m.t0)
    # Values at time n
    n = 1
    xn = m.x0 # View for type stability
    tn = m.t0 
    Sn = copy(S0)
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    isacceptedLHA::Bool = isaccepted(Sn)
    # If x0 is absorbing
    if isabsorbing || isacceptedLHA 
        _resize_trajectory!(full_values, times, transitions, 1)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return SynchronizedTrajectory(Sn, product, values, times, transitions)
    end
    # Alloc of vectors where we stock n+1 values
    vec_x = zeros(Int, m.d)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    Snplus1 = copy(Sn)
    # First we fill the allocated array
    for i = 2:m.estim_min_states
        m.f!(vec_x, l_t, l_tr, xn, tn, m.p)
        tn = l_t[1]
        if tn > m.time_bound
            break
        end
        n += 1
        xn = vec_x
        tr_n = l_tr[1]
        next_state!(Snplus1, A, xn, tn, tr_n, Sn)
        _update_values!(full_values, times, transitions, xn, tn, tr_n, i)
        Sn = Snplus1
        isabsorbing = m.isabsorbing(m.p,xn)
        isacceptedLHA = isaccepted(Snplus1)
        if isabsorbing || isacceptedLHA 
            break
        end
    end
    # If simulation ended before the estimation of states
    if n < m.estim_min_states
        _resize_trajectory!(full_values, times, transitions, n)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return SynchronizedTrajectory(Sn, product, values, times, transitions)
    end
    # Otherwise, buffering system
    size_tmp = 0
    while !isabsorbing && tn <= m.time_bound && !isacceptedLHA
        # Alloc buffer
        _resize_trajectory!(full_values, times, transitions, m.estim_min_states+size_tmp+m.buffer_size)
        i = 0
        while i < m.buffer_size
            i += 1
            m.f!(vec_x, l_t, l_tr, xn, tn, m.p)
            tn = l_t[1]
            if tn > m.time_bound
                i -= 1
                break
            end
            xn = vec_x
            tr_n = l_tr[1]
            next_state!(Snplus1, A, xn, tn, tr_n, Sn)
            _update_values!(full_values, times, transitions, 
                            xn, tn, tr_n, m.estim_min_states+size_tmp+i)
            Sn = Snplus1
            isabsorbing = m.isabsorbing(m.p,xn)
            isacceptedLHA = isaccepted(Snplus1)
            if isabsorbing || isacceptedLHA
                break
            end
        end
        # If simulation ended before the end of buffer
        if i < m.buffer_size
            _resize_trajectory!(full_values, times, transitions, m.estim_min_states+size_tmp+i)
        end
        size_tmp += i
        n += i
    end
    values = full_values[m._g_idx]
    if isbounded(m) && !isaccepted(Sn)
        # Add last value: the convention is the last transition is nothing,
        # the trajectory is bounded
        _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
    end
    return SynchronizedTrajectory(Sn, product, values, times, transitions)
end


function simulate(m::ContinuousTimeModel, n::Int)
end

function set_observed_var!(m::Model, g::Vector{String})
    dobs = length(g)
    _map_obs_var_idx = Dict{String}{Int}()
    _g_idx = zeros(Int, dobs)
    for i = 1:dobs
        _g_idx[i] = m.map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::String => idx in state space )
        _map_obs_var_idx[g[i]] = i
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = _map_obs_var_idx
end

function observe_all!(m::Model)
    g = Vector{String}(undef, m.d)
    _g_idx = collect(1:m.d)
    for var in keys(m.map_var_idx)
        g[m.map_var_idx[var]] = var
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = m.map_var_idx
end

isbounded(m::ContinuousTimeModel) = m.time_bound < Inf
function check_consistency(m::ContinuousTimeModel) 
    @assert m.d == length(m.map_var_idx) 
    @assert m.k == length(m.map_param_idx)
    @assert m.k == length(m.p)
    @assert length(m.g) <= m.d
    @assert length(m._g_idx) == length(m.g)
    @assert m.buffer_size >= 0
    @assert typeof(m.isabsorbing(m.p, m.x0)) == Bool
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
get_observed_var(m::ContinuousTimeModel) = m.g
get_observed_var(sm::SynchronizedModel) = (sm.m).g

# Prior methods
function draw_model!(pm::ParametricModel)
    dict_dist = pm.map_l_param_dist
    for l_name in keys(dict_dist)
        if length(l_name) == 1
            set_param!(get_proba_model(pm.m), l_name[1], rand(dict_dist[l_name]))
        else
            set_param!(get_proba_model(pm.m), l_name, rand(dict_dist[l_name]))
        end
    end
end

function dim_free_param(pm::ParametricModel)
    return 1
end

function draw!(vec_p::Vector{Float64}, pm::ParametricModel)
end

function draw!(mat_p::Matrix{Float64}, pm::ParametricModel)
end

function prior_density!(wl::Vector{Float64}, mat_p::Matrix{Float64}, pm::ParametricModel)
end

