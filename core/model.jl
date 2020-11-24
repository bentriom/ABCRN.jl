
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
    xn = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
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
    # First we fill the allocated array
    mat_x = zeros(Int, 1, m.d)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    for i = 2:m.estim_min_states
        m.f!(mat_x, l_t, l_tr, 1, xn, tn, m.p)
        tn = l_t[1]
        if tn > m.time_bound
            break
        end
        n += 1
        xn = view(mat_x, 1, :)
        # Updating value
        for k = eachindex(full_values) full_values[k][n] = xn[k] end
        times[n] = tn
        transitions[n] = l_tr[1]
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
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{Transition}(undef, m.buffer_size)
    # Alloc buffer values
    tmp_full_values = Vector{Vector{Int}}(undef, m.d)
    for i = eachindex(tmp_full_values) tmp_full_values[i] = zeros(Int, 0) end
    tmp_times = zeros(Float64, 0)
    tmp_transitions = Vector{Transition}(undef, 0)
    tmp_idx = 0
    while !isabsorbing && tn <= m.time_bound
        _resize_trajectory!(tmp_full_values, tmp_times, tmp_transitions, tmp_idx+m.buffer_size)
        i = 0
        while i < m.buffer_size
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            tn = l_t[i]
            if tn > m.time_bound
                i -= 1
                break
            end
            xn = view(mat_x, i, :)
            isabsorbing = m.isabsorbing(m.p,xn)
            if isabsorbing 
                break
            end
        end
        # Update values
        rng_tmp = (tmp_idx+1):(tmp_idx+m.buffer_size)
        for k = eachindex(tmp_full_values) tmp_full_values[k][rng_tmp] = view(mat_x, :, k) end
        tmp_times[rng_tmp] = l_t
        tmp_transitions[rng_tmp] = l_tr
        if i < m.buffer_size
            _resize_trajectory!(tmp_full_values, tmp_times, tmp_transitions, tmp_idx+i)
        end
        tmp_idx += i
        n += i
    end
    # Push the temporary values 
    for k = eachindex(full_values) append!(full_values[k], tmp_full_values[k]) end
    append!(times, tmp_times)
    append!(transitions, tmp_transitions)
    
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
    xn = view(reshape(m.x0, 1, m.d), 1, :) # View for type stability
    tn = m.t0 
    Sn = copy(S0)
    # at time n+1
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    isacceptedLHA::Bool = isaccepted(Sn)
    if isabsorbing || isacceptedLHA
        _resize_trajectory!(full_values, times, transitions, 1)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return SynchronizedTrajectory(Sn, product, values, times, transitions)
    end
    # First we fill the allocated array
    mat_x = zeros(Int, 1, m.d)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    Snplus1 = copy(Sn)
    for i = 2:m.estim_min_states
        m.f!(mat_x, l_t, l_tr, 1, xn, tn, m.p)
        tn = l_t[1]
        if tn > m.time_bound
            break
        end
        n += 1
        xn = view(mat_x, 1, :)
        tr_n = l_tr[1]
        next_state!(Snplus1, A, xn, tn, tr_n, Sn)
        Sn = Snplus1
        # Updating value
        for k = eachindex(full_values) full_values[k][n] = xn[k] end
        times[n] = tn
        transitions[n] = l_tr[1]
        isabsorbing = m.isabsorbing(m.p,xn)
        isacceptedLHA = isaccepted(Snplus1)
        if isabsorbing || isacceptedLHA 
            break
        end
    end
    # If simulation ended
    if n < m.estim_min_states
        _resize_trajectory!(full_values, times, transitions, n)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        return SynchronizedTrajectory(Snplus1, product, values, times, transitions)
    end
    # Otherwise, buffering system
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{Transition}(undef, m.buffer_size)
    # Alloc buffer values
    tmp_full_values = Vector{Vector{Int}}(undef, m.d)
    for i = eachindex(tmp_full_values) tmp_full_values[i] = zeros(Int, 0) end
    tmp_times = zeros(Float64, 0)
    tmp_transitions = Vector{Transition}(undef, 0)
    tmp_idx = 0
    while !isabsorbing && tn <= m.time_bound && !isacceptedLHA
        _resize_trajectory!(tmp_full_values, tmp_times, tmp_transitions, tmp_idx+m.buffer_size)
        i = 0
        while i < m.buffer_size
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            tn = l_t[i]
            if tn > m.time_bound
                i -= 1
                break
            end
            xn = view(mat_x, i, :)
            tr_n = l_tr[i]
            next_state!(Snplus1, A, xn, tn, tr_n, Sn)
            Sn = Snplus1
            isabsorbing = m.isabsorbing(m.p,xn)
            isacceptedLHA = isaccepted(Snplus1)
            if isabsorbing || isacceptedLHA 
                break
            end
        end
        # Update values
        rng_tmp = (tmp_idx+1):(tmp_idx+m.buffer_size)
        for k = eachindex(tmp_full_values) tmp_full_values[k][rng_tmp] = view(mat_x, :, k) end
        tmp_times[rng_tmp] = l_t
        tmp_transitions[rng_tmp] = l_tr
        if i < m.buffer_size
            _resize_trajectory!(tmp_full_values, tmp_times, tmp_transitions, tmp_idx+i)
        end
        tmp_idx += i
        n += i
    end
    # Push the temporary values 
    for k = eachindex(full_values) append!(full_values[k], tmp_full_values[k]) end
    append!(times, tmp_times)
    append!(transitions, tmp_transitions)
    
    values = full_values[m._g_idx]
    if isbounded(m) && !isaccepted(Sn)
        # Add last value: the convention is the last transition is nothing,
        # the trajectory is bounded
        _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
    end
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

