
import Random: rand, rand!
import Distributions: insupport, pdf

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
function simulate(m::ContinuousTimeModel; p::Union{Nothing,AbstractVector{Float64}} = nothing)
    p_sim = m.p
    if p != nothing
        p_sim = p
    end
    # First alloc
    full_values = Vector{Vector{Int}}(undef, m.dim_state)
    for i = eachindex(full_values) full_values[i] = zeros(Int, m.estim_min_states) end
    times = zeros(Float64, m.estim_min_states)
    transitions = Vector{Transition}(undef, m.estim_min_states)
    # Initial values
    for i = eachindex(full_values) full_values[i][1] = m.x0[i] end
    times[1] = m.t0
    transitions[1] = nothing
    # Values at time n
    n = 1
    xn = copy(m.x0)
    tn = m.t0 
    # at time n+1
    isabsorbing::Bool = m.isabsorbing(p_sim,xn)
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
    vec_x = zeros(Int, m.dim_state)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    # First we fill the allocated array
    for i = 2:m.estim_min_states
        m.f!(vec_x, l_t, l_tr, xn, tn, p_sim)
        tn = l_t[1]
        if tn > m.time_bound || vec_x == xn
            isabsorbing = (vec_x == xn)
            break
        end
        n += 1
        copyto!(xn, vec_x)
        _update_values!(full_values, times, transitions, xn, tn, l_tr[1], i)
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
            m.f!(vec_x, l_t, l_tr, xn, tn, p_sim)
            tn = l_t[1]
            if tn > m.time_bound 
                i -= 1
                break
            end
            if vec_x == xn
                isabsorbing = true
                i -= 1
                break
            end
            copyto!(xn, vec_x)
            _update_values!(full_values, times, transitions, 
                            xn, tn, l_tr[1], m.estim_min_states+size_tmp+i)
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

function simulate(product::SynchronizedModel; p::Union{Nothing,AbstractVector{Float64}} = nothing,
                  verbose::Bool = false)
    m = product.m
    A = product.automaton
    p_sim = m.p
    if p != nothing
        p_sim = p
    end
    # First alloc
    full_values = Vector{Vector{Int}}(undef, m.dim_state)
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
    xn = copy(m.x0)
    tn = m.t0 
    Sn = copy(S0)
    isabsorbing::Bool = m.isabsorbing(p_sim,xn)
    isacceptedLHA::Bool = isaccepted(Sn)
    # Alloc of vectors where we stock n+1 values
    vec_x = zeros(Int, m.dim_state)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    Snplus1 = copy(Sn)
    # If x0 is absorbing
    if isabsorbing || isacceptedLHA 
        _resize_trajectory!(full_values, times, transitions, 1)
        values = full_values[m._g_idx]
        if isbounded(m)
            _finish_bounded_trajectory!(values, times, transitions, m.time_bound)
        end
        if isabsorbing && !isacceptedLHA
            next_state!(Snplus1, A, xn, m.time_bound, nothing, Sn; verbose = verbose)
            copyto!(Sn, Snplus1)
        end
        return SynchronizedTrajectory(Sn, product, values, times, transitions)
    end
    # First we fill the allocated array
    for i = 2:m.estim_min_states
        m.f!(vec_x, l_t, l_tr, xn, tn, p_sim)
        tn = l_t[1]
        if tn > m.time_bound || vec_x == xn
            isabsorbing = (vec_x == xn)
            break
        end
        n += 1
        copyto!(xn, vec_x)
        tr_n = l_tr[1]
        next_state!(Snplus1, A, xn, tn, tr_n, Sn; verbose = verbose)
        _update_values!(full_values, times, transitions, xn, tn, tr_n, i)
        copyto!(Sn, Snplus1)
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
        if isabsorbing && !isacceptedLHA
            next_state!(Snplus1, A, xn, m.time_bound, nothing, Sn; verbose = verbose)
            copyto!(Sn, Snplus1)
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
            m.f!(vec_x, l_t, l_tr, xn, tn, p_sim)
            tn = l_t[1]
            if tn > m.time_bound
                i -= 1
                break
            end
            if vec_x == xn
                isabsorbing = true
                i -= 1
                break
            end
            copyto!(xn, vec_x)
            tr_n = l_tr[1]
            next_state!(Snplus1, A, xn, tn, tr_n, Sn; verbose = verbose)
            _update_values!(full_values, times, transitions, 
                            xn, tn, tr_n, m.estim_min_states+size_tmp+i)
            copyto!(Sn, Snplus1)
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
    if isabsorbing && !isacceptedLHA
        next_state!(Snplus1, A, xn, m.time_bound, nothing, Sn; verbose = verbose)
        copyto!(Sn, Snplus1)
    end
    return SynchronizedTrajectory(Sn, product, values, times, transitions)
end
"""
    `volatile_simulate(sm::SynchronizedModel; p, verbose)`

Simulates a model synchronized with an automaton but does not store the values of the simulation
in order to improve performance.
It returns the last state of the simulation `S::StateLHA` not a trajectory `σ::SynchronizedTrajectory`.
"""
function volatile_simulate(product::SynchronizedModel; 
                           p::Union{Nothing,AbstractVector{Float64}} = nothing, verbose::Bool = false)
    m = product.m
    A = product.automaton
    p_sim = m.p
    if p != nothing
        p_sim = p
    end
    S0 = init_state(A, m.x0, m.t0)
    # Values at time n
    n = 1
    xn = copy(m.x0)
    tn = m.t0 
    Sn = copy(S0)
    isabsorbing::Bool = m.isabsorbing(p_sim,xn)
    isacceptedLHA::Bool = isaccepted(Sn)
    # Alloc of vectors where we stock n+1 values
    vec_x = zeros(Int, m.dim_state)
    l_t = Float64[0.0]
    l_tr = Transition[nothing]
    Snplus1 = copy(Sn)
    # If x0 is absorbing
    if isabsorbing || isacceptedLHA 
        if !isacceptedLHA
            next_state!(Snplus1, A, xn, m.time_bound, nothing, Sn; verbose = verbose)
            copyto!(Sn, Snplus1)
        end
        return Sn
    end
    while !isabsorbing && tn <= m.time_bound && !isacceptedLHA
        m.f!(vec_x, l_t, l_tr, xn, tn, p_sim)
        tn = l_t[1]
        if tn > m.time_bound
            i -= 1
            break
        end
        if vec_x == xn
            isabsorbing = true
            break
        end
        copyto!(xn, vec_x)
        tr_n = l_tr[1]
        next_state!(Snplus1, A, xn, tn, tr_n, Sn; verbose = verbose)
        copyto!(Sn, Snplus1)
        isacceptedLHA = isaccepted(Snplus1)
        n += 1
    end
    if isabsorbing && !isacceptedLHA
        next_state!(Snplus1, A, xn, m.time_bound, nothing, Sn; verbose = verbose)
        copyto!(Sn, Snplus1)
    end
    return Sn
end


"""
    `simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})

Simulates the model contained in pm with p_prior values.
It simulates the model by taking the parameters contained in get_proba_model(pm).p and
replace the 1D parameters pm.params with p_prior.
"""
function simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})
    full_p = copy(get_proba_model(pm).p)
    full_p[pm._param_idx] = p_prior
    
    return simulate(pm.m; p = full_p) 
end
"""
    `volatile_simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})

A volatile version of `simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})`.
The model in pm should be of type SynchronizedModel (`typeof(pm.m) <: SynchronizedModel`).
It returns `S::StateLHA`, not a trajectory.
"""
function volatile_simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})
    @assert typeof(pm.m) <: SynchronizedModel
    full_p = copy(get_proba_model(pm).p)
    full_p[pm._param_idx] = p_prior
    
    return volatile_simulate(pm.m; p = full_p) 
end
"""
    `distribute_mean_value_lha(sm::SynchronizedModel, str_var::String, nbr_stim::Int)`

Distribute over workers the computation of the mean value 
of a LHA over `nbr_sim` simulations of the model.
"""
function distribute_mean_value_lha(sm::SynchronizedModel, str_var::String, nbr_sim::Int)
    sum_val = @distributed (+) for i = 1:nbr_sim 
        volatile_simulate(sm)[str_var] 
    end
    return sum_val / nbr_sim
end

function mean_value_lha(sm::SynchronizedModel, str_var::String, nbr_sim::Int)
    sum_val = 0.0 
    for i = 1:nbr_sim 
        sum_val += volatile_simulate(sm)[str_var] 
    end
    return sum_val / nbr_sim
end

function distribute_prob_accept_lha(sm::SynchronizedModel, nbr_sim::Int)
    sum_val = @distributed (+) for i = 1:nbr_sim 
        Int(isaccepted(volatile_simulate(sm))) 
    end
    return sum_val / nbr_sim
end

isbounded(m::Model) = get_proba_model(m).time_bound < Inf
function check_consistency(m::ContinuousTimeModel) 
    @assert m.dim_state == length(m.map_var_idx) 
    @assert m.dim_state == length(m.x0) 
    @assert m.dim_params == length(m.map_param_idx)
    @assert m.dim_params == length(m.p)
    @assert length(m.g) <= m.dim_state
    @assert length(m._g_idx) == length(m.g)
    @assert m.buffer_size >= 0
    @assert typeof(m.isabsorbing(m.p, m.x0)) == Bool
    return true
end

# Set and get Model fields
function set_observed_var!(am::Model, g::Vector{String})
    m = get_proba_model(am)
    dim_obs_state = length(g)
    _map_obs_var_idx = Dict{String}{Int}()
    _g_idx = zeros(Int, dim_obs_state)
    for i = 1:dim_obs_state
        _g_idx[i] = m.map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::String => idx in state space )
        _map_obs_var_idx[g[i]] = i
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = _map_obs_var_idx
end
function observe_all!(am::Model)
    m = get_proba_model(am)
    g = Vector{String}(undef, m.dim_state)
    _g_idx = collect(1:m.dim_state)
    for var in keys(m.map_var_idx)
        g[m.map_var_idx[var]] = var
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = m.map_var_idx
end
function set_param!(am::Model, new_p::Vector{Float64})
    m = get_proba_model(am)
    @assert length(new_p) == m.dim_params
    m.p = new_p
end
function set_param!(am::Model, name_p::String, p_i::Float64) 
    m = get_proba_model(am)
    m.p[m.map_param_idx[name_p]] = p_i
end
function set_param!(am::Model, l_name_p::Vector{String}, p::Vector{Float64}) 
    m = get_proba_model(am)
    @assert length(l_name_p) == length(p)
    for i = eachindex(l_name_p)
        set_param!(m, l_name_p[i], p[i])
    end
end
function set_x0!(am::Model, new_x0::Vector{Int})
    m = get_proba_model(am)
    @assert length(new_x0) == m.dim_state
    m.x0 = new_x0
end
set_time_bound!(am::Model, b::Float64) = (get_proba_model(am).time_bound = b)


get_param(am::Model) = get_proba_model(am).p
function getindex(am::Model, name_p::String)
    m = get_proba_model(am)
    m.p[m.map_param_idx[name_p]]
end
function getproperty(m::ContinuousTimeModel, sym::Symbol)
    if sym == :dim_obs_state
        return length(m.g)
    else
        return getfield(m, sym)
    end
end
function getproperty(pm::ParametricModel, sym::Symbol)
    if sym == :df
        return length(pm.params)
    else
        return getfield(pm, sym)
    end
end

get_proba_model(m::ContinuousTimeModel) = m
get_proba_model(sm::SynchronizedModel) = sm.m
get_proba_model(pm::ParametricModel) = get_proba_model(pm.m)

get_observed_var(m::Model) = get_proba_model(am).g

# Prior methods
"""
    `draw_model!(pm::ParametricModel)`

Draw a parameter from the prior disitribution defined in `pm::ParametricModel`
and save it in the model contained in `pm`.
"""
draw_model!(pm::ParametricModel) = set_param!(get_proba_model(pm), pm.params, rand(pm.distribution))
"""
    `draw!(vec_p, pm)`

Draw a parameter from the prior distribution defined in pm and stores it in vec_p.
"""
draw!(vec_p::AbstractVector{Float64}, pm::ParametricModel) = rand!(pm.distribution, vec_p)
"""
    `draw!(mat_p, pm)`

Draw `size(mat_p)[2]` (number of columns of mat_p) parameters from the prior distribution 
defined in pm and stores it in mat_p.
"""
function draw!(mat_p::AbstractMatrix{Float64}, pm::ParametricModel)
    for i = 1:size(mat_p)[2]
        draw!(view(mat_p, :, i), pm)
    end
end
"""
    `prior_pdf(p_prior, pm)`
 
Computes the density of the prior distribution defined in pm on argument p_prior.
"""
prior_pdf(pm::ParametricModel, p_prior::AbstractVector{Float64}) = pdf(pm.distribution, p_prior)
"""
    `prior_pdf(vec_res, mat_p, pm)`
 
Computes the density of the prior distribution defined in pm on each column
ov mat_p. Stores it in mat_p. (`length(vec_res) == size(mat_p)[2]`)
"""
function prior_pdf!(res_pdf::AbstractVector{Float64}, pm::ParametricModel, mat_p::AbstractMatrix{Float64})
    for i = eachindex(res_pdf)
        res_pdf[i] = prior_pdf(pm, view(mat_p, :, i))
    end
end
fill!(pm::ParametricModel, p_prior::AbstractVector{Float64}) = get_proba_model(pm).p[pm._param_idx] = p_prior
insupport(pm::ParametricModel, p_prior::AbstractVector{Float64}) = insupport(pm.distribution, p_prior)

# to do: simulate(pm), create_res_dir, check_bounds, ajouter un champ _param_idx pour pm.
#

