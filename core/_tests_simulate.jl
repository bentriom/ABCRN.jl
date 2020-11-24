
struct OldTrajectory <: AbstractTrajectory 
    m::ContinuousTimeModel
    values::Matrix{Int}
    times::Vector{Float64}
    transitions::Vector{Union{String,Nothing}}
end 

# File for benchmarking simulation and memory access of the package.

# Trajectories

_get_values_col(σ::OldTrajectory, var::String) = 
@view σ.values[(σ.m)._map_obs_var_idx[var],:] 
_get_values_row(σ::OldTrajectory, var::String) = 
@view σ.values[:,(σ.m)._map_obs_var_idx[var]] 

_get_state_col(σ::OldTrajectory, idx::Int) = 
@view σ.values[:,idx]
_get_state_row(σ::OldTrajectory, idx::Int) = 
@view σ.values[idx,:]

_get_value_col(σ::OldTrajectory, var::String, idx::Int) = 
σ.values[(σ.m)._map_obs_var_idx[var],idx] 
_get_value_row(σ::OldTrajectory, var::String, idx::Int) = 
σ.values[idx,(σ.m)._map_obs_var_idx[var]] 

# Model

function _simulate_col(m::ContinuousTimeModel)
    # trajectory fields
    full_values = zeros(m.d, 0)
    times = zeros(0)
    transitions = Vector{Union{String,Nothing}}(undef,0)
    # values at time n
    n = 0
    xn = @view m.x0[:]
    tn = m.t0 
    tr = [""]
    # at time n+1
    xnplus1 = zeros(Int, m.d)
    tnplus1 = zeros(Float64, 1)
    isabsorbing = (m.isabsorbing(m.p,xn))::Bool
    while !isabsorbing && (tn <= m.time_bound)
        m.f!(xnplus1, tnplus1, tr, xn, tn, m.p)
        full_values = hcat(full_values, xnplus1)
        push!(times, tnplus1[1])
        push!(transitions, tr[1])
        xn, tn = xnplus1, tnplus1[1]
        n += 1
        isabsorbing = m.isabsorbing(m.p,xn)::Bool
    end
    values = @view full_values[m._g_idx,:] 
    if isbounded(m)
        if times[end] > m.time_bound
            values[:, end] = values[:,end-1]
            times[end] = m.time_bound
            transitions[end] = nothing
        end
    end
    return OldTrajectory(m, values, times, transitions)
end

function _simulate_row(m::ContinuousTimeModel)
    # trajectory fields
    full_values = zeros(m.d, 0)
    times = zeros(0)
    transitions = Vector{Union{String,Nothing}}(undef,0)
    # values at time n
    n = 0
    xn = @view m.x0[:]
    tn = m.t0 
    tr = [""]
    # at time n+1
    xnplus1 = zeros(Int, m.d)
    tnplus1 = zeros(Float64, 1)
    isabsorbing = (m.isabsorbing(m.p,xn))::Bool
    while !isabsorbing && (tn <= m.time_bound)
        m.f!(xnplus1, tnplus1, tr, xn, tn, m.p)
        full_values = vcat(full_values, xnplus1)
        push!(times, tnplus1[1])
        push!(transitions, tr[1])
        xn, tn = xnplus1, tnplus1[1]
        n += 1
        isabsorbing = m.isabsorbing(m.p,xn)::Bool
    end
    values = @view full_values[m._g_idx,:] 
    if isbounded(m)
        if times[end] > m.time_bound
            values[:, end] = values[:,end-1]
            times[end] = m.time_bound
            transitions[end] = nothing
        end
    end
    return OldTrajectory(m, values, times, transitions)
end


function _simulate_col_buffer(m::ContinuousTimeModel; buffer_size::Int = 5)
    # trajectory fields
    full_values = zeros(m.d, 0)
    times = zeros(0)
    transitions = Vector{Union{String,Nothing}}(undef,0)
    # values at time n
    n = 0
    xn = @view m.x0[:]
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, m.d, buffer_size)
    l_t = zeros(Float64, buffer_size)
    l_tr = Vector{Union{String,Nothing}}(undef, buffer_size)
    isabsorbing = m.isabsorbing(m.p,xn)::Bool
    while !isabsorbing && (tn <= m.time_bound)
        i = 0
        while i < buffer_size && !isabsorbing && (tn <= m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            xn = @view mat_x[:,i]
            tn = l_t[i]
            isabsorbing = m.isabsorbing(m.p,xn)::Bool
        end
        full_values = hcat(full_values, @view mat_x[:,1:i])
        append!(times, @view l_t[1:i])
        append!(transitions,  @view l_tr[1:i])
        n += i
        isabsorbing = m.isabsorbing(m.p,xn)::Bool
    end
    values = @view full_values[m._g_idx,:] 
    if isbounded(m)
        if times[end] > m.time_bound
            values[:, end] = values[:,end-1]
            times[end] = m.time_bound
            transitions[end] = nothing
        end
    end
    return OldTrajectory(m, values, times, transitions)
end

function _simulate_row_buffer(m::ContinuousTimeModel; buffer_size::Int = 5)
    # trajectory fields
    full_values = zeros(0, m.d)
    times = zeros(0)
    transitions = Vector{Union{String,Nothing}}(undef,0)
    # values at time n
    n = 0
    xn = @view m.x0[:]
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, buffer_size, m.d)
    l_t = zeros(Float64, buffer_size)
    l_tr = Vector{Union{String,Nothing}}(undef, buffer_size)
    isabsorbing = m.isabsorbing(m.p,xn)::Bool
    while !isabsorbing && (tn <= m.time_bound)
        i = 0
        while i < buffer_size && !isabsorbing && (tn <= m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            xn = @view mat_x[i,:]
            tn = l_t[i]
            isabsorbing = m.isabsorbing(m.p,xn)::Bool
        end
        full_values = vcat(full_values, @view mat_x[1:i,:])
        append!(times, @view l_t[1:i])
        append!(transitions,  @view l_tr[1:i])
        n += i
        isabsorbing = m.isabsorbing(m.p,xn)::Bool
    end
    values = @view full_values[:,m._g_idx] 
    if isbounded(m)
        if times[end] > m.time_bound
            values[end,:] = values[end-1,:]
            times[end] = m.time_bound
            transitions[end] = nothing
        end
    end
    return OldTrajectory(m, values, times, transitions)
end

function _simulate_without_view(m::ContinuousTimeModel)
    # trajectory fields
    full_values = Matrix{Int}(undef, 1, m.d)
    full_values[1,:] = m.x0
    times = Float64[m.t0]
    transitions = Union{String,Nothing}[nothing]
    # values at time n
    n = 0
    xn = @view m.x0[:]
    tn = m.t0 
    # at time n+1
    mat_x = zeros(Int, m.buffer_size, m.d)
    l_t = zeros(Float64, m.buffer_size)
    l_tr = Vector{Union{String,Nothing}}(undef, m.buffer_size)
    isabsorbing = m.isabsorbing(m.p,xn)::Bool
    while !isabsorbing && (tn <= m.time_bound)
        i = 0
        while i < m.buffer_size && !isabsorbing && (tn <= m.time_bound)
            i += 1
            m.f!(mat_x, l_t, l_tr, i, @view(xn[:]), tn, m.p)
            xn = mat_x[i,:]
            tn = l_t[i]
            isabsorbing = m.isabsorbing(m.p,@view(xn[:]))::Bool
        end
        full_values = vcat(full_values, mat_x[1:i,:])
        append!(times, l_t[1:i])
        append!(transitions,  l_tr[1:i])
        n += i
        isabsorbing = m.isabsorbing(m.p,@view(xn[:]))::Bool
    end
    if isbounded(m)
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
    values = full_values[:,m._g_idx] 
    return OldTrajectory(m, values, times, transitions)
end

# With trajectory values in Matrix type
function _simulate_27d56(m::ContinuousTimeModel)
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
    l_tr = Vector{Union{String,Nothing}}(undef, m.buffer_size)
    isabsorbing::Bool = m.isabsorbing(m.p,xn)
    # use sizehint! ?
    while !isabsorbing && tn <= m.time_bound
        i = 0
        while i < m.buffer_size && !isabsorbing && tn <= m.time_bound
            i += 1
            m.f!(mat_x, l_t, l_tr, i, xn, tn, m.p)
            tn = l_t[i]
            if tn > m.time_bound
                i -= 1 # 0 is an ok value, 1:0 is allowed
                break
            end
            xn = view(mat_x, i, :)
            isabsorbing = m.isabsorbing(m.p,xn)
        end
        full_values = vcat(full_values, view(mat_x, 1:i, :))
        append!(times, view(l_t, 1:i))
        append!(transitions,  view(l_tr, 1:i))
        n += i
    end
    if isbounded(m)
        #=
        if times[end] > m.time_bound
            full_values[end,:] = full_values[end-1,:]
            times[end] = m.time_bound
            transitions[end] = nothing
        else
        end
        =#
        full_values = vcat(full_values, reshape(full_values[end,:], 1, m.d))
        push!(times, m.time_bound)
        push!(transitions, nothing)
    end
    values = view(full_values, :, m._g_idx)
    return OldTrajectory(m, values, times, transitions)
end


function _simulate_d7458(m::ContinuousTimeModel)
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
    l_tr = Vector{Union{String,Nothing}}(undef, m.buffer_size)
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



function _simulate_last(product::SynchronizedModel)
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


