
# File for benchmarking simulation and memory access of the package.

# Trajectories

_get_values_col(σ::AbstractTrajectory, var::String) = 
@view σ.values[(σ.m)._map_obs_var_idx[var],:] 
_get_values_row(σ::AbstractTrajectory, var::String) = 
@view σ.values[:,(σ.m)._map_obs_var_idx[var]] 

_get_state_col(σ::AbstractTrajectory, idx::Int) = 
@view σ.values[:,idx]
_get_state_row(σ::AbstractTrajectory, idx::Int) = 
@view σ.values[idx,:]

_get_value_col(σ::AbstractTrajectory, var::String, idx::Int) = 
σ.values[(σ.m)._map_obs_var_idx[var],idx] 
_get_value_row(σ::AbstractTrajectory, var::String, idx::Int) = 
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
    return Trajectory(m, values, times, transitions)
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
    return Trajectory(m, values, times, transitions)
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
    l_tr = Vector{String}(undef, buffer_size)
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
    return Trajectory(m, values, times, transitions)
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
    l_tr = Vector{String}(undef, buffer_size)
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
    return Trajectory(m, values, times, transitions)
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
    l_tr = Vector{String}(undef, m.buffer_size)
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
    return Trajectory(m, values, times, transitions)
end

