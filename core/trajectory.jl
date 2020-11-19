
abstract type AbstractTrajectory end
ContinuousObservations = AbstractVector{AbstractTrajectory}

struct Trajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::Matrix{Int}
    times::Vector{Float64}
    transitions::Vector{Union{String,Nothing}}
end

# Operations
function +(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function -(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end

# Top-level Lp distance function
"""
    `dist_lp(l_σ1, l_σ2; verbose, p, str_stat_list, str_stat_trajectory)`   

Function that computes Lp distance between two set of any dimensional trajectories.
...
# Arguments
- `l_σ1::Vector{AbstractTrajectory}` is the first set of trajectories. l_σ2 is the second.
- `verbose::Bool` If `true`, launch a verbose execution of the computation. 
- `str_stat_list::String` allows to set the statistic to apply on the distances of the trajectories. 
It is salso called linkage function in clustering field. Only "mean" is available for now on.
...
"""
function dist_lp(l_σ1::Vector{<:AbstractTrajectory}, l_σ2::Vector{<:AbstractTrajectory};  
                 verbose::Bool = false, p::Int = 1, str_stat_list::String = "mean", str_stat_trajectory::String = "mean")
    if str_stat_list == "mean"
        return dist_lp_mean(l_σ1, l_σ2; verbose = verbose, p = p, str_stat_trajectory = str_stat_trajectory)
    end
    error("Unrecognized statistic in dist_lp")
end

"""
    `dist_lp(σ1, σ2; verbose, p, str_stat)`   

Function that computes Lp distance between two trajectories of any dimension.
It computes Lp distances for each observed variable (contained in `get_obs_var(σ)`). and then apply a statistic on these distances.
Requires `get_obs_var(σ1) == get_obs_var(σ2)`, it is verified if they are simulated from the same model.
...
# Arguments
- `σ1::AbstractTrajectory` is the first trajectory. σ2 is the second.
- `verbose::Bool` If `true`, launch a verbose execution of the computation. 
- `str_stat::String` allows to set the statistic to apply on the distances computed  of the trajectories. 
Only "mean" is available for now on.
...
"""
function dist_lp(σ1::AbstractTrajectory, σ2::AbstractTrajectory;  
                 verbose::Bool = false, p::Int = 1, str_stat::String = "mean")
    if str_stat == "mean"
        return dist_lp_mean(σ1, σ2; verbose = verbose, p = p)
    end
    error("Unrecognized statistic in dist_lp")
end

"""
    `dist_lp(σ1, σ2, var; verbose, p, str_stat)`   

Function that computes Lp distance between two trajectories of any dimension.
It computes Lp distances for each observed variable (contained in `get_obs_var(σ)`). and then apply a statistic on these distances.
Requires `get_obs_var(σ1) == get_obs_var(σ2)`, it is verified if they are simulated from the same model.
...
# Arguments
- `σ1::AbstractTrajectory` is the first trajectory. σ2 is the second.
- `var::String` is an observed variable. Have to be contained in `get_obs_var(σ1)` and `get_obs_var(σ2)`.
- `verbose::Bool` If `true`, launch a verbose execution of the computation. 
...
"""
function dist_lp(σ1::AbstractTrajectory, σ2::AbstractTrajectory, var::String; 
                 verbose::Bool = false, p::Int = 1)
    if !is_bounded(σ1) || !is_bounded(σ2)
        @warn "Lp distance computed on unbounded trajectories. Result should be wrong"
    end
    return dist_lp(σ1[var], times(σ1), σ2[var], times(σ2); verbose = false, p = p)
end

# Distance function. Vectorized version
function dist_lp(x_obs::SubArray{Int,1}, t_x::Vector{Float64}, y_obs::SubArray{Int,1}, t_y::Vector{Float64}; 
                 verbose::Bool = false, p::Int = 1)
    current_y_obs = y_obs[1]
    current_t_y = t_y[2]
    idx = 1
    res = 0.0
    for i = 1:(length(x_obs)-1)
        if verbose
            @show i
            @show (t_x[i], x_obs[i])
            @show current_t_y
            @show t_x[i+1]
        end
        last_t_y = t_x[i]
        while current_t_y < t_x[i+1]
            rect =  abs(current_y_obs - x_obs[i])^p * (current_t_y - last_t_y)
            res += rect
            if verbose
                println("-- in loop :")
                println("-- add : $rect abs($current_y_obs - $(x_obs[i]))^p * ($current_t_y - $(last_t_y)) / $(abs(current_y_obs - x_obs[i])^p) * $(current_t_y - last_t_y)")
                print("-- ")
                @show current_y_obs, current_t_y
            end
            idx += 1
            last_t_y = t_y[idx]
            current_y_obs = y_obs[idx]
            current_t_y = t_y[idx+1]
        end
        last_t_read = max(last_t_y, t_x[i])
        rect = abs(current_y_obs - x_obs[i])^p * (t_x[i+1] - last_t_read)
        res += rect
        if verbose
            println("add : $rect abs($current_y_obs - $(x_obs[i]))^p * ($(t_x[i+1]) - $last_t_read) / $(abs(current_y_obs - x_obs[i])^p) * $(t_x[i+1] - last_t_read)")
            @show t_x[i+1], t_y[idx+1]
            @show y_obs[idx], x_obs[i]
            @show last_t_read
            @show current_y_obs
        end
    end
    return res^(1/p)
end
# For all the observed variables
function dist_lp_mean(σ1::AbstractTrajectory, σ2::AbstractTrajectory;  
                      verbose::Bool = false, p::Int = 1)
    if get_obs_var(σ1) != get_obs_var(σ2) error("Lp distances should be computed with the same observed variables") end
    res = 0.0
    for var in get_obs_var(σ1)
        res += dist_lp(σ1, σ2, var; verbose = verbose, p = p)
    end
    return res / length_obs_var(σ1)
end
# For a list of trajectories
function dist_lp_mean(l_σ1::Vector{<:AbstractTrajectory}, l_σ2::Vector{<:AbstractTrajectory};  
                 verbose::Bool = false, p::Int = 1, str_stat_trajectory::String = "mean")
    res = 0.0
    for σ1 in l_σ1
        for σ2 in l_σ2
            res += dist_lp(σ1, σ2; verbose = verbose, p = p, str_stat = str_stat_trajectory)
        end
    end
    return res / (length(l_σ1)*length(l_σ2))
end
# Get the value at any time
function _f_step(l_x::AbstractArray, l_t::AbstractArray, t::Float64)
    @assert length(l_x) == length(l_t)
    for i = 1:(length(l_t)-1)
        if l_t[i] <= t < l_t[i+1]
            return l_x[i]
        end
    end
    if l_t[end] == t
        return l_x[end]
    end
    return missing
end
# Riemann sum
function _riemann_sum(f::Function, t_begin::Real, t_end::Real, step::Float64)
    res = 0.0
    l_t = collect(t_begin:step:t_end)
    for i in 2:length(l_t)
        res += f(l_t[i]) * (l_t[i] - l_t[i-1])
    end
    return res
end

function check_consistency(σ::AbstractTrajectory)
    test_sizes = length(σ.times) == length(σ.transitions) == size(σ.values)[1]
    test_obs_var = length_obs_var(σ) == length(σ.m.g) == size(σ.values)[2]
    if !test_sizes error("Issue between sizes of values/times/transitions in this trajectory.") end
    if !test_obs_var error("Issue between sizes of observed variables in model and values") end
    return true
end

# Properties of the trajectory
length_states(σ::AbstractTrajectory) = length(σ.times)
length_obs_var(σ::AbstractTrajectory) = size(σ.values)[2]
get_obs_var(σ::AbstractTrajectory) = (σ.m).g
is_bounded(σ::AbstractTrajectory) = σ.transitions[end] == nothing 

# Access to trajectory values
get_var_values(σ::AbstractTrajectory, var::String) = 
    view(σ.values, :, (σ.m)._map_obs_var_idx[var])
get_state(σ::AbstractTrajectory, idx::Int) = 
    view(σ.values, idx, :)
get_value(σ::AbstractTrajectory, var::String, idx::Int) = 
    σ.values[idx,(σ.m)._map_obs_var_idx[var]] 


δ(σ::AbstractTrajectory,t::Float64) = true

states(σ::AbstractTrajectory) = σ.values
times(σ::AbstractTrajectory) = σ.times
transitions(σ::AbstractTrajectory) = σ.transitions

# Get i-th state [i]
getindex(σ::AbstractTrajectory, idx::Int) = get_state(σ, idx)
# Get i-th value of var ["I", idx]
getindex(σ::AbstractTrajectory, var::String, idx::Int) = get_value(σ, var, idx)
# Get the path of a variable ["I"]
getindex(σ::AbstractTrajectory, var::String) = get_var_values(σ, var)

