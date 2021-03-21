
## About distances

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
- `var::Symbol` is an observed variable. Have to be contained in `get_obs_var(σ1)` and `get_obs_var(σ2)`.
- `verbose::Bool` If `true`, launch a verbose execution of the computation. 
...
"""
function dist_lp(σ1::AbstractTrajectory, σ2::AbstractTrajectory, var::VariableModel; 
                 verbose::Bool = false, p::Int = 1)
    if !isbounded(σ1) || !isbounded(σ2)
        @warn "Lp distance computed on unbounded trajectories. Result should be wrong"
    end
    return dist_lp(σ1[var], times(σ1), σ2[var], times(σ2); verbose = false, p = p)
end

# Distance function. Vectorized version
function dist_lp(x_obs::Vector{Int}, t_x::Vector{Float64}, y_obs::Vector{Int}, t_y::Vector{Float64}; 
                 verbose::Bool = false, p::Int = 1)
    current_y_obs = y_obs[1]
    current_t_y = t_y[2]
    nbr_y_obs = length(y_obs)
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
        while current_t_y < t_x[i+1] && idx <= nbr_y_obs
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

function vectorize(σ::AbstractTrajectory, sym_var::Symbol, 
                   timeline::AbstractVector{Float64})
    times_var = times(σ)
    time_end = isbounded(σ) ? times_var[end] : Inf
    @assert timeline[end] <= time_end "Trajectory is bounded and timeline is out of bounds."
    states_var = σ[sym_var]
    nbr_states = length(states_var)
    nbr_points = length(timeline)
    trajectory_points = zeros(nbr_points)
    index_timeline = 1
    for i = eachindex(states_var)
        next_transition_time = (i < nbr_states) ? times_var[i+1] : time_end*1.01
        while index_timeline <= nbr_points && timeline[index_timeline] < next_transition_time
            trajectory_points[index_timeline] = states_var[i]
            index_timeline += 1
        end
        if index_timeline > nbr_points
            break
        end
    end
    return trajectory_points
end

function euclidean_distance(σ::AbstractTrajectory, sym_obs::Symbol,
                            tml::AbstractVector{Float64}, observations::AbstractVector{Float64})
    traj_obs = vectorize(σ, sym_obs, tml)
    diff_obs = observations - traj_obs
    return sqrt(dot(diff_obs, diff_obs))
end

function check_consistency(σ::AbstractTrajectory)
    test_length_var = true
    for i = 1:σ.m.dim_obs_state 
        test_length_i = (length(σ.values[1]) == length(σ.values[i]))
        test_length_var = test_length_var && test_length_i 
    end
    @assert begin (length(σ.times) == length(σ.transitions)) && 
                  (length(σ.times) == length(σ.values[1])) &&
                  test_length_var
            end
    @assert length_obs_var(σ) == σ.m.dim_obs_state
    return true
end

function Base.show(io::IO, σ::Trajectory)
    print(io, "Trajectory of $(σ.m)\n")
    print(io, "- Variable trajectories:\n")
    for obs_var in σ.m.g
        print(io, "* $obs_var: $(σ[obs_var])\n")
    end
    print(io, "- times  = $(times(σ))\n")
    print(io, "- transitions  = $(transitions(σ))")
end
function Base.show(io::IO, σ::SynchronizedTrajectory)
    print(io, "SynchronizedTrajectory\n")
    print(io, "End LHA state:\n")
    print(io, σ.state_lha_end)
    print(io, "\n")
    print(io, "- Model: $(typeof(σ.m)) \n")
    print(io, "- Variable trajectories:\n")
    for obs_var in σ.m.g
        print(io, "* $obs_var: $(σ[obs_var])\n")
    end
    print(io, "- times  = $(times(σ))\n")
    print(io, "- transitions  = $(transitions(σ))")
end

# Properties of the trajectory
length_states(σ::AbstractTrajectory) = length(σ.times)
length_obs_var(σ::AbstractTrajectory) = length(σ.values)
get_obs_var(σ::AbstractTrajectory) = σ.m.g
isbounded(σ::AbstractTrajectory) = σ.transitions[end] == nothing && length_states(σ) >= 2
isaccepted(σ::SynchronizedTrajectory) = isaccepted(σ.state_lha_end)
issteadystate(σ::AbstractTrajectory) = @warn "Unimplemented"

# Access to trajectory values
get_var_values(σ::AbstractTrajectory, var::VariableModel) = σ.values[(σ.m)._map_obs_var_idx[var]]
get_state(σ::AbstractTrajectory, idx::Int) = [σ.values[i][idx] for i = 1:length(σ.values)] # /!\ Creates an array
get_value(σ::AbstractTrajectory, var::VariableModel, idx::Int) = get_var_values(σ, var)[idx]
# Operation σ@t
function get_state_from_time(σ::AbstractTrajectory, t::Float64)
    @assert t >= 0.0
    l_t = times(σ)
    if t == l_t[end] return σ[end] end
    if t > l_t[end]
        if !isbounded(σ)
            return σ[end]
        else 
            error("This trajectory is bounded and you're accessing out of the bounds")
        end
    end
    for i in eachindex(l_t)
        if l_t[i] <= t < l_t[i+1]
            return σ[i]
        end
    end
    error("Unexpected behavior")
end
# Operation σ@t
function get_var_from_time(σ::AbstractTrajectory, sym_var::Symbol, t::Float64)
    @assert t >= 0.0
    l_t = times(σ)
    if t == l_t[end] return σ[sym_var][end] end
    if t > l_t[end]
        if !isbounded(σ)
            return σ[sym_var][end]
        else 
            error("This trajectory is bounded and you're accessing out of the bounds")
        end
    end
    for i in eachindex(l_t)
        if l_t[i] <= t < l_t[i+1]
            return σ[sym_var][i]
        end
    end
    error("Unexpected behavior")
end

function getproperty(σ::SynchronizedTrajectory, sym::Symbol)
    if sym == :sm
        return getfield(σ, :sm)
    elseif sym == :m
        return getfield(σ.sm, :m)
    elseif sym in keys((σ.sm).m.map_var_idx)
        return get_var_values(σ, sym)
    else
        return getfield(σ, sym)
    end
end
function getproperty(σ::Trajectory, sym::Symbol)
    if sym == :m
        return getfield(σ, :m)
    elseif sym in keys((σ.m).map_var_idx)
        return get_var_values(σ, sym)
    else
        return getfield(σ, sym)
    end
end

states(σ::AbstractTrajectory) = σ.values
times(σ::AbstractTrajectory) = σ.times
transitions(σ::AbstractTrajectory) = σ.transitions

# Get i-th state [i]
getindex(σ::AbstractTrajectory, idx::Int) = get_state(σ, idx)
# Get i-th value of var [:I, idx]
getindex(σ::AbstractTrajectory, var::VariableModel, idx::Int) = get_value(σ, var, idx)
# Get the path of a variable [:I]
getindex(σ::AbstractTrajectory, var::VariableModel) = get_var_values(σ, var)
lastindex(σ::AbstractTrajectory) = length_states(σ)

# Operations
function +(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function -(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
δ(σ::AbstractTrajectory,idx::Int) = times(σ)[i+1] - times(σ)[i]

function trajectory_from_csv(csv_file, model::ContinuousTimeModel)
    csv_mat_values, header = readdlm(csv_file, ',', header = true)
    nbr_states = size(csv_mat_values, 1)
    times = zeros(nbr_states)
    values = Vector{Vector{Int}}(undef, length(model.g))
    transitions = fill(nothing, nbr_states)
    for i = eachindex(header)
        model_var = header[i]
        if model_var == "time"
            times = csv_mat_values[:,i]
        elseif model_var == "transitions"
            transitions = csv_mat_values[:,i]
        else
            @assert Symbol(model_var) in model.g "Variable is not observed in the model"
            values[model._map_obs_var_idx[Symbol(model_var)]] = csv_mat_values[:,i]
        end
    end
    return Trajectory(model, values, times, transitions)
end

