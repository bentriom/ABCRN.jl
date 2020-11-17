
abstract type AbstractTrajectory end
ContinuousObservations = AbstractVector{AbstractTrajectory}

struct Trajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::Matrix{Int}
    times::Vector{Float64}
    transitions::Vector{Union{String,Nothing}}
end

function +(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function -(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function δ(σ1::AbstractTrajectory,t::Float64) end

# Properties of the trajectory
get_states_number(σ::AbstractTrajectory) = length(σ.times)
get_obs_variables(σ::AbstractTrajectory) = (σ.m).g

# Access to trajectory values
get_var_values(σ::AbstractTrajectory, var::String) = 
@view σ.values[:,(σ.m)._map_obs_var_idx[var]] 
get_state(σ::AbstractTrajectory, idx::Int) = @view σ.values[idx,:]
get_value(σ::AbstractTrajectory, var::String, idx::Int) = 
σ.values[idx,(σ.m)._map_obs_var_idx[var]] 

# Get var values ["I"]
function getindex(σ::AbstractTrajectory, var::String)
    if var  == "times"
        return σ.times
    elseif var == "transitions"
        return σ.transitions
    else
        return get_var_values(σ, var)
    end
end
# Get i-th state [i]
getindex(σ::AbstractTrajectory, idx::Int) = get_state(σ, idx)
# Get i-th value of var ["I", idx]
function getindex(σ::AbstractTrajectory, var::String, idx::Int)
    if var  == "times"
        return σ.times[idx]
    elseif var == "transitions"
        return σ.transitions[idx]
    else
        return get_value(σ, var, idx)
    end
end

