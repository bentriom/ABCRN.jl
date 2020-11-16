
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
function get_obs_variables(σ::AbstractTrajectory) end

get_values(σ::AbstractTrajectory, var::String) = 
σ.values[(σ.m)._map_obs_var_idx[var],:] 

get_states_number(σ::AbstractTrajectory) =
length(σ.times)

function getindex(σ::AbstractTrajectory, idx::String)
    if idx  == "times"
        return σ.times
    elseif idx == "transitions"
        return σ.transitions
    else
        return get_values(σ, idx)
    end
end

