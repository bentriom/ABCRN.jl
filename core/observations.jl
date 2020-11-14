
abstract type AbstractTrajectory end
ContinuousObservations = AbstractVector{AbstractTrajectory}

struct Trajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::AbstractMatrix{Real}
    times::AbstractMatrix{Real}
    transitions::AbstractVector{Union{String,Missing,Nothing}}
end

struct ObservedTrajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::AbstractMatrix{Real}
    times::AbstractMatrix{Real}
end

function +(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function -(σ1::AbstractTrajectory,σ2::AbstractTrajectory) end
function δ(σ1::AbstractTrajectory,t::Real) end
function get_obs_variables(σ::Trajectory) end
function get_obs_variables(σ::ObservedTrajectory) end
function get_values(σ::AbstractTrajectory, variable::String) end
function get_times(σ::AbstractTrajectory, variable::String) end
function getindex(σ::AbstractTrajectory, idx::String) end

