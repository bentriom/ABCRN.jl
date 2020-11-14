
abstract type Model end 
abstract type ContinuousTimeModel <: Model end 
abstract type DiscreteTimeModel <: Model end 

function check_consistency(m::Model) end
function simulate(m::Model, n::Int; bound::Real = Inf)::AbstractObservations end
function set_param!(m::Model, p::AbstractVector{Real})::Nothing end
function get_param(m::Model)::AbstractVector{Real} end

function load_model(name_model::String)
    include(pathof(@__MODULE__) * "/../../models/" * name_model * ".jl")
end

