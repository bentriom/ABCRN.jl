module MarkovProcesses

import Base: +, -, getfield, getindex

export Model, ContinuousTimeModel, DiscreteTimeModel
export simulate, set_param!, get_param
export load_model
include("model.jl")

export Observations, AbstractTrajectory
export +,-,δ,get_obs_variables,get_states_number
include("observations.jl")

end

