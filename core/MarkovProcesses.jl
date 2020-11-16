module MarkovProcesses

import Base: +, -, getfield, getindex

export Model, ContinuousTimeModel, DiscreteTimeModel
export simulate, set_param!, get_param, set_observed_var!
export is_bounded
export load_model, get_module_path
include("model.jl")

export Observations, AbstractTrajectory, Trajectory
export +,-,δ,get_obs_variables,get_states_number
include("observations.jl")

end

