module MarkovProcesses

import Base: +, -, getfield, getindex

export Model, ContinuousTimeModel, DiscreteTimeModel
export simulate, set_param!, get_param, set_observed_var!
export is_bounded
export load_model, get_module_path
include("model.jl")

export Observations, AbstractTrajectory, Trajectory
export +, -, δ, dist_lp
export get_obs_var, length_states, length_obs_var, is_bounded
include("trajectory.jl")

end

