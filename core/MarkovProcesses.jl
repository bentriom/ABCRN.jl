module MarkovProcesses

import Base: +, -, *
import Base: copy, getfield, getindex, lastindex, setindex!, getproperty, setproperty!

import StaticArrays: SVector

# Common types and constructors
export Observations, AbstractTrajectory, Trajectory
export LHA, StateLHA, Edge
export Model, ContinuousTimeModel, DiscreteTimeModel

# Trajectory related methods
export +, -, δ, dist_lp
export get_obs_var, length_states, length_obs_var, get_state_from_time 
export is_bounded, times, transitions
export check_consistency, is_steadystate

# LHA related methods
export init_state, next_state!, read_trajectory
export load_automaton, get_index, get_value, length_var

# Model related methods
export simulate, set_param!, get_param, set_observed_var!
export is_bounded, check_consistency
export load_model, get_module_path

# Utils
export get_module_path

include("common.jl")

include("trajectory.jl")
include("lha.jl")
include("model.jl")
include("utils.jl")

end

