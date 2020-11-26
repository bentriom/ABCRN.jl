module MarkovProcesses

import Base: +, -, *
import Base: copy, getfield, getindex, lastindex, setindex!, getproperty, setproperty!

import StaticArrays: SVector
import Distributions: Distribution, Product, Uniform, Normal

export Distribution, Product, Uniform, Normal

# Common types and constructors
export Observations, AbstractTrajectory, Trajectory, SynchronizedTrajectory
export Model, ContinuousTimeModel, SynchronizedModel, ModelPrior
export LHA, StateLHA, Edge

# Trajectory related methods
export +, -, δ, dist_lp
export get_obs_var, length_states, length_obs_var, get_state_from_time 
export isbounded, times, transitions
export check_consistency, issteadystate, isaccepted

# LHA related methods
export init_state, next_state!, read_trajectory
export load_automaton, get_index, get_value, length_var, isaccepted

# Model related methods
export simulate, set_param!, get_param, set_observed_var!, observe_all!
export set_time_bound!, getproperty, draw!
export isbounded, isaccepted, check_consistency
export load_model, get_module_path

# Utils
export get_module_path, cosmos_get_values, load_plots

include("common.jl")

include("trajectory.jl")
include("lha.jl")
include("model.jl")
include("utils.jl")

end

