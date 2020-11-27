module MarkovProcesses

import Base: +, -, *
import Base: copy, getfield, getindex, getproperty, lastindex
import Base: setindex!, setproperty!, fill!

import StaticArrays: SVector
import Distributions: Distribution, Product, Uniform, Normal

export Distribution, Product, Uniform, Normal

# Common types and constructors
export Observations, AbstractTrajectory, Trajectory, SynchronizedTrajectory
export Model, ContinuousTimeModel, SynchronizedModel, ParametricModel
export LHA, StateLHA, Edge

# Trajectory related methods
export +, -, δ, dist_lp
export get_obs_var, length_states, length_obs_var, get_state_from_time 
export isbounded, times, transitions
export check_consistency, issteadystate, isaccepted

# LHA related methods
export init_state, next_state!, read_trajectory
export get_index, get_value, length_var, isaccepted

# Model related methods
export simulate, set_param!, set_time_bound!, set_observed_var!, observe_all!
export get_param, getproperty, get_proba_model, get_observed_var
export isbounded, isaccepted, check_consistency
export draw_model!, draw!, fill!, prior_pdf!, prior_pdf, insupport

# Utils
export get_module_path, cosmos_get_values
export load_model, load_automaton, load_plots

# Algorithms
export automaton_abc

include("common.jl")
include("trajectory.jl")
include("lha.jl")
include("model.jl")
include("utils.jl")
include("../algorithms/automaton_abc.jl")

end

