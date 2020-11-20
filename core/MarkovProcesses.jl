module MarkovProcesses

import Base: +, -, *
import Base: copy, getfield, getindex, lastindex, setindex!, getproperty, setproperty!

export Model, ContinuousTimeModel, DiscreteTimeModel
export simulate, set_param!, get_param, set_observed_var!
export is_bounded
export load_model, get_module_path
include("model.jl")

export Observations, AbstractTrajectory, Trajectory
export +, -, δ, dist_lp
export get_obs_var, length_states, length_obs_var, get_state_from_time 
export is_bounded, times, transitions
export check_consistency, is_steadystate
include("trajectory.jl")

export LHA, StateLHA, Edge
export init_state, next_state!, read_trajectory
export load_automaton, get_index, get_value, length_var
include("lha.jl")

end

