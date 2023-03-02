module MarkovProcesses

## Imports
import Base: +, -, *
import Base: copy, getfield, getindex, getproperty, lastindex
import Base: setindex!, setproperty!, fill!, copyto!

import Dates
import DelimitedFiles: readdlm, writedlm
import Distributed: @everywhere, @distributed, @sync, @async, nworkers, nprocs, workers
import Distributed: nworkers, nprocs, workers, remotecall_fetch
import DistributedArrays: DArray, dzeros, convert, localpart
import Distributions: Product, Uniform, Normal, MvNormal, Categorical
import Distributions: Distribution, Univariate, Continuous, 
                      UnivariateDistribution, DiscreteUnivariateDistribution,
                      MultivariateDistribution, product_distribution
import Distributions: insupport, isbounded, ncategories, pdf
import FunctionWrappers: FunctionWrapper
import Logging: @info
using LinearAlgebra
using MacroTools
import NearestNeighbors: KDTree, knn 
import Random: rand, rand!
import ScikitLearn
import ScikitLearn: fit!, predict, get_params
import ScikitLearn.GridSearch: GridSearchCV
import SharedArrays: SharedVector, SharedMatrix, sdata
import StaticArrays: SVector, @SVector
import Statistics: quantile
import StatsBase: mean, median, std, cov, ProbabilityWeights
# Python objects import
import PyCall: PyNULL
const RandomForestClassifier = PyNULL()
const RandomForestRegressor = PyNULL()
function __init__()
    (ScikitLearn.Skcore).import_sklearn()
    copy!(RandomForestClassifier, (ScikitLearn.Skcore.pyimport("sklearn.ensemble")).RandomForestClassifier)
    copy!(RandomForestRegressor, (ScikitLearn.Skcore.pyimport("sklearn.ensemble")).RandomForestRegressor)
end

## Exports
export Distribution, Product, Uniform, Normal
export @everywhere

# Common types and constructors
export SVector, @SVector
export Observations, AbstractTrajectory, Trajectory, SynchronizedTrajectory
export Model, ContinuousTimeModel, SynchronizedModel, ParametricModel
export VariableModel, ParameterModel, Transition, TransitionSet
export LHA, StateLHA, Edge, Location, VariableAutomaton
export InvariantPredicateFunction, CheckConstraintsFunction, UpdateStateFunction

# Trajectory related methods
export +, -, δ, dist_lp, euclidean_distance
export get_obs_var, length_states, length_obs_var
export get_state_from_time, get_var_from_time, vectorize, trajectory_from_csv
export isbounded, states, times, transitions
export check_consistency, issteadystate, isaccepted

# LHA related methods
export init_state, next_state!, read_trajectory
export get_index, get_value, length_var, isaccepted

# Model related methods
export simulate, volatile_simulate, change_simulation_stop_criteria
export distribute_mean_value_lha, mean_value_lha, distribute_prob_accept_lha, probability_var_value_lha 
export number_simulations_smc_chernoff, smc_chernoff
export set_param!, set_x0!, set_time_bound!, set_observed_var!, observe_all!
export get_param, get_x0, getproperty, get_proba_model, get_observed_var
export isbounded, isaccepted, check_consistency
export draw_model!, draw!, fill!, prior_pdf!, prior_pdf, insupport

# Utils
export get_module_path, cosmos_get_values
export load_model, load_automaton, load_plots

# Algorithms
export automaton_abc, abc_smc
export abc_model_choice_dataset, rf_abc_model_choice, posterior_proba_model

# About biochemical networks
export @network_model

include("common.jl")
include("trajectory.jl")
include("lha.jl")
include("model.jl")
include("utils.jl")
include("network_model.jl")
include("../algorithms/abc_smc.jl")
include("../algorithms/abc_model_choice.jl")

end

