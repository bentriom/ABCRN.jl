

abstract type Model end 
abstract type ContinuousTimeModel <: Model end
abstract type SynchronizedModel <: Model end
abstract type AbstractTrajectory end

abstract type LHA end
abstract type Edge end

const VariableModel = Symbol
const ParameterModel = Symbol
const Transition = Union{Symbol,Nothing}
const TransitionSet = Union{Vector{Symbol},Nothing}
const Location = Symbol
const VariableAutomaton = Symbol

function generate_code_model_type_def(model_name::Symbol)
    return quote
        mutable struct $(model_name) <: ContinuousTimeModel
            dim_state::Int # state space dim
            dim_params::Int # parameter space dim
            map_var_idx::Dict{VariableModel,Int} # maps variable str to index in the state space
            _map_obs_var_idx::Dict{VariableModel,Int} # maps variable str to index in the observed state space
            map_param_idx::Dict{ParameterModel,Int} # maps parameter str to index in the parameter space
            transitions::Vector{Transition}
            p::Vector{Float64}
            x0::Vector{Int}
            t0::Float64
            f!::Symbol
            g::Vector{VariableModel} # of dimension dim_obs_state
            _g_idx::Vector{Int} # of dimension dim_obs_state
            isabsorbing::Symbol
            time_bound::Float64
            buffer_size::Int
            estim_min_states::Int
        end
    end
end

struct Trajectory <: AbstractTrajectory
    m::ContinuousTimeModel
    values::Vector{Vector{Int}}
    times::Vector{Float64}
    transitions::Vector{Transition}
end

function generate_code_lha_type_def(lha_name::Symbol, edge_type::Symbol)
    return quote
        struct $(lha_name) <: LHA
            transitions::Vector{Transition}
            locations::Vector{Location} 
            Λ::Dict{Location,Function}
            locations_init::Vector{Location}
            locations_final::Vector{Location}
            map_var_automaton_idx::Dict{VariableAutomaton,Int} # nvar keys : str_var => idx in values
            flow::Dict{Location,Vector{Float64}} # output of length nvar
            map_edges::Dict{Location, Dict{Location,Vector{$(edge_type)}}}
            map_edges_transitions::Dict{Location, Dict{Location,Vector{TransitionSet}}}
            map_edges_check_constraints::Dict{Location, Dict{Location,Vector{Function}}}
            map_edges_update_state::Dict{Location, Dict{Location,Vector{Function}}}
            constants::Dict{Symbol,Float64}
            map_var_model_idx::Dict{VariableModel,Int} # of dim d (of a model)
        end
    end
end

mutable struct StateLHA
    A::LHA
    loc::Location
    values::Vector{Float64}
    time::Float64
end

function generate_code_synchronized_model_type_def(model_name::Symbol, lha_name::Symbol)
    synchronized_model_name = Symbol("$(model_name)SynchronizedWith$(lha_name)")
    return quote
        mutable struct $(synchronized_model_name) <: SynchronizedModel
            m::$(model_name)
            automaton::$(lha_name)
        end

        Base.:*(m::$(model_name), A::$(lha_name)) = $(synchronized_model_name)(m, A)
        Base.:*(A::$(lha_name), m::$(model_name)) = $(synchronized_model_name)(m, A)
    end
end

struct SynchronizedTrajectory <: AbstractTrajectory
    state_lha_end::StateLHA
    sm::SynchronizedModel
    values::Vector{Vector{Int}}
    times::Vector{Float64}
    transitions::Vector{Transition}
end

struct ParametricModel
    m::Model
    params::Vector{ParameterModel}
    distribution::Distribution
    _param_idx::Vector{Int}
end

# Constructors
function generate_code_model_type_constructor(model_name::Symbol)
    return quote
        function $(model_name)(dim_state::Int, dim_params::Int, map_var_idx::Dict{VariableModel,Int}, 
                               map_param_idx::Dict{ParameterModel,Int}, transitions::Vector{<:Transition},
                               p::Vector{Float64}, x0::Vector{Int}, t0::Float64,
                               f!::Symbol, isabsorbing::Symbol;
                               g::Vector{VariableModel} = [var for var in keys(map_var_idx)], time_bound::Float64 = Inf, 
                               buffer_size::Int = 10, estim_min_states::Int = 50)
            dim_obs_state = length(g)
            transitions = convert(Vector{Transition}, transitions)
            _map_obs_var_idx = Dict()
            _g_idx = Vector{Int}(undef, dim_obs_state)
            for i = 1:dim_obs_state
                _g_idx[i] = map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::VariableModel => idx in state space )
                _map_obs_var_idx[g[i]] = i
            end
            if length(methods(getfield(Main, f!))) >= 2
                @warn "You have possibly redefined a function Model.f! used in a previously instantiated model."
            end
            if length(methods(getfield(Main, isabsorbing))) >= 2
                @warn "You have possibly redefined a function Model.isabsorbing used in a previously instantiated model."
            end
            new_model = $(model_name)(dim_state, dim_params, map_var_idx, _map_obs_var_idx, map_param_idx, transitions, 
                                      p, x0, t0, f!, g, _g_idx, isabsorbing, time_bound, buffer_size, estim_min_states)
            @assert check_consistency(new_model)
            return new_model
        end
    end
end

LHA(A::LHA, map_var::Dict{VariableModel,Int}) = 
getfield(Main, Symbol(typeof(A)))(A.transitions, A.locations, A.Λ, A.locations_init, A.locations_final, 
                                  A.map_var_automaton_idx, A.flow, A.map_edges, A.constants, map_var)

function ParametricModel(am::Model, priors::Tuple{ParameterModel,UnivariateDistribution}...)
    m = get_proba_model(am)
    params = ParameterModel[]
    distributions = Distribution{Univariate,Continuous}[]
    _param_idx = zeros(Int, 0)
    for prior in priors
        check_vars = true
        str_p = prior[1]
        distribution = prior[2]
        @assert str_p in keys(m.map_param_idx)
        push!(params, str_p)
        push!(distributions, distribution)
        push!(_param_idx, m.map_param_idx[str_p])
    end
    return ParametricModel(am, params, product_distribution(distributions), _param_idx)
end

function ParametricModel(am::Model, params::Vector{ParameterModel}, distribution::MultivariateDistribution)
    m = get_proba_model(am)
    _param_idx = zeros(Int, 0)
    for str_p in params
        push!(_param_idx, m.map_param_idx[str_p])
    end
    return ParametricModel(am, params, distribution, _param_idx)
end

