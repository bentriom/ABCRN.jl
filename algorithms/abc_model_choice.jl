
struct AbcModelChoiceDataset
    models_indexes::Vector{Int}
    summary_stats_matrix::Matrix
    epsilon::Float64
end

struct RandomForestABC
    reference_table::AbcModelChoiceDataset
    clf
    summary_stats_observations
    estim_model
end

function getproperty(dataset::AbcModelChoiceDataset, sym::Symbol)
    if sym == :X
        return dataset.summary_stats_matrix
    elseif sym == :y
        return dataset.models_indexes
    else
        return getfield(dataset, sym)
    end
end

function abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}},
                                  summary_stats_observations,
                                  summary_stats_func::Function, distance_func::Function,
                                  k::Int, N::Int; dir_results::Union{Nothing,String} = nothing)
    nbr_models = length(models)
    models_prior = Categorical([1/nbr_models for i = 1:nbr_models])
    return abc_model_choice_dataset(models, models_prior, summary_stats_observations, summary_stats_func, distance_func, k, N; dir_results = dir_results)
end

function abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}}, models_prior::DiscreteUnivariateDistribution,
                                  summary_stats_observations,
                                  summary_stats_func::Function, distance_func::Function,
                                  k::Int, N::Int; dir_results::Union{Nothing,String} = nothing)
    @assert length(models) >= 2 "Should contain at least 2 models"
    @assert ncategories(models_prior) == length(models) "Number of categories of models' prior and number of models do not equal"

    models_indexes = zeros(Int, N)
    summary_stats_matrix = zeros(eltype(summary_stats_observations), length(summary_stats_observations), N)
    distances = zeros(N)
    bool_parametric = typeof(models) <: Vector{ParametricModel} 
    for i = 1:N
        current_idx_model = rand(models_prior)
        models_indexes[i] = current_idx_model
        if bool_parametric
            param_model = models[current_idx_model]
            vec_p = rand(param_model.distribution)
            sim = simulate(param_model, vec_p)
        else
            sim = simulate(models[current_idx_model])
        end
        ss_i = summary_stats_func(sim)
        summary_stats_matrix[:,i] = ss_i 
        distances[i] = distance_func(ss_i, summary_stats_observations)
    end
    k_nn = sortperm(distances, alg = QuickSort)[1:k]

    if dir_results != nothing
        dir_results = basename(dir_results) != "" ? dir_results * "/" : dir_results
    end
    return AbcModelChoiceDataset(models_indexes[k_nn], summary_stats_matrix[:,k_nn], distances[k_nn[end]])
end

function rf_abc_model_choice(models::Vector{<:Union{Model,ParametricModel}},
                             summary_stats_observations,
                             summary_stats_func::Function, N_ref::Int;
                             k::Int = N_ref, distance_func::Function = (x,y) -> 1, 
                             hyperparameters_range::Dict = Dict(:n_estimators => [200], :min_samples_leaf => [1],
                                                                :min_samples_split => [2]))
    @assert k <= N_ref
    trainset = abc_model_choice_dataset(models, summary_stats_observations, summary_stats_func, distance_func, k, N_ref)
    gridsearch = GridSearchCV(RandomForestClassifier(oob_score=true), hyperparameters_range)
    fit!(gridsearch, transpose(trainset.X), trainset.y)
    best_rf = gridsearch.best_estimator_
    return RandomForestABC(trainset, best_rf, summary_stats_observations, predict(best_rf, [summary_stats_observations]))
end

# P(m = m^(ss_obs) | ss_obs) estimate
function posterior_proba_model(rf_abc::RandomForestABC)
    oob_votes = rf_abc.clf.oob_decision_function_
    y_pred_oob = argmax.([oob_votes[i,:] for i = 1:size(oob_votes)[1]])
    y_oob_regression = y_pred_oob .!= rf_abc.reference_table.y
    dict_params = Dict()
    for param in ["n_estimators", "min_samples_leaf", "min_samples_split", "oob_score", "n_jobs"]
        dict_params[Symbol(param)] = get_params(rf_abc.clf)[param]
    end
    rf_regressor = RandomForestRegressor(;dict_params...)
    fit!(rf_regressor, transpose(rf_abc.reference_table.X), y_oob_regression)
    return 1 - predict(rf_regressor, [rf_abc.summary_stats_observations])[1]
end

