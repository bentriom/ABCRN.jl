
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

"""
    abc_model_choice_dataset(models, models_prior,
                             summary_stats_observations,
                             summary_stats_func::Function, distance_func::Function,
                             k::Int, N_ref::Int; dir_results::Union{Nothing,String} = nothing)

Creates a reference table for ABC model choice.

The mandatory arguments are:
* `models` is a list of objects inherited from `Model` or `ParametricModel`,
* `models_prior`: the prior over the models (by default: discrete uniform distribution)
* `summary_stats_observations` are the summary statitics of the observations,
* `summary_stats_func::Function`: the function that computes the summary statistics over a model simulation,
* `distance_func`: the distance function over the summary statistics space,
* `N_ref`: the number of samples in the reference table,
* `k`: the k nearest samples from the observations to keep in the reference table (k < N_ref).

The result is a `AbcModelChoiceDataset` with fields:
* `summary_stats_matrix`: the (N_stats, N_ref) features matrix. Accessible via `.X`.
* `models_indexes`: the labels vector. Accessible via `.y`.

If specified, `dir_results` is the directory where the summary statistics matrix and associated models are stored (CSV).
"""
function abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}}, models_prior::DiscreteUnivariateDistribution,
                                  summary_stats_observations,
                                  summary_stats_func::Function, distance_func::Function,
                                  k::Int, N_ref::Int; dir_results::Union{Nothing,String} = nothing)
    if nprocs() == 1
        return _abc_model_choice_dataset(models, models_prior, summary_stats_observations, summary_stats_func, distance_func, k, N_ref; dir_results = dir_results)
    end
    return _distributed_abc_model_choice_dataset(models, models_prior, summary_stats_observations, summary_stats_func, distance_func, k, N_ref; dir_results = dir_results)
end

"""
    abc_model_choice_dataset(models,
                             summary_stats_observations,
                             summary_stats_func::Function, distance_func::Function,
                             k::Int, N_ref::Int; dir_results::Union{Nothing,String} = nothing)

Creates a reference table for ABC model choice with discrete uniform prior distribution over the models.
"""
function abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}},
                                  summary_stats_observations,
                                  summary_stats_func::Function, distance_func::Function,
                                  k::Int, N_ref::Int; dir_results::Union{Nothing,String} = nothing)
    nbr_models = length(models)
    models_prior = Categorical([1/nbr_models for i = 1:nbr_models])
    return abc_model_choice_dataset(models, models_prior, summary_stats_observations, summary_stats_func, distance_func, k, N_ref; dir_results = dir_results)
end

function _abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}}, models_prior::DiscreteUnivariateDistribution,
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
    knn_models_indexes = models_indexes[k_nn]
    knn_summary_stats_matrix = summary_stats_matrix[:,k_nn]
    if dir_results != nothing
        dir_results = basename(dir_results) != "" ? dir_results * "/" : dir_results
        if !isdir(dir_results) mkdir(dir_results) end
        writedlm(dir_results * "X_abc_dataset.csv", knn_summary_stats_matrix, ',')
        writedlm(dir_results * "y_abc_dataset.csv", knn_models_indexes, ',')
        file_cfg = open(dir_results * "config_abc_dataset.out", "w")
        write(file_cfg, "Models: $(models) \n")
        write(file_cfg, "N: $(N) \n")
        write(file_cfg, "k: $(k) \n")
        write(file_cfg, "tolerance epsilon: $(distances[k_nn[end]]) \n")
        close(file_cfg)
    end

    return AbcModelChoiceDataset(knn_models_indexes, knn_summary_stats_matrix, distances[k_nn[end]])
end

function _distributed_abc_model_choice_dataset(models::Vector{<:Union{Model,ParametricModel}}, models_prior::DiscreteUnivariateDistribution,
                                              summary_stats_observations,
                                              summary_stats_func::Function, distance_func::Function,
                                              k::Int, N::Int; dir_results::Union{Nothing,String} = nothing)
    @assert length(models) >= 2 "Should contain at least 2 models"
    @assert ncategories(models_prior) == length(models) "Number of categories of models' prior and number of models do not equal"

    models_indexes = SharedVector{Int}(N)
    summary_stats_matrix = SharedMatrix{eltype(summary_stats_observations)}(length(summary_stats_observations), N)
    distances = SharedVector{Float64}(N)
    bool_parametric = typeof(models) <: Vector{ParametricModel} 
    @sync @distributed for i = 1:N
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
    k_nn = sortperm(sdata(distances), alg = QuickSort)[1:k]
    knn_models_indexes = sdata(models_indexes)[k_nn]
    knn_summary_stats_matrix = sdata(summary_stats_matrix)[:,k_nn]
    if dir_results != nothing
        dir_results = basename(dir_results) != "" ? dir_results * "/" : dir_results
        if !isdir(dir_results) mkdir(dir_results) end
        writedlm(dir_results * "X_abc_dataset.csv", knn_summary_stats_matrix, ',')
        writedlm(dir_results * "y_abc_dataset.csv", knn_models_indexes, ',')
        file_cfg = open(dir_results * "config_abc_dataset.out", "w")
        write(file_cfg, "Models: $(models) \n")
        write(file_cfg, "N: $(N) \n")
        write(file_cfg, "k: $(k) \n")
        write(file_cfg, "tolerance epsilon: $(distances[k_nn[end]]) \n")
        close(file_cfg)
    end

    return AbcModelChoiceDataset(knn_models_indexes, knn_summary_stats_matrix, distances[k_nn[end]])
end

"""
    rf_abc_model_choice(models, summary_stats_observations,
                        summary_stats_func::Function, N_ref::Int;
                        k::Int = N_ref, distance_func::Function = (x,y) -> 1, 
                        hyperparameters_range::Dict)

Run the Random Forest Approximate Bayesian Computation model choice method.

The mandatory arguments are:
* `models` is a list of objects inherited from `Model` or `ParametricModel`,
* `summary_stats_observations` are the summary statitics of the observations
* `N_ref`: the number of samples in the reference table.
* `summary_stats_func::Function`: the function that computes the summary statistics over a model simulation.

The optional arguments are:
* `models_prior`: the prior over the models (by default: discrete uniform distribution)
* `k`: the k nearest samples from the observations to keep in the reference table (by default: k = N_ref)
* `distance_func`: the distance function, has to be defined if k < N_ref
* `hyperparameters_range`: a dict with the hyperparameters range values for the cross validation
    fit of the Random Forest (by default: `Dict(:n_estimators => [200], :min_samples_leaf => [1], :min_samples_split => [2])`).
    See scikit-learn documentation of RandomForestClassifier for the hyperparameters name.

The result is a `RandomForestABC` object with fields:
* `reference_table` an AbcModelChoiceDataset that corresponds to the reference table of the algorithm, 
* `clf` a random forest classifier (PyObject from scikit-learn),
* `summary_stats_observations` are the summary statitics of the observations
* `estim_model` is the underlying model of the observations inferred with the RF-ABC method.

"""
function rf_abc_model_choice(models::Vector{<:Union{Model,ParametricModel}},
                             summary_stats_observations,
                             summary_stats_func::Function, N_ref::Int;
                             models_prior::DiscreteUnivariateDistribution = Categorical([1/length(models) for i = 1:length(models)]),
                             k::Int = N_ref, distance_func::Function = (x,y) -> 1, 
                             hyperparameters_range::Dict = Dict(:n_estimators => [200], :min_samples_leaf => [1],
                                                                :min_samples_split => [2]))
    @assert k <= N_ref
    trainset = abc_model_choice_dataset(models, models_prior, summary_stats_observations, summary_stats_func, distance_func, k, N_ref)
    gridsearch = GridSearchCV(RandomForestClassifier(oob_score=true), hyperparameters_range)
    fit!(gridsearch, transpose(trainset.X), trainset.y)
    best_rf = gridsearch.best_estimator_
    return RandomForestABC(trainset, best_rf, summary_stats_observations, predict(best_rf, [summary_stats_observations]))
end

"""
    posterior_proba_model(rf_abc::RandomForestABC)

Estimates the posterior probability of the model ``P(M = \\widehat{M}(s_{obs}) | s_{obs})`` with the Random Forest ABC method.
"""
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

