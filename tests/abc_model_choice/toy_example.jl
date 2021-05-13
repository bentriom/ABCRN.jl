
using SpecialFunctions
using LinearAlgebra
using Random
using Distributions
using MarkovProcesses
using ScikitLearn
@sk_import metrics: (accuracy_score, classification_report)

global n = 20

struct Model1 <: Model end
struct Model2 <: Model end
struct Model3 <: Model end
import MarkovProcesses: simulate

function simulate(m::Model1)
    param = rand(Exponential(1))
    return rand(Exponential(param), n)
end
function simulate(m::Model2)
    param = rand(Normal())
    return rand(LogNormal(param,1), n)
end
function simulate(m::Model3)
    param = rand(Exponential(1))
    return rand(Gamma(2,1/param), n)
end

m1, m2, m3 = Model1(), Model2(), Model3()
lh_m1(s) = exp(log(gamma(n+1)) - (n+1)*log(1+s[1]))
lh_m2(s) = exp(-s[2]^2/(2n*(n+1)) - (s[3]^2)/2 + (s[2]^2)/(2n) - s[2]) * (2pi)^(-n/2)*(n+1)^(-1/2)
lh_m3(s) = exp(s[2])*gamma(2n+1)/gamma(2)^n * (1+s[1])^(-2n-1)

ss_func(y) = [sum(y), sum(log.(y)), sum(log.(y).^2)]
dist_l2(s_sim,s_obs) = sqrt(dot(s_sim,s_obs))

observations = simulate(m3)
ss_observations = ss_func(observations)
models = [m1, m2, m3]
abc_testset = abc_model_choice_dataset(models, ss_observations, ss_func, dist_l2, 1000, 1000)

list_lh = [lh_m1, lh_m2, lh_m3]
prob_model(ss, list_lh, idx_model) = list_lh[idx_model](ss) / sum([list_lh[i](ss) for i = eachindex(list_lh)])
prob_model(ss, idx_model) = prob_model(ss, list_lh, idx_model)

#=
using Plots
p = plot(title="Trainset")
colors = ["black", "red", "green"]
begin_idx = 1
for i = 1:3
    models_i = findall(x->x==i, abc_trainset.models_indexes)
    nbr_obs = length(models_i)
    end_idx = begin_idx + nbr_obs - 1
    lh = list_lh[i]
    @show i
    @show nbr_obs
    @show begin_idx:end_idx
    scatter!(p, begin_idx:end_idx, prob_model.(abc_trainset.summary_stats_vector[models_i], 3), 
                color = colors[i], markersize = 3.0, markershape = :cross, label = "Model $i")
    global begin_idx = end_idx + 1
end
savefig("set.svg")
=#

grid = Dict(:n_estimators => [500], :min_samples_leaf => [1], :min_samples_split => [2])
res_rf_abc = rf_abc_model_choice(models, ss_observations, ss_func, 29000; hyperparameters_range = grid)
@show posterior_proba_model(res_rf_abc) 
println(classification_report(y_true = abc_testset.y, y_pred = predict(res_rf_abc.clf, abc_testset.X)))
@show accuracy_score(abc_testset.y, predict(res_rf_abc.clf, abc_testset.X))

