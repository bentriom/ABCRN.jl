# From Pudlo: Reliable ABC model choice, 2016, Appendix B

using ARFIMA
using Random
using LinearAlgebra
using BiochemNetABC
using Distributions
using ScikitLearn
@sk_import metrics: (accuracy_score, classification_report)
using StatsBase: autocor

struct MA1 <: Model end
struct MA2 <: Model end
import BiochemNetABC: simulate

global N_tml = 100
global σ = 1.0

struct TriangleDist <: ContinuousMultivariateDistribution end
function Distributions.rand(d::TriangleDist)
    θ1 = rand(Uniform(-2, 2))
    θ2 = (θ1 < 0) ? rand(Uniform(-θ1-1,1)) : rand(Uniform(θ1-1,1))
    return [θ1, θ2]
end
Distributions.rand!(d::TriangleDist, p::AbstractVector) = p[:] = rand(d)
Distributions.length(d::TriangleDist) = 2
Distributions.pdf(d::TriangleDist, p::AbstractVector) = 1/8
function simulate(m::MA1)
    θ1 = rand(Uniform(-1, 1))
    x = zeros(100)
    ϵtm1 = rand(Normal(0,σ^2))
    x[1] = ϵtm1
    for t = 2:100
        ϵt = rand(Normal(0,σ^2))
        x[t] = ϵt - θ1*ϵtm1
        ϵtm1 = ϵt
    end
    return x
end
function simulate(m::MA2)
    θ1, θ2 = rand(TriangleDist())
    x = zeros(100)
    ϵtm1 = rand(Normal(0,σ^2))
    ϵtm2 = rand(Normal(0,σ^2))
    x[1] = ϵtm2
    x[2] = ϵtm1 - θ1*ϵtm2
    for t = 3:100
        ϵt = rand(Normal(0,σ^2))
        x[t] = ϵt - θ1*ϵtm1 - θ2*ϵtm2
        ϵtm2 = ϵtm1
        ϵtm1 = ϵt
    end
    return x
end
#=
function simulate(m::MA1)
    θ1 = rand(Uniform(-1, 1))
    return arma(N_tml, σ, nothing, SVector(θ1))
end
function simulate(m::MA2)
    θ = rand(TriangleDist())
    return arma(N_tml, σ, nothing, SVector(θ[1],θ[2]))
end
=#

m1, m2 = MA1(), MA2()
models = [m1, m2]

ss_func(y) = autocor(y, 1:7)
dist_l2(s_sim,s_obs) = norm(s_sim-s_obs)

observations = simulate(m1)
ss_observations = ss_func(observations)
abc_testset = abc_model_choice_dataset(models, ss_observations, ss_func, dist_l2, 10000, 10000)

grid = Dict(:n_estimators => [300], :min_samples_leaf => [1], :min_samples_split => [2], :n_jobs => [8])
res_rf_abc = rf_abc_model_choice(models, ss_observations, ss_func, 10000; hyperparameters_range = grid)
@show posterior_proba_model(res_rf_abc)
X_testset = transpose(abc_testset.X)
println(classification_report(y_true = abc_testset.y, y_pred = predict(res_rf_abc.clf, X_testset)))
@show accuracy_score(abc_testset.y, predict(res_rf_abc.clf, X_testset))

return true

