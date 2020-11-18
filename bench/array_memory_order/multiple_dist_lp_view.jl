
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/core/_tests_simulate.jl")

str_model = ARGS[1]
load_model(str_model)
if str_model == "SIR"
    model = SIR
    SIR.time_bound = 150
elseif str_model == "ER"
    model = ER
    ER.time_bound = 10.0
else
    error("Unrecognized model")
end


function compute_dist(m::Model)
    nbr_sim = 10000
    for i = 1:nbr_sim 
        σ1, σ2 = simulate(m), simulate(m)
        dist_lp(σ1, σ2)
    end
end

function compute_dist_without_view(m::Model)
    nbr_sim = 10000
    for i = 1:nbr_sim 
        σ1, σ2 = simulate(m), simulate(m)
        dist_lp(σ1, σ2)
    end
end

@timev simulate(model)
b1 = @benchmark compute_dist($model) 
@show minimum(b1), mean(b1), maximum(b1)

@timev _simulate_without_view(model)
b1_without_view = @benchmark compute_dist_without_view($model) 
@show minimum(b1_without_view), mean(b1_without_view), maximum(b1_without_view)

