
using BenchmarkTools
import BenchmarkTools: mean
using ABCRN
include(get_module_path() * "/src/_tests_simulate.jl")

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

nbr_sim = 10000

@timev simulate(model)
b1 = @benchmark for i = 1:$(nbr_sim) simulate($model) end
@show minimum(b1), mean(b1), maximum(b1)

@timev _simulate_without_view(model)
b1_without_view = @benchmark for i = 1:$(nbr_sim) _simulate_without_view($model) end
@show minimum(b1_without_view), mean(b1_without_view), maximum(b1_without_view)

