
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/src/_tests_simulate.jl")

BenchmarkTools.DEFAULT_PARAMETERS.samples = 20000
if ARGS[1] == "SIR"
    bound_time = 200.0
    l_var = [:S, :I, :R]

    load_model("_bench_perf_test/SIR_col_buffer")
    SIR_col_buffer.time_bound = bound_time
    set_observed_var!(SIR_col_buffer, l_var)
    load_model("_bench_perf_test/SIR_row_buffer")
    SIR_row_buffer.time_bound = bound_time
    set_observed_var!(SIR_row_buffer, l_var)

    model_col_buffer = SIR_col_buffer
    model_row_buffer = SIR_row_buffer
elseif ARGS[1] == "ER"
    l_var = [:E,:S,:ES,:P]
    bound_time = 20.0

    load_model("_bench_perf_test/ER_col_buffer")
    ER_col_buffer.time_bound = bound_time
    set_observed_var!(ER_col_buffer, l_var)
    load_model("_bench_perf_test/ER_row_buffer")
    ER_row_buffer.time_bound = bound_time
    set_observed_var!(ER_row_buffer, l_var)

    model_col_buffer = ER_col_buffer
    model_row_buffer = ER_row_buffer
else
    error("Unavailable model")
end

println("Col buffer:")

function access_trajectory_col(m::Model)
    res = 0
    n_sim = 100
    for k = 1:n_sim
        σ = _simulate_col_buffer(m)
        t = _get_values_col(σ, :I)
        res += t[end-1]
    end
    return res
end
# Bench
@timev access_trajectory_col(model_col_buffer) 
b1_col = @benchmark access_trajectory_col($model_col_buffer)
@show minimum(b1_col), mean(b1_col), maximum(b1_col)

println("Row buffer:")

function access_trajectory_row(m::Model)
    res = 0
    n_sim = 100
    for k = 1:n_sim
        σ = _simulate_row_buffer(m)
        t = _get_values_row(σ, :I)
        res += t[end-1]
    end
    return res
end
# Bench
@timev access_trajectory_row(model_row_buffer) 
b1_row = @benchmark access_trajectory_row($model_row_buffer)
@show minimum(b1_row), mean(b1_row), maximum(b1_row)

