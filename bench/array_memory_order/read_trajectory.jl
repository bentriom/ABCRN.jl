
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/core/_tests_simulate.jl")

BenchmarkTools.DEFAULT_PARAMETERS.samples = 20000
if ARGS[1] == "SIR"
    bound_time = 200.0
    l_var = ["S", "I", "R"]

    load_model("_bench_perf_test/SIR_col_buffer")
    SIR_col_buffer.time_bound = bound_time
    set_observed_var!(SIR_col_buffer, l_var)
    load_model("_bench_perf_test/SIR_row_buffer")
    SIR_row_buffer.time_bound = bound_time
    set_observed_var!(SIR_row_buffer, l_var)

    model_col_buffer = SIR_col_buffer
    model_row_buffer = SIR_row_buffer
elseif ARGS[1] == "ER"
    l_var = ["E","S","ES","P"]
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

function read_trajectory_col(m::Model)
    res = 0
    σ = _simulate_col_buffer(m)
    n_states = length_states(σ)
    n_read = 100000
    for k = 1:n_read
        for i = 1:n_states
            res += _get_value_col(σ, "I", i)
        end
    end
    return res
end
# Bench
@timev read_trajectory_col(model_col_buffer) 
b1_col = @benchmark read_trajectory_col($model_col_buffer)
@show minimum(b1_col), mean(b1_col), maximum(b1_col)

println("Row buffer:")

function read_trajectory_row(m::Model)
    res = 0
    σ = _simulate_row_buffer(m)
    n_states = length_states(σ)
    n_read = 100000
    for k = 1:n_read
        for i = 1:n_states
            res += _get_value_row(σ, "I", i)
        end
    end
    return res
end
# Bench
@timev read_trajectory_row(model_row_buffer) 
b1_row = @benchmark read_trajectory_row($model_row_buffer)
@show minimum(b1_row), mean(b1_row), maximum(b1_row)

