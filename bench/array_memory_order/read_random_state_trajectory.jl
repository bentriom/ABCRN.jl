
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

function random_trajectory_value_col(m::ContinuousTimeModel)
    σ = _simulate_col_buffer(m)
    n_states = length_states(σ)
    nb_rand = 1000
    res = 0
    for i = 1:nb_rand
        a = _get_state_col(σ, rand(1:n_states))
        res += a[2]
    end
    return res
end
# Bench
@timev random_trajectory_value_col(model_col_buffer) 
b1_col_buffer = @benchmark random_trajectory_value_col($model_col_buffer)
@show minimum(b1_col_buffer), mean(b1_col_buffer), maximum(b1_col_buffer)

println("Row buffer:")

function random_trajectory_value_row(m::ContinuousTimeModel)
    σ = _simulate_row_buffer(m)
    n_states = length_states(σ)
    nb_rand = 1000
    res = 0
    for i = 1:nb_rand
        a = _get_state_row(σ, rand(1:n_states))
        res += a[2]
    end
    return res
end
# Bench
@timev random_trajectory_value_row(model_row_buffer) 
b1_row_buffer = @benchmark random_trajectory_value_row($model_row_buffer)
@show minimum(b1_row_buffer), mean(b1_row_buffer), maximum(b1_row_buffer)

