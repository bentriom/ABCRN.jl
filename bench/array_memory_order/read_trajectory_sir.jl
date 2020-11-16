
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/core/_tests_simulate.jl")

bound_time = 200.0
l_var = ["S", "I", "R"]

println("Col buffer:")

MarkovProcesses.load_model("_bench_perf_test/SIR_col_buffer")
SIR_col_buffer.time_bound = bound_time
set_observed_var!(SIR_col_buffer, l_var)
function read_trajectory_col(m::Model)
    res = 0
    σ = _simulate_col_buffer(m)
    n_states = get_states_number(σ)
    n_read = 100000
    for k = 1:n_read
        for i = 1:n_states
            res += _get_value_col(σ, "I", i)
        end
    end
    return res
end
# Bench
@timev read_trajectory_col(SIR_col_buffer) 
b1_col = @benchmark read_trajectory_col($SIR_col_buffer)
@show minimum(b1_col), mean(b1_col), maximum(b1_col)

println("Row buffer:")

MarkovProcesses.load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = bound_time
set_observed_var!(SIR_row_buffer, l_var)
function read_trajectory_row(m::Model)
    res = 0
    σ = _simulate_row_buffer(m)
    n_states = get_states_number(σ)
    n_read = 100000
    for k = 1:n_read
        for i = 1:n_states
            res += _get_value_row(σ, "I", i)
        end
    end
    return res
end
# Bench
@timev read_trajectory_row(SIR_row_buffer) 
b1_row = @benchmark read_trajectory_row($SIR_row_buffer)
@show minimum(b1_row), mean(b1_row), maximum(b1_row)

