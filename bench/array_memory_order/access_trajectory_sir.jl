
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
function access_trajectory_col(m::Model)
    res = 0
    n_sim = 100
    for k = 1:n_sim
        σ = _simulate_col_buffer(m)
        t = _get_values_col(σ, "I")
        res += t[end-1]
    end
    return res
end
# Bench
@timev access_trajectory_col(SIR_col_buffer) 
b1_col = @benchmark access_trajectory_col($SIR_col_buffer)
@show minimum(b1_col), mean(b1_col), maximum(b1_col)

println("Row buffer:")

MarkovProcesses.load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = bound_time
set_observed_var!(SIR_row_buffer, l_var)
function access_trajectory_row(m::Model)
    res = 0
    n_sim = 100
    for k = 1:n_sim
        σ = _simulate_row_buffer(m)
        t = _get_values_row(σ, "I")
        res += t[end-1]
    end
    return res
end
# Bench
@timev access_trajectory_row(SIR_row_buffer) 
b1_row = @benchmark access_trajectory_row($SIR_row_buffer)
@show minimum(b1_row), mean(b1_row), maximum(b1_row)

