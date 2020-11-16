
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/core/_tests_simulate.jl")

l_var = ["S", "I", "R"]
bound_time = 200.0

println("Col buffer:")

load_model("_bench_perf_test/SIR_col_buffer")
SIR_col_buffer.time_bound = bound_time
set_observed_var!(SIR_col_buffer, l_var)
function random_trajectory_value_col(m::ContinuousTimeModel)
    σ = _simulate_col_buffer(m)
    n_states = get_states_number(σ)
    nb_rand = 1000
    res = 0
    for i = 1:nb_rand
        a = _get_state_col(σ, rand(1:n_states))
        res += a[2]
    end
    return res
end
# Bench
@timev random_trajectory_value_col(SIR_col_buffer) 
b1_col_buffer = @benchmark random_trajectory_value_col($SIR_col_buffer)
@show minimum(b1_col_buffer), mean(b1_col_buffer), maximum(b1_col_buffer)

println("Row buffer:")

load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = bound_time
set_observed_var!(SIR_row_buffer, l_var)
function random_trajectory_value_row(m::ContinuousTimeModel)
    σ = _simulate_row_buffer(m)
    n_states = get_states_number(σ)
    nb_rand = 1000
    res = 0
    for i = 1:nb_rand
        a = _get_state_row(σ, rand(1:n_states))
        res += a[2]
    end
    return res
end
# Bench
@timev random_trajectory_value_row(SIR_row_buffer) 
b1_row_buffer = @benchmark random_trajectory_value_row($SIR_row_buffer)
@show minimum(b1_row_buffer), mean(b1_row_buffer), maximum(b1_row_buffer)

