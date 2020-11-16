
using BenchmarkTools
import BenchmarkTools: mean
using MarkovProcesses
include(get_module_path() * "/core/_tests_simulate.jl")

BenchmarkTools.DEFAULT_PARAMETERS.samples = 20000
l_var = ["S", "I", "R"]
bound_time = 200.0

load_model("_bench_perf_test/SIR_col")
SIR_col.time_bound = bound_time
set_observed_var!(SIR_col, l_var)

println("Col")
b1_col = @benchmark _simulate_col($SIR_col)
@timev _simulate_col(SIR_col)
@show minimum(b1_col), mean(b1_col), maximum(b1_col)

load_model("_bench_perf_test/SIR_col_buffer")
SIR_col_buffer.time_bound = bound_time
set_observed_var!(SIR_col_buffer, l_var)

println("Col + buffer:")
b1_col_buffer = @benchmark _simulate_col_buffer($SIR_col_buffer)
@timev _simulate_col_buffer(SIR_col_buffer)
@show minimum(b1_col_buffer), mean(b1_col_buffer), maximum(b1_col_buffer)

println("Col + buffer_10:")
b1_col_buffer_10 = @benchmark _simulate_col_buffer($SIR_col_buffer; buffer_size = 10)
@timev _simulate_col_buffer(SIR_col_buffer; buffer_size = 10)
@show minimum(b1_col_buffer_10), mean(b1_col_buffer_10), maximum(b1_col_buffer_10)

load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = bound_time
set_observed_var!(SIR_row_buffer, l_var)

println("Row + buffer:")
b1_row_buffer = @benchmark _simulate_row_buffer($SIR_row_buffer)
@timev _simulate_row_buffer(SIR_row_buffer)
@show minimum(b1_row_buffer), mean(b1_row_buffer), maximum(b1_row_buffer)

println("Row + buffer_10:")
b1_row_buffer_10 = @benchmark _simulate_row_buffer($SIR_row_buffer; buffer_size = 10)
@timev _simulate_row_buffer(SIR_row_buffer; buffer_size = 10)
@show minimum(b1_row_buffer_10), mean(b1_row_buffer_10), maximum(b1_row_buffer_10)


