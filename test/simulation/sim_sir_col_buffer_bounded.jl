
using MarkovProcesses 
include(get_module_path() * "/src/_tests_simulate.jl")
using PyPlot

load_model("_bench_perf_test/SIR_col_buffer")
SIR_col_buffer.time_bound = 100.0

σ = _simulate_col_buffer(SIR_col_buffer)
plt.figure()
plt.step(times(σ), _get_values_col(σ,:I), "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/test/simulation/res_pics/sim_sir_col_buffer_bounded.svg")
plt.close()

return true

