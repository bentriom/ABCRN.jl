
using MarkovProcesses 
include(get_module_path() * "/src/_tests_simulate.jl")
using PyPlot

load_model("_bench_perf_test/ER_row_buffer")
ER_row_buffer.time_bound = 10.0

σ = _simulate_row_buffer(ER_row_buffer)
plt.figure()
plt.step(times(σ), _get_values_row(σ,:P), "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/test/simulation/res_pics/sim_er_row_buffer_bounded.svg")
plt.close()

return true
