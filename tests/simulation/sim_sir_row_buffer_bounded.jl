
using MarkovProcesses 
include(get_module_path() * "/core/_tests_simulate.jl")
using PyPlot

load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = 100.0

σ = _simulate_row_buffer(SIR_row_buffer)
plt.figure()
plt.step(times(σ), _get_values_row(σ,"I"), "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/tests/simulation/res_pics/sim_sir_row_buffer_bounded.png")
plt.close()

return true

