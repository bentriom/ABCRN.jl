
using ABCRN 
include(get_module_path() * "/src/_tests_simulate.jl")
using Plots

load_model("_bench_perf_test/SIR_row_buffer")
SIR_row_buffer.time_bound = 100.0

σ = _simulate_row_buffer(SIR_row_buffer)
plot(times(σ), _get_values_row(σ,:I), linetype=:steppost, linewidth=1.0)
savefig(get_module_path() * "/test/simulation/res_pics/sim_sir_row_buffer_bounded.svg")

return true

