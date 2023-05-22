
using ABCRN 
include(get_module_path() * "/src/_tests_simulate.jl")
using Plots

load_model("_bench_perf_test/SIR_col_buffer")
SIR_col_buffer.time_bound = 100.0

σ = _simulate_col_buffer(SIR_col_buffer)
plot(times(σ), _get_values_col(σ,:I), linetype=:steppost, linewidth=1.0)
savefig(get_module_path() * "/test/simulation/res_pics/sim_sir_col_buffer_bounded.svg")

return true

