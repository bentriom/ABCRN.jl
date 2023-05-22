
using ABCRN 
include(get_module_path() * "/src/_tests_simulate.jl")
using Plots

load_model("_bench_perf_test/ER_row_buffer")
ER_row_buffer.time_bound = 10.0

σ = _simulate_row_buffer(ER_row_buffer)
plot(times(σ), _get_values_row(σ,:P), linetype=:steppost, linewidth=1.0)
savefig(get_module_path() * "/test/simulation/res_pics/sim_er_row_buffer_bounded.svg")

return true

