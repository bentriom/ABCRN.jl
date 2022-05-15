
using MarkovProcesses 
using Plots

load_model("ER")

σ = simulate(ER)
plot(times(σ), σ[:P], linetype=:steppost, linewidth=1.0)
savefig(get_module_path() * "/test/simulation/res_pics/sim_er.png")

return true

