
using ABCRN 
using Plots

load_model("SIR")

σ = simulate(SIR)
plot(times(σ), σ[:I], linetype=:steppost, linewidth=1.0)
savefig(get_module_path() * "/test/simulation/res_pics/sim_sir.svg")

return true

