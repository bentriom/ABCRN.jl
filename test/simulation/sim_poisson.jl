
using ABCRN
load_plots()

load_model("poisson")
σ = simulate(poisson)
#plot(σ; filename = get_module_path() * "/test/simulation/res_pics/sim_poisson.svg")
plot(σ; filename = get_module_path() * "/test/simulation/res_pics/sim_poisson.svg")

return true

