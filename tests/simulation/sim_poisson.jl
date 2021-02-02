
using MarkovProcesses
load_plots()

load_model("poisson")
σ = simulate(poisson)
MarkovProcesses.plot(σ; filename = get_module_path() * "/tests/simulation/res_pics/sim_poisson.svg")

return true

