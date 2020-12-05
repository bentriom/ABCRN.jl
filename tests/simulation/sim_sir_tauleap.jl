
using MarkovProcesses
load_plots()

load_model("SIR_tauleap")
σ = simulate(SIR_tauleap)
MarkovProcesses.plot(σ; filename = get_module_path() * "/tests/simulation/res_pics/sim_sir_tauleap.png")

return true

