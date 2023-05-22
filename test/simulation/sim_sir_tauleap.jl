
using ABCRN
load_plots()

load_model("SIR_tauleap")
σ = simulate(SIR_tauleap)
plot(σ; filename = get_module_path() * "/test/simulation/res_pics/sim_sir_tauleap.svg")

return true

