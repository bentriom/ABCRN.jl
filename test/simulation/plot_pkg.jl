
using ABCRN

load_model("SIR")
σ = simulate(SIR)

load_plots()
plot(σ; plot_transitions = true)

return true

