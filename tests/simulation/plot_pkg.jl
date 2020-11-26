
using MarkovProcesses

load_model("SIR")
σ = simulate(SIR)

load_plots()
plot(σ; plot_transitions = true)

