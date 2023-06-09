
using BiochemNetABC

load_model("SIR")
σ = simulate(SIR)
σ[:I]
σ[:I,2]
σ[3]

return true

