
using ABCRN

load_model("SIR")
SIR.time_bound = 100.0
σ1, σ2 = simulate(SIR), simulate(SIR)

dist_lp(σ1, σ2)

return true

