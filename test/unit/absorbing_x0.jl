
using BiochemNetABC
load_model("SIR")
new_x0 = [95, 0, 0]
SIR.x0 = new_x0

σ = simulate(SIR)
check_consistency(σ)

return length_states(σ) == 1 && length(times(σ)) == 1 && length(transitions(σ)) == 1

