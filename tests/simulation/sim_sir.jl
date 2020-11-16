
using MarkovProcesses 
using PyPlot

load_model("SIR")

σ = simulate(SIR)
plt.figure()
plt.step(σ["times"], σ["I"], "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/tests/simulation/res_pics/sim_sir.png")
plt.close()

return true

