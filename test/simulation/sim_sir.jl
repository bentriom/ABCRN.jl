
using MarkovProcesses 
using PyPlot

load_model("SIR")

σ = simulate(SIR)
plt.figure()
plt.step(times(σ), σ[:I], "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/test/simulation/res_pics/sim_sir.svg")
plt.close()

return true
