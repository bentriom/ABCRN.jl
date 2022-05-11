
using MarkovProcesses 
using PyPlot

load_model("ER")

σ = simulate(ER)
plt.figure()
plt.step(times(σ), σ[:P], "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig(get_module_path() * "/test/simulation/res_pics/sim_er.png", dpi=480)
plt.close()

return true

