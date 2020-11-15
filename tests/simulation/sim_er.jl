
using MarkovProcesses 
using PyPlot

load_model("ER")

σ = simulate(ER)
plt.figure()
plt.step(σ["times"], σ["P"], "ro--", marker="x", where="post", linewidth=1.0)
plt.savefig("sim_er.png")
plt.close()

