
using MarkovProcesses 
using PyPlot

load_model("sir")

σ = simulate(SIR)
plt.figure()
plot(σ["S,times"], σ["S,values"])
plt.savefig("sim_sir.png")
plt.close()

