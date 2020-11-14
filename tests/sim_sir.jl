
using MarkovProcesses 
using PyPlot

include("models/sir.jl")

σ = simulate(SIR)
plt.figure()
plot(σ["S,times"], σ["S,values"])
plt.savefig("sim_sir.png")
plt.close()

