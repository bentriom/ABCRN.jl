
using MarkovProcesses

load_model("SIR")
load_model("SIR_tauleap")
load_model("ER")
load_model("poisson")

simulate(SIR)
simulate(ER)
simulate(SIR_tauleap)
simulate(poisson)

return true

