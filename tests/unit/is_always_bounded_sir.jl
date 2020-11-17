
using MarkovProcesses

test = true
load_model("SIR")
SIR.time_bound = 120.0
for i = 1:1000
    σ = simulate(SIR)
    global test = test && is_bounded(σ)
end

return test

