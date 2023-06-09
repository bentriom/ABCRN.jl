
using BiochemNetABC

test = true
load_model("SIR")
SIR.time_bound = 120.0
for i = 1:1000
    local σ = simulate(SIR)
    global test = test && isbounded(σ)
end

return test

