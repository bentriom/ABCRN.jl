
using MarkovProcesses
load_model("SIR")
SIR.time_bound = 100.0
load_model("ER")
ER.time_bound = 10.0

test_all = true
nb_sim = 1000

for i = 1:nb_sim
    σ = simulate(SIR)
    σ2 = simulate(ER)
    try
        global test_all = test_all && MarkovProcesses.check_consistency(σ) && MarkovProcesses.check_consistency(σ2)
    catch err
        @show err
        @show length(σ.values), length(σ.times), length(σ.transitions)
        println()
        @show σ.values
        @show σ.times
        @show σ.transitions 
        throw(err)
    end
end

return test_all

