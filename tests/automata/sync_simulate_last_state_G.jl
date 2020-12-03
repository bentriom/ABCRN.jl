
using MarkovProcesses

load_model("ER")
load_automaton("automaton_G")
ER.time_bound = 2.0
x1, x2, t1, t2 = 0.0, Inf, 0.0, 2.0 

A_G = create_automaton_G(ER, x1, x2, t1, t2, "P") # <: LHA
sync_ER = A_G * ER

function test_last_state(A::LHA, m::SynchronizedModel)
    σ = simulate(m)
    test = (get_state_from_time(σ, (σ.state_lha_end).time)[1] == (σ.state_lha_end)["n"]) && ((σ.state_lha_end)["d"] == 0)
    return test
end

test_all = true
nbr_sim = 10000
for i = 1:nbr_sim
    local test = test_last_state(A_G, sync_ER)
    global test_all = test_all && test
end

return test_all

