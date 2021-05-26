
using MarkovProcesses

load_model("SIR")
load_automaton("automaton_F")
SIR.time_bound = 120.0
x1, x2, t1, t2 = 0.0, Inf, 100.0, 120.0 

A_F = create_automaton_F(SIR, x1, x2, t1, t2, :I) # <: LHA
sync_SIR = A_F * SIR

function test_last_state(A::LHA, m::SynchronizedModel)
    σ = simulate(m)
    test = (get_state_from_time(σ, (σ.state_lha_end).time)[1] == (σ.state_lha_end)[:n]) && ((σ.state_lha_end)[:d] == 0)
    if !test
        @show σ.state_lha_end
        @show times(σ)[end]
        @show σ[end]
        @show times(σ)[end-1]
        @show σ[end-1]
        @show get_state_from_time(σ, (σ.state_lha_end).time)[1] 
        @show (σ.state_lha_end)[:n], (σ.state_lha_end)[:d]
        error("Ouch")
    end
    return test
end

test_all = true
nbr_sim = 10000
for i = 1:nbr_sim
    local test = test_last_state(A_F, sync_SIR)
    global test_all = test_all && test
end

return test_all

