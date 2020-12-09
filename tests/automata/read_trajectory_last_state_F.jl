
using MarkovProcesses

load_model("SIR")
load_automaton("automaton_F")
SIR.time_bound = 120.0
observe_all!(SIR)
x1, x2, t1, t2 = 0.0, Inf, 100.0, 120.0 

A_F = create_automaton_F(SIR, x1, x2, t1, t2, "I") # <: LHA

function test_last_state(A::LHA, m::ContinuousTimeModel)
    σ = simulate(m)
    Send = read_trajectory(A, σ)
    test = (get_state_from_time(σ, Send.time)[2] == Send[:n]) && (Send[:d] == 0)
    if !test
        @show Send
        @show get_state_from_time(σ, Send.time)
        error("tkt")
    end
    return test
end

test_all = true
nbr_sim = 10000
for i = 1:nbr_sim
    local test = test_last_state(A_F, SIR)
    global test_all = test_all && test
end

return test_all

