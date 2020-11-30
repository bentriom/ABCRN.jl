
using MarkovProcesses

load_model("ER")
load_automaton("automaton_G")
observe_all!(ER)
ER.time_bound = 2.0
x1, x2, t1, t2 = 0.0, Inf, 0.0, 2.0

A_G = create_automaton_G(ER, x1, x2, t1, t2, "P") # <: LHA

function test_last_state(A::LHA, m::ContinuousTimeModel)
    σ = simulate(m)
    Send = read_trajectory(A, σ)
    test = (get_state_from_time(σ, Send.time)[4] == Send["n"]) && (Send["d"] == 0)
    if !test
        @show Send
        @show get_state_from_time(σ, Send.time)
        error("bad")
    end
    return test
end

test_all = true
nbr_sim = 10000
for i = 1:nbr_sim
    test = test_last_state(A_G, ER)
    global test_all = test_all && test
end

return test_all

