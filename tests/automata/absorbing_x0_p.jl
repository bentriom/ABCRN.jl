
using MarkovProcesses

load_model("ER")
observe_all!(ER)
load_automaton("automaton_F")
load_automaton("automaton_G")
load_automaton("automaton_G_and_F")
A_F_R1 = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, "P") 
A_G_R5 = create_automaton_G(ER, 50.0, 100.0, 0.0, 0.8, "E") 
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
x3, x4, t3, t4 = 30.0, 100.0, 0.8, 0.9
A_G_F_R6 = create_automaton_G_and_F(ER, x1, x2, t1, t2, "E",
                                    x3, x4, t3, t4, "P") 

set_param!(ER, [0.0,0.0,0.0])
test_all = true
for A in [A_F_R1, A_G_R5, A_G_F_R6]
    local sync = A*ER
    local σ, σ2 = simulate(sync), simulate(sync)
    global test_all = test_all && isaccepted(σ) && isaccepted(σ2)
end

set_param!(ER, [1.0,1.0,1.0])
set_x0!(ER, [0, 0, 0, 5])
test_all = true
for A in [A_F_R1, A_G_R5, A_G_F_R6]
    local sync = A*ER
    local σ, σ2 = simulate(sync), simulate(sync)
    global test_all = test_all && isaccepted(σ) && isaccepted(σ2)
end

return true

