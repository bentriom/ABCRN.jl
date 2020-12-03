
using MarkovProcesses
# Values x1, x2  t1, t2
str_model = "ER"
load_model(str_model)
model = ER
observe_all!(ER)
ER.buffer_size = 100
load_automaton("automaton_G")
width = 0.5
level = 0.95
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
A_G = create_automaton_G(model, x1, x2, t1, t2, "P")  
set_param!(ER, "k1", 0.2)
set_param!(ER, "k2", 40.0)

test_all = true
sync_ER = ER*A_G
nb_sim = 1000
for i = 1:nb_sim
    local σ = simulate(sync_ER) 
    global test_all = test_all && isaccepted(σ)
    if !isaccepted(σ)
        @show σ
        error("Ouch")
    end
end

return test_all

