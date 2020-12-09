
using MarkovProcesses

load_model("SIR") # p = 0.0012, 0.05 
load_model("ER") # p = 1; 1; 1 

test_all = true

new_param_ER = [1.0, 2.0, 1.5]
set_param!(ER, new_param_ER)
test_all = test_all && (ER.p == new_param_ER)

k2 = 4.0
set_param!(ER, :k2, k2)
test_all = test_all && (ER.p == [1.0, 4.0, 1.5])

k1 = 0.5
set_param!(ER, :k1, k1)
test_all = test_all && (ER.p == [0.5, 4.0, 1.5])

k3 = 10.0
set_param!(ER, :k3, 10.0)
test_all = test_all && (ER.p == [0.5, 4.0, 10.0])

l_k = [20.0, 10.0]
l_name = [:k3, :k2]
set_param!(ER,  l_name, l_k)
test_all = test_all && (ER.p == [0.5, 10.0, 20.0])

new_param_SIR = [0.0013, 0.08]
set_param!(SIR, new_param_SIR)
test_all = test_all && (SIR.p == new_param_SIR)

kr = 0.06
set_param!(SIR, :kr, kr)
test_all = test_all && (SIR.p == [0.0013, 0.06])

ki = 0.011
set_param!(SIR, :ki, ki)
test_all = test_all && (SIR.p == [0.011, 0.06])

l_name = [:kr, :ki]
l_p = [0.03, 0.015]
set_param!(SIR, l_name, l_p)
test_all = test_all && (SIR.p == [0.015, 0.03])

l_name = [:ki, :kr]
l_p = [0.03, 0.015]
set_param!(SIR, l_name, l_p)
test_all = test_all && (SIR.p == [0.03, 0.015])

l_name = [:kr]
l_p = [0.02]
set_param!(SIR, l_name, l_p)
test_all = test_all && (SIR.p == [0.03, 0.02])

return test_all

