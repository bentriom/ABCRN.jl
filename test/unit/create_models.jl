
using MarkovProcesses

load_model("SIR")
load_model("ER")

test_all = true

new_p = [0.01, 0.2]
new_SIR = create_SIR(new_p)
test_all = test_all && new_SIR.p == new_p && new_SIR !== SIR

p_SIR = [0.0012, 0.05]
new_SIR = create_SIR(p_SIR)
test_all = test_all && new_SIR.p == p_SIR && new_SIR !== SIR

new_p = [0.01, 0.2, 20.0]
new_ER = create_ER(new_p)
test_all = test_all && new_ER.p == new_p && new_ER !== ER

p_ER = [1.0, 1.0, 1.0]
new_ER = create_ER(p_ER)
test_all = test_all && new_ER.p == p_ER && new_ER !== ER

return true

