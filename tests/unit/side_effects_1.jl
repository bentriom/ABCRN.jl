
using MarkovProcesses

load_model("SIR")
load_model("ER")

test = SIR.map_var_idx !== ER.map_var_idx &&
       SIR.map_param_idx !== ER.map_var_idx &&
       SIR.l_transitions !== ER.l_transitions &&
       SIR.g !== ER.g &&
       SIR._g_idx !== ER._g_idx &&
       SIR.f! != ER.f! &&
       SIR.isabsorbing != ER.isabsorbing

return test

