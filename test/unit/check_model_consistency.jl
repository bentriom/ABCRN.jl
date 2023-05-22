
using ABCRN

load_model("SIR")
load_model("ER")

return check_consistency(ER) && check_consistency(SIR)

