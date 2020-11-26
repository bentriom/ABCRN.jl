
using MarkovProcesses

load_model("ER")
load_model("SIR")

observe_all!(ER)
observe_all!(SIR)

return (ER.g == ["E", "S", "ES", "P"] && SIR.g == ["S", "I", "R"])

