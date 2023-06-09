
using BiochemNetABC

load_model("SIR")
load_model("SIR_tauleap")
load_model("ER")
load_model("poisson")
load_model("intracellular_viral_infection")
load_model("repressilator")
load_model("doping_3way_oscillator")
load_model("square_wave_oscillator")

simulate(SIR)
simulate(ER)
simulate(SIR_tauleap)
simulate(poisson)
simulate(intracellular_viral_infection)
simulate(repressilator)
simulate(doping_3way_oscillator)
simulate(square_wave_oscillator)

return true

