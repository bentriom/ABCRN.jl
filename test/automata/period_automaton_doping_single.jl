
using BiochemNetABC
load_plots()

load_model("doping_3way_oscillator")
load_automaton("period_automaton")
set_time_bound!(doping_3way_oscillator, 0.5)
set_x0!(doping_3way_oscillator, [:A, :B, :C], fill(333, 3))
set_x0!(doping_3way_oscillator, [:DA, :DB, :DC], fill(10, 3))

A_per = create_period_automaton(doping_3way_oscillator, 300.0, 360.0, 5, :A; ref_mean_tp = 0.01)
sync_doping = doping_3way_oscillator * A_per

σ = simulate(sync_doping)

test = (σ.state_lha_end[:n] == 5.0 && isapprox(σ.state_lha_end[:mean_tp], 0.01, atol = 0.011))

return test

