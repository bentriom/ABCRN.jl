
using ABCRN
load_plots()
path_pics = get_module_path() * "/test/simulation/res_pics/"

load_model("doping_3way_oscillator")
observe_all!(doping_3way_oscillator)
load_automaton("period_automaton")
A_per = create_period_automaton(doping_3way_oscillator, 300.0, 360.0, 5, :A)  

sync_doping = doping_3way_oscillator * A_per
set_time_bound!(sync_doping, 0.1)
set_x0!(doping_3way_oscillator, [:A, :B, :C, :DA, :DB, :DC], [333, 333, 333, 10, 10, 10])
σ = simulate(sync_doping)
plot(σ; A = A_per, filename = path_pics * "traj_full_doping.png")
plot(σ, :A; A = A_per, filename = path_pics * "traj_doping_A.png")
plot_periodic_trajectory(A_per, σ, :A, filename = path_pics * "traj_doping_period_lha.png")

return true

