
using MarkovProcesses
load_plots()
path_pics = get_module_path() * "/tests/simulation/res_pics/"

load_model("repressilator")
observe_all!(repressilator)
load_automaton("period_automaton")
A_per = create_period_automaton(repressilator, 100.0, 200.0, 5, :P1)

sync_repressilator = repressilator * A_per
σ = simulate(sync_repressilator)
plot(σ; A = A_per, filename = path_pics * "traj_full_repressilator.png")
plot(σ, :P1; A = A_per, filename = path_pics * "traj_repressilator_A.png")
plot_periodic_trajectory(A_per, σ, :P1, filename = path_pics * "traj_repressilator_period_lha.png")

return true

