
using BenchmarkTools
using pygmalion
using MarkovProcesses

println("Pygmalion:")

str_m = "enzymatic_reaction"
str_d = "abc_er"
pygmalion.load_model(str_m)
l_var = ["E", "S", "ES", "P"]
str_oml = "E,S,ES,P,R,time"
ll_om = split(str_oml, ",")
x0 = State(95.0, 95.0, 0.0, 0.0, 0.0, 0.0)
p_true = Parameters(1.0, 1.0, 1.0)
u = Control(10.0)
tml = 1:400
g_all = create_observation_function([ObserverModel(str_oml, tml)]) 
so = pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true)
function random_trajectory_value_pyg(so::SystemObservation)
    n_states = get_number_of_observations(so, "P")
    return to_vec(so, "P", rand(1:n_states))
end
function random_trajectory_state_pyg(so::SystemObservation)
    n_states = get_number_of_observations(so, "P")
    return to_vec(so, so.oml, rand(1:n_states))
end
# Bench
@timev random_trajectory_value_pyg(so)
b1_pyg = @benchmark random_trajectory_value_pyg($so)
@show minimum(b1_pyg), mean(b1_pyg), maximum(b1_pyg)

println("MarkovProcesses:")

MarkovProcesses.load_model("ER")
ER.time_bound = 10.0
σ = MarkovProcesses.simulate(ER)
function random_trajectory_value(σ::AbstractTrajectory)
    n_states = get_states_number(σ)
    return σ["P"][rand(1:n_states)]
end
# Bench
@timev random_trajectory_value(σ) 
b1 = @benchmark random_trajectory_value($σ)
@show minimum(b1), mean(b1), maximum(b1)

