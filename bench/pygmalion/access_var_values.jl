
using BenchmarkTools
using pygmalion
using MarkovProcesses

println("Pygmalion:")

str_m = "enzymatic_reaction"
str_d = "abc_er"
pygmalion.load_model(str_m)
str_oml = "P,R,time"
ll_om = split(str_oml, ",")
x0 = State(95.0, 95.0, 0.0, 0.0, 0.0, 0.0)
p_true = Parameters(1.0, 1.0, 1.0)
u = Control(10.0)
tml = 1:400
g_all = create_observation_function([ObserverModel(str_oml, tml)]) 
so = pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true)
function read_trajectory_v1(so::SystemObservation)
    vals = to_vec(so, "P")
    n_states = length(vals)
    res = 0.0
    for i = 1:n_states
        res += vals[i]
    end
    return res
end
function read_trajectory_v2(so::SystemObservation)
    idx_P = om_findfirst("P", so.oml)
    n_states = length(so.otll[idx_P]) 
    res = 0.0
    for i = 1:n_states
        res += so.otll[idx_P][i][2][1]
    end
    return res
end
# Bench
@timev read_trajectory_v1(so)
b1_pyg = @benchmark read_trajectory_v1($so) 
@show minimum(b1_pyg), mean(b1_pyg), maximum(b1_pyg)
@timev read_trajectory_v2(so)
b2_pyg = @benchmark read_trajectory_v2($so) 
@show minimum(b2_pyg), mean(b2_pyg), maximum(b2_pyg)

println("MarkovProcesses:")

MarkovProcesses.load_model("ER")
ER.time_bound = 10.0
σ = MarkovProcesses.simulate(ER)
function read_trajectory(σ::AbstractTrajectory)
    n_states = get_states_number(σ)
    res = 0
    for i = 1:n_states
        res += (σ["P"])[i]
    end
    return res
end
# Bench
@timev read_trajectory(σ) 
b1 = @benchmark read_trajectory($σ)
@show minimum(b1), mean(b1), maximum(b1)

