
using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.samples = 20000

nb_sim = 1000

println("Pygmalion:")

using pygmalion
str_m = "enzymatic_reaction"
str_d = "abc_er"
pygmalion.load_model(str_m)
str_oml = "P,R,time"
ll_om = split(str_oml, ",")
x0 = State(95.0, 95.0, 0.0, 0.0, 0.0, 0.0)
p_true = Parameters(0.2, 40.0, 1.0)
u = Control(20.0)
tml = 1:400
g_all = create_observation_function([ObserverModel(str_oml, tml)]) 
@timev pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true)

b1_pyg = @benchmark for i=1:$(nb_sim) pygmalion.simulate($f, $g_all, $x0, $u, $p_true; on = nothing, full_timeline = true) end
#b2_pyg = @benchmark for i=1:$(nb_sim) pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true) end
@show minimum(b1_pyg), mean(b1_pyg), maximum(b1_pyg)
#@show minimum(b2_pyg), mean(b2_pyg), maximum(b2_pyg)

println("MarkovProcesses:")

using MarkovProcesses
MarkovProcesses.load_model("ER")
ER.time_bound = 20.0
set_param!(ER, "k1", 0.2)
set_param!(ER, "k2", 40.0)
@timev MarkovProcesses.simulate(ER)

println("Default buffer size=10")
b1 = @benchmark for i=1:$(nb_sim) MarkovProcesses.simulate($ER) end
#b2 = @benchmark for i=1:$(nb_sim) MarkovProcesses.simulate(ER) end
@show minimum(b1), mean(b1), maximum(b1)
#@show minimum(b2), mean(b2), maximum(b2)
println("Buffer size 100 / min_states 8000")
ER.estim_min_states = 8000 
ER.buffer_size = 100
b1_100 = @benchmark for i=1:$(nb_sim) MarkovProcesses.simulate($ER) end
#b2_100 = @benchmark for i=1:$(nb_sim) MarkovProcesses.simulate(ER) end
@show minimum(b1_100), mean(b1_100), maximum(b1_100)
#@show minimum(b2_100), mean(b2_100), maximum(b2_100)

