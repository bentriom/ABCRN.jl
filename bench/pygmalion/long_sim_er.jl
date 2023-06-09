
using BenchmarkTools

println("Pygmalion:")

using pygmalion
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
b1_pyg = @benchmark pygmalion.simulate($f, $g_all, $x0, $u, $p_true; on = nothing, full_timeline = true)
b2_pyg = @benchmark pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true)

@timev pygmalion.simulate(f, g_all, x0, u, p_true; on = nothing, full_timeline = true)
@show minimum(b1_pyg), mean(b1_pyg), maximum(b1_pyg)
@show minimum(b2_pyg), mean(b2_pyg), maximum(b2_pyg)

println("BiochemNetABC:")

using BiochemNetABC
BiochemNetABC.load_model("ER")
ER.time_bound = 10.0
b1 = @benchmark BiochemNetABC.simulate($ER)
b2 = @benchmark BiochemNetABC.simulate(ER)

@timev BiochemNetABC.simulate(ER)
@show minimum(b1), mean(b1), maximum(b1)
@show minimum(b2), mean(b2), maximum(b2)

