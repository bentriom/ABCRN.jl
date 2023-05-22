
using BenchmarkTools
@everywhere using ABCRN
using Profile

load_model("ER")
observe_all!(ER)
set_param!(ER, [:k1, :k2], [0.2, 40.0])
set_time_bound!(ER, 0.9)

load_automaton("automaton_G_and_F")
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
x3, x4, t3, t4 = 30.0, 100.0, 0.8, 0.9
A_G_F = create_automaton_G_and_F(ER, x1, x2, t1, t2, :E,
                                 x3, x4, t3, t4, :P)
sync_ER = ER * A_G_F

println("Sim ER")
b_sim = @benchmark simulate(ER)
@btime simulate(ER)
@show minimum(b_sim), mean(b_sim), maximum(b_sim)

println("Sim synchronized ER")
b_sync_sim = @benchmark simulate(sync_ER)
@btime simulate(sync_ER)
@show minimum(b_sync_sim), mean(b_sync_sim), maximum(b_sync_sim)

println("Sim volatile synchronized ER")
b_vol_sync_sim = @benchmark volatile_simulate(sync_ER)
@btime volatile_simulate(sync_ER)
@show minimum(b_vol_sync_sim), mean(b_vol_sync_sim), maximum(b_vol_sync_sim)

Profile.clear_malloc_data() 
volatile_simulate(sync_ER)

