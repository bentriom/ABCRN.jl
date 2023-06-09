
using Profile
using Statistics
using BenchmarkTools
using BiochemNetABC
import LinearAlgebra: dot
import Distributions: Uniform

load_automaton("euclidean_distance_automaton")
load_automaton("euclidean_distance_automaton_2")
load_automaton("abc_euclidean_distance_automaton")
load_model("repressilator")
tb = 210.0
tml_obs = 0:10.0:210.0
set_param!(repressilator, [:α, :β, :n, :α0], [200.0, 2.0, 2.0, 0.0])
set_time_bound!(repressilator, tb)
y_obs = vectorize(simulate(repressilator), :P1, tml_obs)

println("Vectorize:")
b_vectorize = @benchmark begin 
    σ = simulate($(repressilator)) 
    euclidean_distance(σ, :P1, tml_obs, y_obs)
end
#=
@btime begin
    σ = simulate($(repressilator)) 
    euclidean_distance(σ, :P1, tml_obs, y_obs)
end
=#
@show minimum(b_vectorize), mean(b_vectorize), maximum(b_vectorize)

println("Automaton with 1 loc")
aut1 = create_euclidean_distance_automaton(repressilator, tml_obs, y_obs, :P1)
sync1 = repressilator * aut1
b_sim_aut1 = @benchmark (σ = simulate($(sync1)))
#@btime (σ = simulate($(sync1)))
@show minimum(b_sim_aut1), mean(b_sim_aut1), maximum(b_sim_aut1)
b_vol_sim_aut1 = @benchmark (σ = volatile_simulate($(sync1)))
#@btime (σ = volatile_simulate($(sync1)))
@show minimum(b_vol_sim_aut1), mean(b_vol_sim_aut1), maximum(b_vol_sim_aut1)

println("ABC reject automaton with 1 loc")
aut1_abc = create_abc_euclidean_distance_automaton(repressilator, tml_obs, y_obs, :P1)
aut1_abc.ϵ = Inf
sync1_abc = repressilator * aut1_abc
b_sim_aut1_abc = @benchmark (σ = simulate($(sync1_abc)))
#@btime (σ = simulate($(sync1_abc)))
@show minimum(b_sim_aut1_abc), mean(b_sim_aut1_abc), maximum(b_sim_aut1_abc)
b_vol_sim_aut1_abc = @benchmark (σ = volatile_simulate($(sync1_abc)))
#@btime (σ = volatile_simulate($(sync1_abc)))
@show minimum(b_vol_sim_aut1_abc), mean(b_vol_sim_aut1_abc), maximum(b_vol_sim_aut1_abc)

#=
println("Memory test")
Profile.clear_malloc_data()
σ = volatile_simulate(sync1)
exit()
=#

println("Automaton with nbr_obs loc")
aut2 = create_euclidean_distance_automaton_2(repressilator, tml_obs, y_obs, :P1)
sync2 = repressilator * aut2
b_sim_aut2 = @benchmark (σ = simulate($(sync2)))
#@btime (σ = simulate($(sync2)))
@show minimum(b_sim_aut2), mean(b_sim_aut2), maximum(b_sim_aut2)
b_vol_sim_aut2 = @benchmark (σ = volatile_simulate($(sync2)))
#@btime (σ = volatile_simulate($(sync2)))
@show minimum(b_vol_sim_aut2), mean(b_vol_sim_aut2), maximum(b_vol_sim_aut2)

