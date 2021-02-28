
using Profile
using Statistics
using BenchmarkTools
@everywhere using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

load_automaton("euclidean_distance_automaton")
load_automaton("euclidean_distance_automaton_2")
load_model("SIR")
tb = 7.0*12
tml_obs = 0:7:tb
set_time_bound!(SIR, tb)
y_obs = vectorize(simulate(SIR), :I, tml_obs)

println("Vectorize:")
b_vectorize = @benchmark (σ = simulate($(SIR)); euclidean_distance(σ, :I, tml_obs, y_obs)) 
@btime (σ = simulate($(SIR)); euclidean_distance(σ, :I, tml_obs, y_obs)) 
@show minimum(b_vectorize), mean(b_vectorize), maximum(b_vectorize)

println("Automaton with 1 loc")
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
sync1 = SIR * aut1
b_sim_aut1 = @benchmark (σ = simulate($(sync1)))
@btime (σ = simulate($(sync1)))
@show minimum(b_sim_aut1), mean(b_sim_aut1), maximum(b_sim_aut1)
b_vol_sim_aut1 = @benchmark (σ = volatile_simulate($(sync1)))
@btime (σ = volatile_simulate($(sync1)))
@show minimum(b_vol_sim_aut1), mean(b_vol_sim_aut1), maximum(b_vol_sim_aut1)

#=
println("Memory test")
Profile.clear_malloc_data()
σ = volatile_simulate(sync1)
=#

println("Automaton with nbr_obs loc")
aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs, :I)
sync2 = SIR * aut2
b_sim_aut2 = @benchmark (σ = simulate($(sync2)))
@btime (σ = simulate($(sync2)))
@show minimum(b_sim_aut2), mean(b_sim_aut2), maximum(b_sim_aut2)
b_vol_sim_aut2 = @benchmark (σ = volatile_simulate($(sync2)))
@btime (σ = volatile_simulate($(sync2)))
@show minimum(b_vol_sim_aut2), mean(b_vol_sim_aut2), maximum(b_vol_sim_aut2)

