
using Statistics
using BenchmarkTools
using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

load_automaton("euclidean_distance_automaton")
load_model("repressilator")
tb = 210.0
tml_obs = 0:10.0:210.0
set_param!(repressilator, [:α, :β, :n, :α0], [200.0, 2.0, 2.0, 0.0])
set_time_bound!(repressilator, tb)

b_sim = @benchmark σ = simulate($(repressilator)) 
@btime σ = simulate($(repressilator)) 
@show mean(b_sim).time, mean(b_sim).memory

