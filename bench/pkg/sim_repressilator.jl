
using Statistics
using BenchmarkTools
using BiochemNetABC
import LinearAlgebra: dot
import Distributions: Uniform

load_model("repressilator")
observe_all!(repressilator)
N_periods, ref_mean_tp = 8, 20.0
set_param!(repressilator, [:α, :β, :n, :α0], [200.0, 2.0, 2.0, 0.0])
set_time_bound!(repressilator, (N_periods)*ref_mean_tp)

b_sim = @benchmark σ = simulate($(repressilator)) 
@show mean(b_sim)
@show mean(b_sim).time, mean(b_sim).memory

load_automaton("period_automaton")
L, H = 20.0, 100.0
A_per = create_period_automaton(repressilator, L, H, N_periods, :P1; ref_mean_tp = ref_mean_tp, error_func = :max_mean_var_relative_error)
sync_repressilator = repressilator * A_per

b_period = @benchmark σ = simulate($(sync_repressilator)) 
@show mean(b_period)
@show mean(b_period).time, mean(b_period).memory

b_vol_period = @benchmark σ = volatile_simulate($(sync_repressilator)) 
@show mean(b_vol_period)
@show mean(b_vol_period).time, mean(b_vol_period).memory

