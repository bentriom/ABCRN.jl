
using Plots
using BiochemNetABC
import LinearAlgebra: dot

# ER model
load_model("ER")
set_param!(ER, [:k1, :k2, :k3], [0.2, 40.0, 1.0])

# Observations
timeline = 0.0:0.1:2.0
observations = [vectorize(simulate(ER), :P, timeline)]
epsilon = 0.1 * sqrt(dot(observations[1], observations[1]))

load_automaton("abc_euclidean_distance_automaton")
aut = create_abc_euclidean_distance_automaton(ER, timeline, observations[1], :P)
sync_ER = ER * aut
pm_sync_ER = ParametricModel(sync_ER, (:k1, Uniform(0.0, 100.0)), (:k2, Uniform(0.0, 100.0)))

path_res = "abc_eucl_aut/"
res_abc = automaton_abc(pm_sync_ER; nbr_particles = 100, tolerance = epsilon)

#=
samples_abc_post = res_abc.mat_p_end
samples_weights = res_abc.weights
histogram2d(samples_abc_post[1,:], samples_abc_post[2,:], weights = samples_weights, normalize = :density)
savefig(path_res * "/histogram.svg")
=#

return true

