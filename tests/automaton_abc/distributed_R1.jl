
using MarkovProcesses
using Distributed
addprocs(2)
module_path = get_module_path()
@everywhere module_path = $module_path
@everywhere push!(LOAD_PATH, "$(module_path)/core")
@everywhere using MarkovProcesses
#=
@everywhere begin
    path_module = $(path_module)
    push!(LOAD_PATH, path_module)
    using MarkovProcesses
    load_model("ER")
    load_automaton("automaton_F")
end
=#
load_model("ER")
load_automaton("automaton_F")
A_F_R1 = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P)
sync_ER = A_F_R1*ER 
pm_sync_ER = ParametricModel(sync_ER, (:k3, Uniform(0.0, 100.0)))
nbr_pa = 404

r = automaton_abc(pm_sync_ER; nbr_particles = nbr_pa)

rmprocs(workers())

test = size(r.mat_p_end)[1] == pm_sync_ER.df &&
       size(r.mat_p_end)[2] == nbr_pa &&
       length(r.vec_dist) == nbr_pa &&
       length(r.weights) == nbr_pa &&
       r.epsilon == 0.0

return test

#histogram(r.mat_p_end', weights = r.weights, normalize = :density)
#savefig("R1_hist.svg")

