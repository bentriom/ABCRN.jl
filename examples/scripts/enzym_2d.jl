
using DelimitedFiles
using Distributed
@everywhere using MarkovProcesses
using Dates
using Plots

@assert length(ARGS) > 0 "No arguments entered. Please specify at least one argument that specifies the experiment: R4, R5, R6; e.g. julia enzym_2d.jl R5."
@assert length(ARGS) <= 2 "Too much arguments"
exp = ARGS[1]
@assert exp in ["R4", "R5", "R6"] "Wrong name of experiment"
nbr_exp = parse(Int, exp[2]) - 3
date_run = Dates.format(now(), "Y-m-d_HH:MM:SS")
path_results = "./results_$(exp)/"
if length(ARGS) == 2
    path_results = ARGS[2]
end
if !isdir(path_results) mkdir(path_results) end
exps_p_star_k1_k2 = [[0.01, 60.0], [0.25, 40.0], [0.50, 40.0]]
exps_prob_p_star_k1_k2 = [0.8634, 1.0, 0.363553]

samples_abc_post = readdlm(path_results * "mat_p_end.csv", ',')
samples_weights = readdlm(path_results * "weights_end.csv", ',')[:,1]
#=
# Chemical reaction network model
load_model("ER")
observe_all!(ER)
# Choice of the automaton
load_automaton("automaton_F")
load_automaton("automaton_G")
load_automaton("automaton_G_and_F")
dict_automata = Dict()
dict_automata["R4"] =  create_automaton_F(ER, 5.0, 15.0, 8.0, 10.0, :P)
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
x3, x4, t3, t4 = 30.0, 100.0, 0.8, 0.9
dict_automata["R5"] = create_automaton_G(ER, x1, x2, t1, t2, :E) 
dict_automata["R6"] = create_automaton_G_and_F(ER, x1, x2, t1, t2, :E, x3, x4, t3, t4, :P)
aut = dict_automata[exp]
# Synchronized model
sync_ER = aut * ER 
pm_sync_ER = ParametricModel(sync_ER, (:k1, Uniform(0.0, 100.0)), (:k2, Uniform(0.0, 100.0)))
nbr_pa = 1000

r = automaton_abc(pm_sync_ER; nbr_particles = nbr_pa, dir_results = path_results)
samples_abc_post = r.mat_p_end
samples_weights = r.weights

# Histogram
histogram2d(samples_abc_post[1,:], samples_abc_post[2,:], weights = samples_weights, normalize = :density)
savefig(path_results * "histogram.svg")
=#
## Satisfaction function

# Optimal bandwidth with BoundedKDE
using BoundedKDE
d_exp_kernel = Dict("R4" => "chen99",
                    "R5" => "chen99",
                    "R6" => "gaussian")
kernel_kde = d_exp_kernel[exp]
#bw_lbound = fill.([0.01, 0.05, 0.01], 2)
#bw_ubound = fill.([1.0, 0.5, 1.0], 2)
bw_lbound = [fill(0.1, 2), [0.5, 0.5], fill(0.1, 2)]
bw_ubound = [fill(1.0, 2), [1.5, 1.5], fill(1.0, 2)]
coeff_bw_R5 = [0.025704814875365627, 0.07660069455807658] # -0.0166403765
estim_abc_post = MultivariateKDE(samples_abc_post; kernel = kernel_kde, weights = samples_weights, lower_bound = [0.0, 0.0], upper_bound = [2.0, 100.0])
bw_ref = asymptotic_bandwidth(estim_abc_post)
@show kernel_kde
@show lscv(estim_abc_post, bw_ref .* coeff_bw_R5)
lscv_bandwidth = select_bandwidth(estim_abc_post, bw_lbound[nbr_exp], bw_ubound[nbr_exp], fill(5, 2); maxevals_int = typemax(Int))
#=
lscv_bandwidth = minimize_lscv(estim_abc_post, bw_lbound[nbr_exp], bw_ubound[nbr_exp]; 
                               coeff_init_bw = coeff_bw_R5, rt = 0.5, f_tol = 1E-5, x_tol = 1E-4,
                               verbosity = 2, maxevals_int = 200)
=#
#lscv_bandwidth = [0.004308901764294055, 0.07706748767996473] 
@show lscv_bandwidth, lscv_bandwidth / asymptotic_bandwidth(estim_abc_post)
estim_abc_post = change_bandwidth(estim_abc_post, lscv_bandwidth)
#estim_abc_post = change_bandwidth(estim_abc_post, opt_bw)
pdf_estim_abc_post(x) = pdf(estim_abc_post, x)

# Estimation of the constant with a probability estimated by Statistical Model Checking
constant = exps_prob_p_star_k1_k2[nbr_exp] / pdf_estim_abc_post(exps_p_star_k1_k2[nbr_exp])

# Plot of satisfaction probability function
prob_func(x,y) = pdf_estim_abc_post([x,y]) * constant
@show prob_func(exps_p_star_k1_k2[nbr_exp]...), exps_prob_p_star_k1_k2[nbr_exp]
xaxis = 0:0.1:1.5
yaxis = 0:1.0:100.0
p = plot(title = "Multivariate KDE", dpi = 480, background_color_legend = :transparent)
plot!(p, xaxis, yaxis, prob_func, st = :surface, c = :coolwarm, camera = (30, 45))
savefig(path_results * "estim_abc_satisfaction_prob_function.svg")
x_MC = readdlm("/home/moud/results_last_automata/estim_satisfaction_func_MC/$(exp)/grid_X.csv", ',')
y_MC = readdlm("/home/moud/results_last_automata/estim_satisfaction_func_MC/$(exp)/grid_Y.csv", ',')
z_MC = readdlm("/home/moud/results_last_automata/estim_satisfaction_func_MC/$(exp)/satisfaction_func.csv", ',')
p = plot(title = "ABC MC", dpi = 480, background_color_legend = :transparent)
plot!(p, [x_MC...], [y_MC...], [z_MC...], st = :surface, c = :coolwarm, camera = (30, 45), label = "MC spf")
savefig(path_results * "estim_MC_satisfaction_prob_function.svg")

