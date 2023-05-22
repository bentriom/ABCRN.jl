
using Distributed
@everywhere using ABCRN
using Dates
using Plots
using DelimitedFiles

@assert length(ARGS) > 0 "No arguments entered. Please specify at least one argument that specifies the experiment: R1, R2, R3; e.g. julia enzym_1d.jl R3."
@assert length(ARGS) <= 2 "Too much arguments"
exp = ARGS[1]
@assert exp in ["R1", "R2", "R3", "R5"] "Wrong name of experiment"
nbr_exp = parse(Int, exp[2])
date_run = Dates.format(now(), "Y-m-d_HH:MM:SS")
path_results = "./results_$(exp)_$(date_run)/"
if length(ARGS) == 2
    path_results = ARGS[2]
end
if !isdir(path_results) mkdir(path_results) end
exps_p_star_k3 = [90.0, 25.0, 10.0]
exps_prob_p_star_k3 = [1.0, 1.0, 0.9997]
#=
samples_abc_post = readdlm(path_results * "mat_p_end.csv", ',')
samples_weights = readdlm(path_results * "weights_end.csv", ',')[:,1]
=#

# Chemical reaction network model
load_model("ER")
observe_all!(ER)
# Choice of the automaton
load_automaton("automaton_F")
load_automaton("automaton_G")
dict_automata = Dict()
dict_automata["R1"] = create_automaton_F(ER, 50.0, 75.0, 0.025, 0.05, :P) 
dict_automata["R2"] = create_automaton_F(ER, 50.0, 75.0, 0.05, 0.075, :P) 
dict_automata["R3"] = create_automaton_F(ER, 25.0, 50.0, 0.05, 0.075, :P)
aut = dict_automata[exp]
# Synchronized model
sync_ER = aut * ER 
pm_sync_ER = ParametricModel(sync_ER, (:k3, Uniform(0.0, 100.0)))
nbr_pa = 1500

r = automaton_abc(pm_sync_ER; nbr_particles = nbr_pa, dir_results = path_results)
samples_abc_post = r.mat_p_end[1,:]
samples_weights = r.weights

# Histogram
histogram(samples_abc_post, weights = r.weights, normalize = :density)
savefig(path_results * "histogram.svg")

## Satisfaction function

# Optimal bandwidth with KernelEstimator
using KernelEstimator
observations = samples_abc_post
observations_scaled = observations ./ 100
d_exp_kernel_func = Dict("R1" => betakernel,
                         "R2" => betakernel,
                         "R3" => gaussiankernel)
kernel_func = d_exp_kernel_func[exp]
opt_bw = bwlscv(observations_scaled, betakernel)
@show opt_bw
pdf_estim_abc_post(x::Float64) = kerneldensity(observations, xeval = [x], lb = 0.0, ub = 100.0, kernel = betakernel, h = opt_bw)[1]

# Optimal bandwidth with BoundedKDE
using BoundedKDE
d_exp_kernel = Dict("R1" => "chen99",
                    "R2" => "chen99",
                    "R3" => "gaussian")
kernel_kde = d_exp_kernel[exp]
bw_lbound = [0.02, 0.02, 0.05]
bw_ubound = [1.2, 1.2, 1.2]
@show kernel_kde
estim_abc_post = UnivariateKDE(samples_abc_post; kernel = kernel_kde, weights = r.weights, lower_bound = 0.0, upper_bound = 100.0)
#lscv_bandwidth = select_bandwidth(estim_abc_post, bw_lbound[nbr_exp], bw_ubound[nbr_exp], 10; verbose = true)
lscv_bandwidth = minimize_lscv(estim_abc_post, bw_lbound[nbr_exp], bw_ubound[nbr_exp]; verbose = true)
@show lscv_bandwidth, lscv_bandwidth / asymptotic_bandwidth(estim_abc_post)
estim_abc_post = change_bandwidth(estim_abc_post, lscv_bandwidth)
#estim_abc_post = change_bandwidth(estim_abc_post, opt_bw)
pdf_estim_abc_post_pkg(x) = pdf(estim_abc_post, x)

# Estimation of the constant with a probability estimated by Statistical Model Checking
constant = exps_prob_p_star_k3[nbr_exp] / pdf_estim_abc_post(exps_p_star_k3[nbr_exp])
constant_pkg = exps_prob_p_star_k3[nbr_exp] / pdf_estim_abc_post_pkg(exps_p_star_k3[nbr_exp])

# Plot of satisfaction probability function
prob_func(x) = pdf_estim_abc_post(x) * constant
prob_func_pkg(x) = pdf_estim_abc_post_pkg(x) * constant_pkg
xaxis = 0:0.1:100
plot(xaxis, prob_func.(xaxis), label = "ABC spf KernelEstimator", dpi = 480)
plot!(xaxis, prob_func_pkg.(xaxis), label = "ABC spf BoundedKDE", dpi = 480)
y_MC = readdlm("/home/moud/plot_R1-3/estim_MC/$(exp)/satisfaction_func.csv", ',')[:,1]
inf_x, sup_x = 0.0, 100.0
x_MC = inf_x:((sup_x-inf_x)/(length(y_MC)-1)):sup_x
plot!(x_MC, y_MC, label = "MC spf")
savefig(path_results * "satisfaction_prob_function.svg")

