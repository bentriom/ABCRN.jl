
using BiochemNetABC

load_plots()
load_model("ER")

pm_ER = ParametricModel(ER, (:k1, Uniform(0.0,100.0)), (:k2, Uniform(0.0,100.0)))

prior_p = [0.2, 40.0]

σ = simulate(pm_ER, prior_p)
plot(σ; filename = get_module_path() * "/test/simulation/res_pics/sim_pm_er_long.svg")

return true

