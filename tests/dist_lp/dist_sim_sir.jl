
using MarkovProcesses
load_model("SIR")
time_bound = 100.0
SIR.time_bound = time_bound

test_all = true

for p = 1:2
    nb_sim = 10
    for i = 1:nb_sim
        σ1 = simulate(SIR)
        σ2 = simulate(SIR)
        d = dist_lp(σ1, σ2, "I"; p = p)
        d2 = dist_lp(σ1, σ2; p = p)

        f_x(t::Float64) = MarkovProcesses._f_step(σ1["I"], σ1["times"], t)
        f_y(t::Float64) = MarkovProcesses._f_step(σ2["I"], σ2["times"], t)
        diff_f(t) = abs(f_x(t) - f_y(t))^p
        int_riemann = MarkovProcesses._riemann_sum(diff_f, SIR.t0, SIR.time_bound, 1E-3)
        res_int_riemann = int_riemann^(1/p)

        test = isapprox(d, res_int_riemann; atol = 1E-1)
        test2 = isapprox(d2, res_int_riemann; atol = 1E-1)
        #@show d, res_int_riemann
        #@show test
        
        global test_all = test_all && test && test2
    end
end

return test_all

