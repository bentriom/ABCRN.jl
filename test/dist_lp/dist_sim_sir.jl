
using ABCRN
load_model("SIR")
time_bound = 100.0
SIR.time_bound = time_bound

test_all = true

for p = 1:2
    local nb_sim = 10
    for i = 1:nb_sim
        let σ1, σ2, test, test2
            σ1 = simulate(SIR)
            σ2 = simulate(SIR)
            d = dist_lp(σ1, σ2, :I; p = p)
            d2 = dist_lp(σ1, σ2; p = p)

            f_x(t::Float64) = ABCRN._f_step(σ1[:I], times(σ1), t)
            f_y(t::Float64) = ABCRN._f_step(σ2[:I], times(σ2), t)
            diff_f(t) = abs(f_x(t) - f_y(t))^p
            int_riemann = ABCRN._riemann_sum(diff_f, SIR.t0, SIR.time_bound, 1E-3)
            res_int_riemann = int_riemann^(1/p)

            test = isapprox(d, res_int_riemann; atol = 1E-1)
            test2 = isapprox(d2, res_int_riemann; atol = 1E-1)
            #@show d, res_int_riemann
            #@show test

            global test_all = test_all && test && test2
        end
    end
end

return test_all

