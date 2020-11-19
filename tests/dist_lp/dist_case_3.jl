
using MarkovProcesses
import QuadGK: quadgk
load_model("SIR")

test_all = true
for p = 1:2
    x_obs =[5, 6, 5, 4, 3, 2, 1, 1]
    t_x = [0.0, 3.10807, 4.29827, 4.40704, 5.67024, 7.1299, 11.2763, 20.0]
    values = zeros(length(x_obs), 1)
    values[:,1] = x_obs
    l_tr = Vector{Nothing}(nothing, length(x_obs))
    σ1 = Trajectory(SIR, values, t_x, l_tr)

    y_obs =[5, 4, 5, 4, 3, 4, 3, 2, 3, 2, 1, 2, 3, 4, 3, 4, 4]
    t_y = [0.0, 0.334082, 1.21012, 1.40991, 1.58866, 2.45879, 2.94545, 4.66746, 5.44723, 5.88066, 7.25626, 11.4036, 13.8373, 17.1363, 17.8193, 18.7613, 20.0]
    values = zeros(length(y_obs), 1)
    values[:,1] = y_obs
    l_tr = Vector{Nothing}(nothing, length(y_obs))
    σ2 = Trajectory(SIR, values, t_y, l_tr)

    f_x(t::Real) = MarkovProcesses._f_step(x_obs, t_x, t)
    f_y(t::Real) = MarkovProcesses._f_step(y_obs, t_y, t)
    diff_f(t) = abs(f_x(t) - f_y(t))^p
    int, err = quadgk(diff_f, 0.0, 20.0, rtol=1e-10)
    int_riemann = MarkovProcesses._riemann_sum(diff_f, 0.0, 20.0, 1E-5)
    int_riemann = int_riemann^(1/p)

    res1 = dist_lp(@view(x_obs[:]), t_x, @view(y_obs[:]), t_y; p=p)
    res2 = dist_lp(σ1,σ2; p=p)
    res1_bis = dist_lp(@view(y_obs[:]), t_y, @view(x_obs[:]), t_x; p=p)
    res2_bis = dist_lp(σ2,σ1; p=p)
    test_1 = isapprox(res1, int_riemann; atol = 1E-3) 
    test_1 = test_1 && isapprox(res2, int_riemann; atol = 1E-3)
    test_1_bis = isapprox(res1_bis, int_riemann; atol = 1E-3) 
    test_1_bis = test_1_bis && isapprox(res2_bis, int_riemann; atol = 1E-3)

    test_2 = res1 == res2 == res1_bis == res2_bis
    
    global test_all = test_all && test_1 && test_1_bis && test_2
end

return test_all

