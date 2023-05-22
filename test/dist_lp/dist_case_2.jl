
using ABCRN
import QuadGK: quadgk
load_model("SIR")

test_all = true
p=2
for p = 1:2
    let x_obs, y_obs, t_x, t_y, σ1, σ2, test_1, test_1_bis, test_2
        res = (4-1)^p * 0.5 + (4-2)^p * 0.1 + 0 + (5-3)^p * 0.1 + (5-1)^p * (1.4 - 0.8) + (3-1)^p * (3.0-1.4)
        res = res^(1/p)

        y_obs = [1, 2, 5, 3, 3]
        t_y = [0.0, 0.5, 0.7, 1.4, 3.0]
        values = [zeros(length(y_obs))]
        values[1] = y_obs
        l_tr = Vector{Nothing}(nothing, length(y_obs))
        σ2 = Trajectory(SIR, values, t_y, l_tr)

        x_obs = [4, 2, 3, 1, 1]
        t_x = [0.0, 0.6, 0.7, 0.8, 3.0]
        values = [zeros(length(x_obs))]
        values[1] = x_obs
        l_tr = Vector{Nothing}(nothing, length(x_obs))
        σ1 = Trajectory(SIR, values, t_x, l_tr)

        test_1 = isapprox(dist_lp(x_obs, t_x, y_obs, t_y; p=p), res; atol = 1E-10) 
        test_1 = test_1 && isapprox(dist_lp(σ1,σ2; p=p), res; atol = 1E-10)
        test_1_bis = isapprox(dist_lp(y_obs, t_y, x_obs, t_x; p=p), res; atol = 1E-10) 
        test_1_bis = test_1_bis && isapprox(dist_lp(σ2,σ1; p=p), res; atol = 1E-10)

        f_x(t::Real) = ABCRN._f_step(x_obs, t_x, t)
        f_y(t::Real) = ABCRN._f_step(y_obs, t_y, t)
        diff_f(t) = abs(f_x(t) - f_y(t))^p
        int, err = quadgk(diff_f, 0.0, 3.0)
        res_int = int^(1/p)

        test_2 = isapprox(res, res_int; atol = err)
        global test_all = test_all && test_1 && test_1_bis && test_2
    end
end

return test_all

