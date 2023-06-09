
using BiochemNetABC
import QuadGK: quadgk
load_model("SIR")

# Case 1 : 6.4

test_all = true

let x_obs, y_obs, t_x, t_y, σ1, σ2
    x_obs = [2, 4, 3, 3]
    t_x = [0.0, 0.5, 1.2, 2.2]
    t_y_bis = [0.0, 0.5, 1.2, 2.2]
    values = [zeros(length(x_obs))]
    values[1] = x_obs
    l_tr = Vector{Nothing}(nothing, length(x_obs))
    σ1 = Trajectory(SIR, values, t_x, l_tr)

    y_obs = [6, 6]
    t_y = [0.0, 2.2]
    t_x_bis = [0.0, 2.2]
    values = [zeros(length(y_obs))]
    values[1] = y_obs
    l_tr = Vector{Nothing}(nothing, length(y_obs))
    σ2 = Trajectory(SIR, values, t_y, l_tr)

    test_1 = dist_lp(x_obs, t_x, y_obs, t_y; p=1) == 6.4
    test_1 = test_1 && dist_lp(σ1, σ2; p=1) == 6.4

    f_x(t::Real) = BiochemNetABC._f_step(x_obs, t_x, t)
    f_y(t::Real) = BiochemNetABC._f_step(y_obs, t_y, t)
    diff_f(t) = abs(f_x(t) - f_y(t))
    int, err = quadgk(diff_f, 0.0, 2.2)

    test_2 = isapprox(6.4, int; atol=err)

    # Case 1 bis : inverse of case 1 

    test_1_bis = dist_lp(y_obs, t_y, x_obs, t_x; p=1) == 6.4 
    test_1_bis = test_1_bis && dist_lp(σ2, σ1; p=1) == 6.4

    global test_all = test_1 && test_1_bis && test_2
end

return test_all

