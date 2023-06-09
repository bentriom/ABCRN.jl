
using BiochemNetABC
load_model("ER")

k1 = ER[:k1]
dist_mv_unif = product_distribution(Uniform.([2.5,6.0], [3.5,7.0]))
pm = ParametricModel(ER, [:k3,:k2], dist_mv_unif)
draw_model!(pm)
test1 = 2.5 <= ER[:k3] <= 3.5 && 6.0 <= ER[:k2] <= 7.0 && pm.df == 2

p_drawn = copy(ER.p)
p_drawn_bis = copy((pm.m).p)
draw_model!(pm)
test2 = 2.5 <= ER[:k3] <= 3.5 && 6.0 <= ER[:k2] <= 7.0 && pm.df == 2 && p_drawn == p_drawn_bis && ER.p != p_drawn

vec_p = zeros(2)                       
draw!(vec_p, pm)
test3 = 2.5 <= vec_p[1] <= 3.5 && 6.0 <= vec_p[2] <= 7.0 && pm.df == 2

fill!(pm, vec_p)
test4 = ER.p[pm._param_idx] == vec_p && ER.p[[3,2]] == vec_p && ER.p[1] == k1 

return test1 && test2 && test3 && test4

