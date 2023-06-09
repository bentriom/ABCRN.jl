
using BenchmarkTools
using BiochemNetABC

load_model("ER")
long_p = [0.2, 40.0, 1.0] 
xn = ER.x0
tn = ER.t0
vec_x = zeros(Int, ER.dim_state)
l_t = [0.0]
l_tr = Transition[nothing]
b1 = @benchmark getfield(ER, :f!)(vec_x, l_t, l_tr, xn, tn, ER.p)
@show minimum(b1), mean(b1), maximum(b1)

