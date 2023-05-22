
using Profile
using Statistics
using BenchmarkTools
using ABCRN

ER = @network_model begin
    R1: (E+S => ES, k1*E*S)
    R2: (ES => E+S, k2*ES)
    R3: (ES => E+P, k3*ES)
end "ER macro"
long_p = [0.2, 40.0, 1.0] 
xn = ER.x0
tn = ER.t0
vec_x = zeros(Int, ER.dim_state)
l_t = [0.0]
l_tr = Transition[nothing]

b_ER = @benchmark $(getfield(Main, ER.f!))(vec_x, l_t, l_tr, xn, tn, long_p)
@btime $(getfield(Main, ER.f!))(vec_x, l_t, l_tr, xn, tn, long_p)
@show minimum(b_ER), mean(b_ER), maximum(b_ER)

load_model("repressilator")
xn = repressilator.x0
tn = repressilator.t0
vec_x = zeros(Int, repressilator.dim_state)
l_t = [0.0]
l_tr = Transition[nothing]
b_repressilator = @benchmark $(getfield(Main, repressilator.f!))(vec_x, l_t, l_tr, xn, tn, repressilator.p)
@btime $(getfield(Main, repressilator.f!))(vec_x, l_t, l_tr, xn, tn, repressilator.p)
@show minimum(b_repressilator), mean(b_repressilator), maximum(b_repressilator)

println("Memory test")
Profile.clear_malloc_data()
@btime $(getfield(Main, repressilator.f!))(vec_x, l_t, l_tr, xn, tn, repressilator.p)

