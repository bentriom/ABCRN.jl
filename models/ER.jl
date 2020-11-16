
import StaticArrays: SVector, SMatrix, @SMatrix

d=4
k=3
dict_var = Dict("E" => 1, "S" => 2, "ES" => 3, "P" => 4)
dict_p = Dict("k1" => 1, "k2" => 2, "k3" => 3)
l_tr = ["R1","R2","R3"]
p = [1.0, 1.0, 1.0]
x0 = [100, 100, 0, 0]
t0 = 0.0

function f!(xnplus1::Vector{Int}, tnplus1::Vector{Float64}, tr::Vector{String}, 
           xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    a1 = p[1] * xn[1] * xn[2]
    a2 = p[2] * xn[3]
    a3 = p[3] * xn[3]
    l_a = SVector(a1, a2, a3)
    asum = sum(l_a)
    l_nu = @SMatrix [-1 1 1;
                     -1 1 0;
                     1 -1 -1;
                     0 0 1]
    u1, u2 = rand(), rand()
    tau = - log(u1) / asum
    b_inf = 0.0
    b_sup = a1
    reaction = 0
    for i = 1:3
        if b_inf < asum*u2 < b_sup
            reaction = i
            break
        end
        b_inf += l_a[i]
        b_sup += l_a[i+1]
    end
 
    nu = @view l_nu[:,reaction]
    xnplus1[1] = xn[1]+nu[1] 
    xnplus1[2] = xn[2]+nu[2]
    xnplus1[3] = xn[3]+nu[3]
    xnplus1[4] = xn[4]+nu[4]
    tnplus1[1] = tn + tau
    tr[1] = "R$(reaction)"
end
is_absorbing_er(p::Vector{Float64},xn::Vector{Int}) = 
    (p[1]*xn[1]*xn[2] + (p[2]+p[3])*xn[3]) === 0.0
g = ["P"]

ER = CTMC(d,k,dict_var,dict_p,l_tr,p,x0,t0,f!,is_absorbing_er; g=g)
export ER

