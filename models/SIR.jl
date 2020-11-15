
import StaticArrays: SVector, SMatrix, @SMatrix

State = SVector{3, Int}
Parameters = SVector{2, Real}

d=3
k=2
dict_var = Dict("S" => 1, "I" => 2, "R" => 3)
dict_p = Dict("ki" => 1, "kr" => 2)
l_tr = ["R1","R2"]
p = Parameters(0.0012, 0.05)
x0 = State(95, 5, 0)
t0 = 0.0
function f(xn::State, tn::Real, p::Parameters)
    a1 = p[1] * xn[1] * xn[2]
    a2 = p[2] * xn[2]
    l_a = SVector(a1, a2)
    asum = sum(l_a)
    # column-major order
    l_nu = @SMatrix [-1 0;
                     1 -1;
                     0 1]
    
    u1, u2 = rand(), rand()
    tau = - log(u1) / asum
    b_inf = 0.0
    b_sup = a1
    reaction = 0
    for i = 1:2 
        if b_inf < asum*u2 < b_sup
            reaction = i
            break
        end
        b_inf += l_a[i]
        b_sup += l_a[i+1]
    end
 
    nu = l_nu[:,reaction]
    xnplus1 = State(xn[1]+nu[1], xn[2]+nu[2], xn[3]+nu[3])
    tnplus1 = tn + tau
    transition = "R$(reaction)"

    return xnplus1, tnplus1, transition
end
is_absorbing_sir(p::Parameters,xn::State) = (p[1]*xn[1]*xn[2] + p[2]*xn[2]) == 0.0
g = SVector("I")

SIR = CTMC(d,k,dict_var,dict_p,l_tr,p,x0,t0,f,is_absorbing_sir; g=g)
export SIR

