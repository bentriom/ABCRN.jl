
import StaticArrays: SVector, SMatrix, @SVector, @SMatrix

d=3
k=2
dict_var = Dict("S" => 1, "I" => 2, "R" => 3)
dict_p = Dict("ki" => 1, "kr" => 2)
l_tr_SIR = ["R1","R2"]
p_SIR = [0.0012, 0.05]
x0_SIR = [95, 5, 0]
t0_SIR = 0.0
function SIR_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{Union{Nothing,String}},
                xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    a1 = p[1] * xn[1] * xn[2]
    a2 = p[2] * xn[2]
    l_a = (a1, a2)
    asum = sum(l_a)
    if asum == 0.0
        copyto!(xnplus1, xn)
        return nothing
    end
    # column-major order
    nu_1 = (-1, 1, 0)
    nu_2 = (0, -1, 1)
    l_nu = (nu_1, nu_2)
    l_str_R = ("R1","R2")

    u1 = rand()
    u2 = rand()
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
 
    nu = l_nu[reaction]
    for i = 1:3
        xnplus1[i] = xn[i]+nu[i]
    end
    l_t[1] = tn + tau
    l_tr[1] = l_str_R[reaction]
end
isabsorbing_SIR(p::Vector{Float64}, xn::Vector{Int}) = (p[1]*xn[1]*xn[2] + p[2]*xn[2]) === 0.0
g_SIR = ["I"]

SIR = ContinuousTimeModel(d,k,dict_var,dict_p,l_tr_SIR,p_SIR,x0_SIR,t0_SIR,SIR_f!,isabsorbing_SIR; g=g_SIR)

function create_SIR(new_p::Vector{Float64})
    SIR_new = deepcopy(SIR)
    @assert length(SIR_new.p) == length(new_p)
    set_param!(SIR_new, new_p)
    return SIR_new
end

export SIR, create_SIR

