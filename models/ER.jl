
import StaticArrays: SVector, SMatrix, @SVector, @SMatrix

d=4
k=3
dict_var_ER = Dict(:E => 1, :S => 2, :ES => 3, :P => 4)
dict_p_ER = Dict(:k1 => 1, :k2 => 2, :k3 => 3)
l_tr_ER = [:R1,:R2,:R3]
p_ER = [1.0, 1.0, 1.0]
x0_ER = [100, 100, 0, 0]
t0_ER = 0.0
function ER_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{<:Transition},
               xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    @inbounds a1 = p[1] * xn[1] * xn[2]
    @inbounds a2 = p[2] * xn[3]
    @inbounds a3 = p[3] * xn[3]
    l_a = (a1, a2, a3)
    asum = sum(l_a)
    if asum == 0.0
        copyto!(xnplus1, xn)
        return nothing
    end
    nu_1 = (-1, -1, 1, 0)
    nu_2 = (1, 1, -1, 0)
    nu_3 = (1, 0, -1, 1) 
    l_nu = (nu_1, nu_2, nu_3)
    l_str_R = (:R1, :R2, :R3)

    u1 = rand()
    u2 = rand()
    tau = - log(u1) / asum
    b_inf = 0.0
    b_sup = a1
    reaction = 0
    for i = 1:3
        if b_inf < asum*u2 < b_sup
            reaction = i
            break
        end
        @inbounds b_inf += l_a[i]
        @inbounds b_sup += l_a[i+1]
    end
 
    nu = l_nu[reaction]
    for i = 1:4
        @inbounds xnplus1[i] = xn[i]+nu[i]
    end
    @inbounds l_t[1] = tn + tau
    @inbounds l_tr[1] = l_str_R[reaction]
end
isabsorbing_ER(p::Vector{Float64},xn::Vector{Int}) = 
    @inbounds(p[1]*xn[1]*xn[2] + (p[2]+p[3])*xn[3] === 0.0)
g_ER = [:P]

ER = ContinuousTimeModel(d,k,dict_var_ER,dict_p_ER,l_tr_ER,p_ER,x0_ER,t0_ER,ER_f!,isabsorbing_ER; g=g_ER,name="ER pkg")

function create_ER(new_p::Vector{Float64})
    ER_new = deepcopy(ER)
    @assert length(ER_new.p) == length(new_p)
    set_param!(ER_new, new_p)
    return ER_new
end

export ER, create_ER

