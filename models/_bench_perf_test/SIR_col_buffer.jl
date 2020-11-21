
import StaticArrays: SVector, SMatrix, @SMatrix

d=3
k=2
dict_var = Dict("S" => 1, "I" => 2, "R" => 3)
dict_p = Dict("ki" => 1, "kr" => 2)
l_tr = ["R1","R2"]
p = [0.0012, 0.05]
x0 = [95, 5, 0]
t0 = 0.0
function SIR_col_buffer_f!(mat_x::Matrix{Int}, l_t::Vector{Float64}, l_tr::Vector{String}, idx::Int,
           xn::AbstractVector{Int}, tn::Float64, p::Vector{Float64})
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
 
    nu = @view l_nu[:,reaction] # macro for avoiding a copy
    for i = 1:3
        mat_x[i,idx] = xn[i]+nu[i]
    end
    l_t[idx] = tn + tau
    l_tr[idx] = "R$(reaction)"
end
isabsorbing_SIR_col_buffer(p::Vector{Float64}, xn::AbstractVector{Int}) = (p[1]*xn[1]*xn[2] + p[2]*xn[2]) === 0.0
g = ["I"]

SIR_col_buffer = ContinuousTimeModel(d,k,dict_var,dict_p,l_tr,p,x0,t0,SIR_col_buffer_f!,isabsorbing_SIR_col_buffer; g=g)
export SIR_col_buffer

