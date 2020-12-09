
import StaticArrays: SVector, SMatrix, @SVector, @SMatrix
import Distributions: Poisson, rand

d=1
k=1
dict_var_poisson = Dict(:N => 1)
dict_p_poisson = Dict(:λ => 1)
l_tr_poisson = [:R]
p_poisson = [5.0]
x0_poisson = [0]
t0_poisson = 0.0
function poisson_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{<:Transition},
                    xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    
    u1 = rand()
    tau = (-log(u1)/p[1])
    xnplus1[1] += 1
    l_t[1] = tn + tau
    l_tr[1] = :R
end
isabsorbing_poisson(p::Vector{Float64}, xn::Vector{Int}) = p[1] === 0.0
g_poisson = [:N]

poisson = ContinuousTimeModel(d,k,dict_var_poisson,dict_p_poisson,l_tr_poisson,
                                  p_poisson,x0_poisson,t0_poisson,
                                  poisson_f!,isabsorbing_poisson; 
                                  g=g_poisson, time_bound=1.0, name="Poisson process pkg")
function create_poisson(new_p::Vector{Float64})
    poisson_new = deepcopy(poisson)
    @assert length(poisson_new.p) == length(new_p)
    set_param!(poisson_new, new_p)
    return poisson_new
end

export poisson, create_poisson

