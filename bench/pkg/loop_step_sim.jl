
using StaticArrays
using BenchmarkTools
using ABCRN
Transition = Union{String,Nothing}

load_model("ER")
long_p = [0.2, 40.0, 1.0] 
ER.p = long_p
xn = view(reshape(ER.x0, 1, ER.d), 1, :) # View for type stability
xn = ER.x0 
tn = ER.t0
mat_x = zeros(Int, 1, ER.d)
vec_x = zeros(Int, ER.d)
l_t = Float64[0.0]
l_tr = Transition[nothing]

function _update_values!(values::Vector{Vector{Int}}, times::Vector{Float64}, transitions::Vector{Transition},
                         xn::Vector{Int}, tn::Float64, tr_n::Transition, idx::Int)
    for k = eachindex(values) values[k][idx] = xn[k] end
    times[idx] = tn
    transitions[idx] = tr_n
end

function ER_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{Union{Nothing,String}},
               xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    a1 = p[1] * xn[1] * xn[2]
    a2 = p[2] * xn[3]
    a3 = p[3] * xn[3]
    l_a = SVector(a1, a2, a3)
    asum = sum(l_a)
    nu_1 = SVector(-1, -1, 1, 0)
    nu_2 = SVector(1, 1, -1, 0)
    nu_3 = SVector(1, 0, -1, 1) 
    l_nu = SVector(nu_1, nu_2, nu_3)
    l_str_R = SVector(:R1, :R2, :R3)

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
        b_inf += l_a[i]
        b_sup += l_a[i+1]
    end
 
    nu = l_nu[reaction]
    for i = 1:4
        xnplus1[i] = xn[i]+nu[i]
    end
    l_t[1] = tn + tau
    l_tr[1] = l_str_R[reaction]
end
ER_isabsorbing(p::Vector{Float64},xn::Vector{Int}) = 
    (p[1]*xn[1]*xn[2] + (p[2]+p[3])*xn[3]) === 0.0

function model_loop_f(m::ContinuousTimeModel, nb_steps::Int)
    # Alloc 
    #mat_full_values = zeros(Int, m.dim_state, nb_steps)
    full_values = Vector{Vector{Int}}(undef, m.dim_state)
    for i = eachindex(full_values) full_values[i] = zeros(Int, nb_steps) end
    times = zeros(Float64, nb_steps)
    transitions = Vector{Transition}(undef, nb_steps)
    n = 1
    #xn = view(reshape(m.x0, 1, m.dim_state), 1, :) # View for type stability
    loc_xn = m.x0 # View for type stability
    loc_tn = m.t0 
    loc_vec_x = zeros(Int, m.dim_state)
    loc_l_t = Float64[0.0]
    loc_l_tr = Transition[nothing]
    time_bound = m.time_bound
    for i = 2:nb_steps
        getfield(Main, ER.f!)(loc_vec_x, loc_l_t, loc_l_tr, loc_xn, loc_tn, m.p)
        tn = l_t[1]
        if tn > m.time_bound 
            break
        end
        n += 1
        xn = loc_vec_x
        # Updating value
        _update_values!(full_values, times, transitions, xn, tn, l_tr[1], i)
        isabsorbing = getfield(Main, m.isabsorbing)(m.p,xn)
        if isabsorbing 
            break
        end
        i += 1
    end
end

function loop_f(x0::Vector{Int}, t0::Float64, p::Vector{Float64}, d::Int, time_bound::Float64, nb_steps::Int)
    # Alloc 
    #mat_full_values = zeros(Int, m.dim_state, nb_steps)
    full_values = Vector{Vector{Int}}(undef, 4)
    for i = eachindex(full_values) full_values[i] = zeros(Int, nb_steps) end
    times = zeros(Float64, nb_steps)
    transitions = Vector{Transition}(undef, nb_steps)
    n = 1
    #xn = view(reshape(m.x0, 1, m.dim_state), 1, :) # View for type stability
    loc_xn = x0 # View for type stability
    loc_tn = t0 
    loc_vec_x = zeros(Int, d)
    loc_l_t = Float64[0.0]
    loc_l_tr = Transition[nothing]
    for i = 2:nb_steps
        getfield(Main, ER.f!)(loc_vec_x, loc_l_t, loc_l_tr, loc_xn, loc_tn, p)
        tn = l_t[1]
        if tn > time_bound 
            break
        end
        n += 1
        xn = loc_vec_x
        # Updating value
        _update_values!(full_values, times, transitions, xn, tn, l_tr[1], i)
        isabsorbing = ER_isabsorbing(p,xn)
        if isabsorbing 
            break
        end
        i += 1
    end
end

nb = 6000
b_step = @benchmark getfield(Main, ER.f!)(vec_x, l_t, l_tr, xn, tn, ER.p)
b_step_loc = @benchmark getfield(Main, ER.f!)(vec_x, l_t, l_tr, xn, tn, ER.p)
b_loop = @benchmark model_loop_f(ER, nb)
b_loop_loc = @benchmark loop_f(ER.x0, ER.t0, ER.p, ER.d, ER.time_bound, nb)

