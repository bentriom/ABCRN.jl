
d=3
k=7
dict_var_square = Dict(:A => 1, :HIGH => 2, :LOW => 3)
dict_params_square = Dict(:p1 => 1, :w1 => 2, :p2 => 3, :w2 => 4, :p3 => 5, :w3 => 6, :d => 7)
l_tr_square = [:Tu, :t1, :t2, :t3]
p_square = [1, 0, 0, 0, 0, 0, 5.0]
x0_square = [1, 0, 1]
t0_square = 0.0
function square_wave_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{Transition},
                        xn::Vector{Int}, tn::Float64, p::Vector{Float64})
    if p[1] == 0 && p[3] == 0 && p[5] == 0
        copyto!(xnplus1, xn)
        return nothing
    end
    # If in state LOW
    if xn[3] == 1
        xnplus1[1] = 10
        xnplus1[2] = 1
        xnplus1[3] = 0
        l_t[1] += p[7]
        l_tr[1] = :Tu
    # If in state HIGH
    else
        possible_transitions = (:t1, :t2, :t3)
        l_priority = (p[1], p[3], p[5])
        l_weights = (p[2], p[4], p[6])
        max_priority = maximum(l_priority)
        l_idx = zeros(Int, 0)
        for i = eachindex(l_priority)
            if l_priority[i] >= max_priority
                push!(l_idx, i)
            end
        end
        if length(l_idx) == 1
            transition = l_idx[1]
        else
            wsum = sum(l_weights[l_idx])
            b_inf = 0.0
            b_sup = l_weights[l_idx[1]]
            transition = 0
            u = rand()
            for i in l_idx
                if b_inf < wsum*u < b_sup
                    transition = i
                    break
                end
                @inbounds b_inf += l_weights[i]
                @inbounds b_sup += l_weights[i]
            end
        end
        xnplus1[1] = 1
        xnplus1[2] = 0
        xnplus1[3] = 1
        l_t[1] += 2^(transition-1) * p[7]
        l_tr[1] = possible_transitions[transition]
    end
end
isabsorbing_square_wave(p::Vector{Float64}, xn::Vector{Int}) = (p[1] == 0 && p[3] == 0 && p[5] == 0)
g_square_wave = [:A, :HIGH, :LOW]

square_wave_oscillator = ContinuousTimeModel(d, k, dict_var_square, dict_params_square, l_tr_square, 
                                             p_square, x0_square, t0_square, square_wave_f!, isabsorbing_square_wave; 
                                             g=g_square_wave, name="square wave oscillator pkg", time_bound = 105.0)

export square_wave_oscillator

