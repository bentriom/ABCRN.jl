
function _resize_trajectory!(values::Vector{Vector{Int}}, times::Vector{Float64}, 
                             transitions::Vector{Transition}, size::Int)
    for i = eachindex(values) resize!((values[i]), size) end
    resize!(times, size)
    resize!(transitions, size)
end

function _finish_bounded_trajectory!(values::Vector{Vector{Int}}, times::Vector{Float64}, 
                                    transitions::Vector{Transition}, time_bound::Float64)
    for i = eachindex(values) push!((values[i]), (values[i][end])) end
    push!(times, time_bound)
    push!(transitions, nothing)
end

function _update_values!(values::Vector{Vector{Int}}, times::Vector{Float64}, transitions::Vector{Transition},
                         xn::Vector{Int}, tn::Float64, tr_n::Transition, idx::Int)
    for k = eachindex(values) values[k][idx] = xn[k] end
    (times[idx] = tn)
    (transitions[idx] = tr_n)
end

function generate_code_simulation(model_name::Symbol, f!::Symbol, isabsorbing::Symbol)

    return quote
        import MarkovProcesses: simulate
        
        """
        `simulate(m)`

        Simulates a model. If `m::SynchronizedModel`, the simulation is synchronized with a 
        Linear Hybrid Automaton.
        """
        function simulate(m::$(model_name); p::Union{Nothing,AbstractVector{Float64}} = nothing)
            p_sim = getfield(m, :p)
            if p != nothing
                p_sim = p
            end
            time_bound = getfield(m, :time_bound)
            buffer_size = getfield(m, :buffer_size)
            estim_min_states = getfield(m, :estim_min_states)
            # First alloc
            full_values = Vector{Vector{Int}}(undef, getfield(m, :dim_state))
            for i = eachindex(full_values) full_values[i] = zeros(Int, estim_min_states) end
            times = zeros(Float64, estim_min_states)
            transitions = Vector{Transition}(undef, estim_min_states)
            # Initial values
            for i = eachindex(full_values) full_values[i][1] = getfield(m, :x0)[i] end
            times[1] = getfield(m, :t0)
            transitions[1] = nothing
            # Values at time n
            n = 1
            xn = copy(getfield(m, :x0))
            tn = getfield(m, :t0) 
            # at time n+1
            isabsorbing::Bool = $(isabsorbing)(p_sim,xn)
            # If x0 is absorbing
            if isabsorbing
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, 1)
                values = full_values[getfield(m, :_g_idx)]
                if isbounded(m)
                    MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
                end
                return Trajectory(m, values, times, transitions)
            end
            # Alloc of vectors where we store n+1 values
            vec_x = zeros(Int, getfield(m, :dim_state))
            l_t = Float64[0.0]
            l_tr = Transition[nothing]
            # First we fill the allocated array
            for i = 2:estim_min_states
                $(f!)(vec_x, l_t, l_tr, xn, tn, p_sim)
                tn = l_t[1]
                isabsorbing = vec_x == xn
                #isabsorbing = $(isabsorbing)(p_sim,xn)
                if isabsorbing || tn > time_bound
                    break
                end
                # n+1 values are now n values
                n += 1
                copyto!(xn, vec_x)
                MarkovProcesses._update_values!(full_values, times, transitions, xn, tn, l_tr[1], i)
            end
            # If simulation ended before the estimation of states
            if n < estim_min_states
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, n)
                values = full_values[getfield(m, :_g_idx)]
                if isbounded(m)
                    MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
                end
                return Trajectory(m, values, times, transitions)
            end
            # Otherwise, buffering system
            size_tmp = 0
            while !isabsorbing && tn <= time_bound
                # Alloc buffer
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, estim_min_states+size_tmp+buffer_size)
                i = 0
                while i < buffer_size
                    i += 1
                    $(f!)(vec_x, l_t, l_tr, xn, tn, p_sim)
                    tn = l_t[1]
                    isabsorbing = vec_x == xn
                    #isabsorbing = $(isabsorbing)(p_sim,xn)
                    if isabsorbing || tn > time_bound
                        i -= 1
                        break
                    end
                    # n+1 values are now n values
                    copyto!(xn, vec_x)
                    MarkovProcesses._update_values!(full_values, times, transitions, 
                                    xn, tn, l_tr[1], estim_min_states+size_tmp+i)
                end
                # If simulation ended before the end of buffer
                if i < buffer_size
                    MarkovProcesses._resize_trajectory!(full_values, times, transitions, estim_min_states+size_tmp+i)
                end
                size_tmp += i
                n += i
            end
            values = full_values[getfield(m, :_g_idx)]
            if isbounded(m)
                # Add last value: the convention is the last transition is nothing,
                # the trajectory is bounded
                MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
            end
            return Trajectory(m, values, times, transitions)
        end
    end
end

function generate_code_synchronized_simulation(model_name::Symbol, lha_name::Symbol, 
                                               edge_type::Symbol, f!::Symbol, isabsorbing::Symbol)
    
    return quote
        import MarkovProcesses: simulate, volatile_simulate

        function simulate(m::$(model_name), A::$(lha_name), product::SynchronizedModel,
                          p_sim::AbstractVector{Float64}, verbose::Bool)
            x0 = getfield(m, :x0)
            t0 = getfield(m, :t0)
            time_bound = getfield(m, :time_bound)
            buffer_size = getfield(m, :buffer_size)
            estim_min_states = getfield(m, :estim_min_states)
            edge_candidates = Vector{$(edge_type)}(undef, 2)
            # First alloc
            full_values = Vector{Vector{Int}}(undef, getfield(m, :dim_state))
            for i = eachindex(full_values) full_values[i] = zeros(Int, estim_min_states) end
            times = zeros(Float64, estim_min_states)
            transitions = Vector{Transition}(undef, estim_min_states)
            # Initial values
            for i = eachindex(full_values) full_values[i][1] = x0[i] end
            times[1] = t0
            transitions[1] = nothing
            S = init_state(A, x0, t0)
            ptr_loc_state = [S.loc]
            values_state = S.values
            ptr_time_state = [S.time]
            # Values at time n
            n = 1
            xn = copy(x0)
            tn = copy(t0) 
            next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                        x0, t0, nothing, x0, p_sim, edge_candidates; verbose = verbose)
            isabsorbing::Bool = $(isabsorbing)(p_sim,xn)
            isacceptedLHA::Bool = isaccepted(ptr_loc_state[1], A)
            # Alloc of vectors where we stock n+1 values
            vec_x = zeros(Int, getfield(m, :dim_state))
            l_t = Float64[0.0]
            l_tr = Transition[nothing]
            # If x0 is absorbing
            if isabsorbing || isacceptedLHA 
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, 1)
                values = full_values[getfield(m, :_g_idx)]
                if isbounded(m)
                    MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
                end
                if isabsorbing && !isacceptedLHA
                    next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                                xn, time_bound, nothing, xn, p_sim, edge_candidates; verbose = verbose)
                end
                setfield!(S, :loc, ptr_loc_state[1])
                setfield!(S, :time, ptr_time_state[1])
                return SynchronizedTrajectory(S, product, values, times, transitions)
            end
            # First we fill the allocated array
            for i = 2:estim_min_states
                $(f!)(vec_x, l_t, l_tr, xn, tn, p_sim)
                tn = l_t[1]
                isabsorbing = vec_x == xn
                if isabsorbing || tn > time_bound
                    break
                end
                next_state!(A, ptr_loc_state, values_state, ptr_time_state, vec_x, l_t[1], l_tr[1], xn, p_sim, edge_candidates; verbose = verbose)
                # n+1 values are now n values
                n += 1
                copyto!(xn, vec_x)
                tr_n = l_tr[1]
                MarkovProcesses._update_values!(full_values, times, transitions, xn, tn, tr_n, i)
                isacceptedLHA = isaccepted(ptr_loc_state[1], A)
                if isacceptedLHA 
                    break
                end
            end
            # If simulation ended before the estimation of states
            if n < estim_min_states
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, n)
                values = full_values[getfield(m, :_g_idx)]
                if isbounded(m)
                    MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
                end
                if isabsorbing && !isacceptedLHA
                    next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                                xn, time_bound, nothing, xn, p_sim, edge_candidates; verbose = verbose)
                end
                setfield!(S, :loc, ptr_loc_state[1])
                setfield!(S, :time, ptr_time_state[1])
                return SynchronizedTrajectory(S, product, values, times, transitions)
            end
            # Otherwise, buffering system
            size_tmp = 0
            while !isabsorbing && tn <= time_bound && !isacceptedLHA
                # Alloc buffer
                MarkovProcesses._resize_trajectory!(full_values, times, transitions, estim_min_states+size_tmp+buffer_size)
                i = 0
                while i < buffer_size
                    i += 1
                    $(f!)(vec_x, l_t, l_tr, xn, tn, p_sim)
                    tn = l_t[1]
                    isabsorbing = vec_x == xn
                    if isabsorbing || tn > time_bound
                        i -= 1
                        break
                    end
                    next_state!(A, ptr_loc_state, values_state, ptr_time_state, vec_x, l_t[1], l_tr[1], xn, p_sim, edge_candidates; verbose = verbose)
                    # n+1 values are now n values
                    copyto!(xn, vec_x)
                    tr_n = l_tr[1]
                    MarkovProcesses._update_values!(full_values, times, transitions, 
                                                    xn, tn, tr_n, estim_min_states+size_tmp+i)
                    isacceptedLHA = isaccepted(ptr_loc_state[1], A)
                    if isacceptedLHA
                        break
                    end
                end
                # If simulation ended before the end of buffer
                if i < buffer_size
                    MarkovProcesses._resize_trajectory!(full_values, times, transitions, estim_min_states+size_tmp+i)
                end
                size_tmp += i
                n += i
            end
            values = full_values[getfield(m, :_g_idx)]
            if isbounded(m) && !isaccepted(ptr_loc_state[1], A)
                # Add last value: the convention is that if the last transition is nothing,
                # the trajectory is bounded
                MarkovProcesses._finish_bounded_trajectory!(values, times, transitions, time_bound)
            end
            if isabsorbing && !isacceptedLHA
                next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                            xn, time_bound, nothing, xn, p_sim, edge_candidates; verbose = verbose)
            end
            setfield!(S, :loc, ptr_loc_state[1])
            setfield!(S, :time, ptr_time_state[1])
            return SynchronizedTrajectory(S, product, values, times, transitions)
        end
        
        function volatile_simulate(m::$(model_name), A::$(lha_name), p_sim::AbstractVector{Float64}, 
                                   epsilon::Float64, verbose::Bool)
            if $(Meta.quot(lha_name)) == :ABCEuclideanDistanceAutomaton A.ϵ = epsilon end
            x0 = getfield(m, :x0)
            t0 = getfield(m, :t0)
            time_bound = getfield(m, :time_bound)
            S = init_state(A, x0, t0)
            ptr_loc_state = [S.loc]
            values_state = S.values
            ptr_time_state = [S.time]
            edge_candidates = Vector{$(edge_type)}(undef, 2)
            # Values at time n
            n = 1
            xn = copy(x0)
            tn = copy(t0) 
            next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                        x0, t0, nothing, x0, p_sim, edge_candidates; verbose = verbose)
            isabsorbing::Bool = $(isabsorbing)(p_sim,xn)
            isacceptedLHA::Bool = isaccepted(ptr_loc_state[1], A)
            # Alloc of vectors where we stock n+1 values
            vec_x = zeros(Int, getfield(m, :dim_state))
            l_t = Float64[0.0]
            l_tr = Transition[nothing]
            # If x0 is absorbing
            if isabsorbing || isacceptedLHA 
                if !isacceptedLHA
                    next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                                xn, time_bound, nothing, xn, p_sim, edge_candidates; verbose = verbose)
                end
                setfield!(S, :loc, ptr_loc_state[1])
                setfield!(S, :time, ptr_time_state[1])
                return S
            end
            while !isabsorbing && tn <= time_bound && !isacceptedLHA
                $(f!)(vec_x, l_t, l_tr, xn, tn, p_sim)
                if l_t[1] > time_bound
                    tn = l_t[1]
                    break
                end
                isabsorbing = vec_x == xn
                if isabsorbing
                    break
                end
                next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                            vec_x, l_t[1], l_tr[1], xn, p_sim, edge_candidates; verbose = verbose)
                # n+1 values are now n values
                copyto!(xn, vec_x)
                tn = l_t[1]
                tr_n = l_tr[1]
                isacceptedLHA = isaccepted(ptr_loc_state[1], A)
                n += 1
            end
            if isabsorbing && !isacceptedLHA
                next_state!(A, ptr_loc_state, values_state, ptr_time_state, 
                            xn, time_bound, nothing, xn, p_sim, edge_candidates; verbose = verbose)
            end
            setfield!(S, :loc, ptr_loc_state[1])
            setfield!(S, :time, ptr_time_state[1])
            return S
        end
    end
end

"""
`volatile_simulate(sm::SynchronizedModel; p, verbose)`

Simulates a model synchronized with an automaton but does not store the values of the simulation
in order to improve performance.
It returns the last state of the simulation `S::StateLHA` not a trajectory `σ::SynchronizedTrajectory`.
"""
function volatile_simulate(product::SynchronizedModel; 
                           p::Union{Nothing,AbstractVector{Float64}} = nothing, epsilon::Float64 = 0.0, verbose::Bool = false)
    m = product.m
    A = product.automaton
    p_sim = getfield(m, :p)
    if p != nothing
        p_sim = p
    end
    S = volatile_simulate(m, A, p_sim, epsilon, verbose)
    return S
end

function simulate(product::SynchronizedModel; 
    p::Union{Nothing,AbstractVector{Float64}} = nothing, verbose::Bool = false)
    m = getfield(product, :m)
    A = getfield(product, :automaton)
    p_sim = getfield(m, :p)
    if p != nothing
        p_sim = p
    end
    σ = simulate(m, A, product, p_sim, verbose)
    return σ
end

"""
    `simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})

Simulates the model contained in pm with p_prior values.
It simulates the model by taking the parameters contained in get_proba_model(pm).p and
replace the 1D parameters pm.params with p_prior.
"""
function simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})
    full_p = copy(get_proba_model(pm).p)
    full_p[pm._param_idx] = p_prior
    
    return simulate(pm.m; p = full_p) 
end
"""
    `volatile_simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})

A volatile version of `simulate(pm::ParametricModel, p_prior::AbstractVector{Float64})`.
The model in pm should be of type SynchronizedModel (`typeof(pm.m) <: SynchronizedModel`).
It returns `S::StateLHA`, not a trajectory.
"""
function volatile_simulate(pm::ParametricModel, p_prior::AbstractVector{Float64};
                           epsilon::Union{Nothing,Float64} = nothing)
    @assert typeof(pm.m) <: SynchronizedModel
    full_p = copy(get_proba_model(pm).p)
    full_p[pm._param_idx] = p_prior
    
    return volatile_simulate(pm.m; p = full_p, epsilon = epsilon)
end
"""
    `distribute_mean_value_lha(sm::SynchronizedModel, sym_var::Symbol, nbr_stim::Int)`

Distribute over workers the computation of the mean value 
of a LHA over `nbr_sim` simulations of the model.
"""
function distribute_mean_value_lha(sm::SynchronizedModel, sym_var::VariableAutomaton, nbr_sim::Int; with_accepts::Bool = false)
    sum_val = @distributed (+) for i = 1:nbr_sim
        S = volatile_simulate(sm)
        if !with_accepts
            S[sym_var]
        else
            [S[sym_var], isaccepted(S)]
        end
    end
    return sum_val / nbr_sim
end

function distribute_mean_value_lha(sm::SynchronizedModel, sym_var::Vector{VariableAutomaton}, nbr_sim::Int; with_accepts::Bool = false)
    sum_val = @distributed (+) for i = 1:nbr_sim 
        S = volatile_simulate(sm)
        if !with_accepts
            S[sym_var]
        else
            vcat(S[sym_var], isaccepted(S))
        end
    end
    return sum_val / nbr_sim
end

function mean_value_lha(sm::SynchronizedModel, sym_var::VariableAutomaton, nbr_sim::Int)
    sum_val = 0.0 
    for i = 1:nbr_sim 
        sum_val += volatile_simulate(sm)[sym_var] 
    end
    return sum_val / nbr_sim
end
"""
    `distribute_var_value_lha(sm::SynchronizedModel, nbr_sim::Int, value = 0, sym_var = :d)

Compute the probability that the variable `sym_var` is equal to  `value`
of a LHA over `nbr_sim` simulations of the model.
"""
function probability_var_value_lha(sm::SynchronizedModel, nbr_sim::Int; 
                                   value::Float64 = 0.0, sym_var::VariableAutomaton = :d)
    sum_val = @distributed (+) for i = 1:nbr_sim
        S = volatile_simulate(sm)
        S[sym_var] == value
    end
    return sum_val / nbr_sim
end

number_simulations_smc_chernoff(approx::Float64, conf::Float64) = log(2/(1-conf)) / (2*approx^2)
function smc_chernoff(sm::SynchronizedModel; approx::Float64 = 0.01, confidence::Float64 = 0.99)
    @assert typeof(sm.automaton) <: Union{AutomatonF,AutomatonG,AutomatonGandF}
    nbr_sim = number_simulations_smc_chernoff(approx, confidence) 
    nbr_sim = convert(Int, trunc(nbr_sim)+1)
    return probability_var_value_lha(sm, nbr_sim)
end

function distribute_prob_accept_lha(sm::SynchronizedModel, nbr_sim::Int)
    sum_val = @distributed (+) for i = 1:nbr_sim 
        Int(isaccepted(volatile_simulate(sm))) 
    end
    return sum_val / nbr_sim
end

function Base.show(io::IO, m::ContinuousTimeModel)
    print(io, "$(typeof(m)) <: ContinuousTimeModel model\n")
    print(io, "- variables :\n")
    for (var, idx) in m.map_var_idx
        print(io, "* $var (index = $idx in state space)\n")
    end
    print(io, "- parameters :\n")
    for (param, idx) in m.map_param_idx
        print(io, "* $param (index = $idx in parameter space)\n")
    end
    print(io, "- transitions : $(join(m.transitions,','))\n")
    print(io, "- observed variables :\n")
    for i in eachindex(m.g)
        print(io, "* $(m.g[i]) (index = $i in observed state space, index = $(m._g_idx[i]) in state space)\n")
    end
    print(io, "p = $(m.p)\n")
    print(io, "x0 = $(m.x0)\n")
    print(io, "t0 = $(m.t0)\n")
    print(io, "time bound = $(m.time_bound)")
end

isbounded(m::Model) = get_proba_model(m).time_bound < Inf
function check_consistency(m::ContinuousTimeModel) 
    @assert m.dim_state == length(m.map_var_idx) 
    @assert m.dim_state == length(m.x0) 
    @assert m.dim_params == length(m.map_param_idx)
    @assert m.dim_params == length(m.p)
    @assert length(m.g) <= m.dim_state
    @assert length(m._g_idx) == length(m.g)
    @assert m.buffer_size >= 0
    #@assert typeof(getfield(Main, m.isabsorbing)(m.p, m.x0)) == Bool
    return true
end

# Set and get Model fields
function set_observed_var!(am::Model, g::Vector{VariableModel})
    m = get_proba_model(am)
    dim_obs_state = length(g)
    _map_obs_var_idx = Dict{VariableModel}{Int}()
    _g_idx = zeros(Int, dim_obs_state)
    for i = 1:dim_obs_state
        _g_idx[i] = m.map_var_idx[g[i]] # = ( (g[i] = i-th obs var)::VariableModel => idx in state space )
        _map_obs_var_idx[g[i]] = i
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = _map_obs_var_idx
end
function observe_all!(am::Model)
    m = get_proba_model(am)
    g = Vector{VariableModel}(undef, m.dim_state)
    _g_idx = collect(1:m.dim_state)
    for var in keys(m.map_var_idx)
        g[m.map_var_idx[var]] = var
    end
    m.g = g
    m._g_idx = _g_idx
    m._map_obs_var_idx = m.map_var_idx
end
function set_param!(am::Model, new_p::Vector{Float64})
    m = get_proba_model(am)
    @assert length(new_p) == m.dim_params "New parameter vector hasn't the same dimension of parameter space"
    m.p = new_p
end
function set_param!(am::Model, name_p::ParameterModel, p_i::Float64) 
    m = get_proba_model(am)
    m.p[m.map_param_idx[name_p]] = p_i
end
function set_param!(am::Model, l_name_p::Vector{ParameterModel}, p::Vector{Float64}) 
    m = get_proba_model(am)
    @assert length(l_name_p) == length(p) "Parameter names vector and parameter values haven't the same dimensions"
    for i = eachindex(l_name_p)
        set_param!(m, l_name_p[i], p[i])
    end
end
function set_x0!(am::Model, l_name_var::Vector{VariableModel}, x0::Vector{Int})
    m = get_proba_model(am)
    @assert length(l_name_var) == length(x0) "State names vector state values haven't the same dimensions"
    for i = eachindex(l_name_var)
        set_x0!(m, l_name_var[i], x0[i])
    end
end
function set_x0!(am::Model, name_var::VariableModel, var_i::Int) 
    m = get_proba_model(am)
    m.x0[m.map_var_idx[name_var]] = var_i
end
function set_x0!(am::Model, new_x0::Vector{Int})
    m = get_proba_model(am)
    @assert length(new_x0) == m.dim_state "New x0 vector hasn't the dimension of state space"
    m.x0 = new_x0
end
set_time_bound!(am::Model, b::Float64) = (get_proba_model(am).time_bound = b)

get_param(am::Model) = get_proba_model(am).p
get_x0(am::Model) = get_proba_model(am).x0
function getindex(am::Model, name_p::ParameterModel)
    m = get_proba_model(am)
    m.p[m.map_param_idx[name_p]]
end
function getproperty(m::ContinuousTimeModel, sym::Symbol)
    if sym == :dim_obs_state
        return length(m.g)
    else
        return getfield(m, sym)
    end
end
function getproperty(pm::ParametricModel, sym::Symbol)
    if sym == :df
        return length(pm.params)
    else
        return getfield(pm, sym)
    end
end

get_proba_model(m::ContinuousTimeModel) = m
get_proba_model(sm::SynchronizedModel) = sm.m
get_proba_model(pm::ParametricModel) = get_proba_model(pm.m)

get_observed_var(m::Model) = get_proba_model(am).g

# Prior methods
"""
    `draw_model!(pm::ParametricModel)`

Draw a parameter from the prior disitribution defined in `pm::ParametricModel`
and save it in the model contained in `pm`.
"""
draw_model!(pm::ParametricModel) = set_param!(get_proba_model(pm), pm.params, rand(pm.distribution))
"""
    `draw!(vec_p, pm)`

Draw a parameter from the prior distribution defined in pm and stores it in vec_p.
"""
draw!(vec_p::AbstractVector{Float64}, pm::ParametricModel) = rand!(pm.distribution, vec_p)
"""
    `draw!(mat_p, pm)`

Draw `size(mat_p)[2]` (number of columns of mat_p) parameters from the prior distribution 
defined in pm and stores it in mat_p.
"""
function draw!(mat_p::AbstractMatrix{Float64}, pm::ParametricModel)
    for i = 1:size(mat_p)[2]
        draw!(view(mat_p, :, i), pm)
    end
end
"""
    `prior_pdf(p_prior, pm)`
 
Computes the density of the prior distribution defined in pm on argument p_prior.
"""
prior_pdf(pm::ParametricModel, p_prior::AbstractVector{Float64}) = pdf(pm.distribution, p_prior)
"""
    `prior_pdf(vec_res, mat_p, pm)`
 
Computes the density of the prior distribution defined in pm on each column
ov mat_p. Stores it in mat_p. (`length(vec_res) == size(mat_p)[2]`)
"""
function prior_pdf!(res_pdf::AbstractVector{Float64}, pm::ParametricModel, mat_p::AbstractMatrix{Float64})
    for i = eachindex(res_pdf)
        res_pdf[i] = prior_pdf(pm, view(mat_p, :, i))
    end
end
fill!(pm::ParametricModel, p_prior::AbstractVector{Float64}) = get_proba_model(pm).p[pm._param_idx] = p_prior
insupport(pm::ParametricModel, p_prior::AbstractVector{Float64}) = insupport(pm.distribution, p_prior)

# to do: simulate(pm), create_res_dir, check_bounds, ajouter un champ _param_idx pour pm.
#

