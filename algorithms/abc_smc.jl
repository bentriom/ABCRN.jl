
import StatsBase: mean, median, std, cov, ProbabilityWeights
import Statistics: quantile
import NearestNeighbors: KDTree, knn 
import Distributions: MvNormal, Categorical
import Random: rand!

import Distributed: @sync, @async, nworkers, nprocs, workers
import DistributedArrays: DArray, dzeros, convert, localpart
using Distributed
using LinearAlgebra
using DelimitedFiles
using Logging

include("_utils_abc.jl")

struct ResultAbc
	mat_p_end::Matrix{Float64}
	mat_cov::Matrix{Float64}
    nbr_sim::Int64
    exec_time::Float64
    vec_dist::Vector{Float64}
    epsilon::Float64
    weights::Vector{Float64}
    l_ess::Vector{Float64}
end

function automaton_abc(pm::ParametricModel; 
                       nbr_particles::Int = 100, tolerance::Float64 = 0.0, alpha::Float64 = 0.75, kernel_type = "mvnormal", 
                       NT::Float64 = nbr_particles/2, duration_time::Float64 = Inf, dir_results::Union{Nothing,String} = nothing,
                       bound_sim::Int = typemax(Int), sym_var_aut::VariableAutomaton = :d, save_iterations::Bool = false, 
                       init_mat_p::Union{Matrix{Float64},Nothing} = nothing, 
                       init_weights::Union{Vector{Float64},Nothing} = nothing, init_vec_dist::Union{Vector{Float64},Nothing} = nothing,
                       verbose::Int = 0) 
    @assert typeof(pm.m) <: SynchronizedModel "Automaton-ABC is defined for synchronized models only"
    @assert 0 < nbr_particles
    @assert 0.0 < alpha < 1.0
    @assert kernel_type in ["mvnormal", "knn_mvnormal"]
    test_init = (init_mat_p == nothing && init_weights == nothing) || 
                (size(init_mat_p, 2) == length(init_weights) == length(init_vec_dist))
    @assert test_init "All initialisation data should be specified and have the correct dimension"
    if dir_results != nothing
        dir_results = basename(dir_results) != "" ? dir_results * "/" : dir_results 
        if !isdir(dir_results) mkdir(dir_results) end
        file_cfg = open(dir_results * "config_abc.out", "w")
        write(file_cfg, "ParametricModel : $(pm) \n")
        write(file_cfg, "Number of particles : $(nbr_particles) \n")
        write(file_cfg, "Final tolerance : $(tolerance) \n")
        write(file_cfg, "alpha : $(alpha) \n")
        write(file_cfg, "kernel type : $(kernel_type) \n")
        close(file_cfg)
    end
    if nprocs() == 1
        return _abc_smc(pm, nbr_particles, tolerance, alpha, kernel_type, NT, duration_time, bound_sim, dir_results, sym_var_aut, save_iterations, init_mat_p, init_weights, init_vec_dist)
    end
    return _distributed_abc_smc(pm, nbr_particles, tolerance, alpha, kernel_type, NT, duration_time, bound_sim, dir_results, sym_var_aut, save_iterations, init_mat_p, init_weights, init_vec_dist)
end

"""
    `abc_smc(pm::ParametricModel, l_obs, func_dist; nbr_particles, alpha, kernel_type, NT
                                               duration_tiùe, bound_sim, sym_var_aut, verbose)`

Run the ABC-SMC algorithm with the pm parametric model. 
`func_dist(l_sim, l_obs)` is the distance function between simulations and observation, 
it corresponds to \$\rho(\eta(y_sim), \eta(y_exp))\$.
`l_obs::Vector{<:T2}` is a collection of observations.
`dist` must have a signature of the form `func_dist(l_sim::Vector{T1}, l_obs::Vector{T2})`.
If pm is defined on a ContinuousTimeModel, then `T1` should verify `T1 <: Trajectory`.

!!! Distance function and distributed ABC

    If you use `abc_smc` with multiple workers, `dist` has to be defined 
    on every workers by using @everywhere.

"""
function abc_smc(pm::ParametricModel, l_obs::AbstractVector, func_dist::Function; 
                 nbr_particles::Int = 100, tolerance::Float64 = 0.0, alpha::Float64 = 0.75, kernel_type = "mvnormal", 
                 NT::Float64 = nbr_particles/2, duration_time::Float64 = Inf, dir_results::Union{Nothing,String} = nothing,
                 bound_sim::Int = typemax(Int), sym_var_aut::VariableAutomaton = :d, save_iterations::Bool = false, 
                 init_mat_p::Union{Matrix{Float64},Nothing} = nothing, 
                 init_weights::Union{Vector{Float64},Nothing} = nothing, init_vec_dist::Union{Vector{Float64},Nothing} = nothing,
                 verbose::Int = 0) 
    @assert 0 < nbr_particles
    @assert 0.0 < alpha < 1.0
    @assert kernel_type in ["mvnormal", "knn_mvnormal"]
    test_init = (init_mat_p == nothing && init_weights == nothing) || 
                (size(init_mat_p, 2) == length(init_weights) == length(init_vec_dist))
    @assert test_init "All initialisation data should be specified and have the correct dimension"
    if dir_results != nothing
        dir_results = basename(dir_results) != "" ? dir_results * "/" : dir_results 
        if !isdir(dir_results) mkdir(dir_results) end
        file_cfg = open(dir_results * "config_abc.out", "w")
        write(file_cfg, "Configuration of ABC algorithm\n")
        write(file_cfg, "ParametricModel : $(pm) \n")
        write(file_cfg, "Number of particles : $(nbr_particles) \n")
        write(file_cfg, "Final tolerance : $(tolerance) \n")
        write(file_cfg, "alpha : $(alpha) \n")
        write(file_cfg, "kernel type : $(kernel_type) \n")
        close(file_cfg)
    end
    if nprocs() == 1
        return _abc_smc(pm, nbr_particles, tolerance, alpha, kernel_type, NT, duration_time, bound_sim, dir_results, sym_var_aut, save_iterations, init_mat_p, init_weights, init_vec_dist; l_obs = l_obs, func_dist = func_dist)
    end
    return _distributed_abc_smc(pm, nbr_particles, tolerance, alpha, kernel_type, NT, duration_time, bound_sim, dir_results, sym_var_aut, save_iterations, init_mat_p, init_weights, init_vec_dist; l_obs = l_obs, func_dist = func_dist)
end


# To code: 
# Pkg related: draw!, prior_density!

function _abc_smc(pm::ParametricModel, nbr_particles::Int, tolerance::Float64, alpha::Float64, 
                  kernel_type::String, NT::Float64, duration_time::Float64, 
                  bound_sim::Int, dir_results::Union{Nothing,String}, sym_var_aut::VariableAutomaton, save_iterations::Bool,
                  init_mat_p::Union{Matrix{Float64},Nothing}, 
                  init_weights::Union{Vector{Float64},Nothing}, init_vec_dist::Union{Vector{Float64},Nothing}; 
                  l_obs::Union{Nothing,AbstractVector} = nothing, func_dist::Union{Nothing,Function} = nothing)
    @info "ABC PMC with $(nworkers()) processus and $(Threads.nthreads()) threads"
    begin_time = time()
    nbr_p = pm.df
    last_epsilon = tolerance
    # Init. Iteration 1
    t = 1
    epsilon = Inf
    mat_p_old = zeros(nbr_p, nbr_particles)
    vec_dist = zeros(nbr_particles)
    wl_old = zeros(nbr_particles)
    @info "Step 1 : Init"
    if init_mat_p == nothing
        _init_abc_automaton!(mat_p_old, vec_dist, pm, sym_var_aut; l_obs = l_obs, func_dist = func_dist)
        prior_pdf!(wl_old, pm, mat_p_old)
        normalize!(wl_old, 1)
    else
        mat_p_old = init_mat_p
        wl_old = init_weights
        vec_dist = init_vec_dist
    end
    ess = effective_sample_size(wl_old)
    l_ess = zeros(0)
    l_ess = push!(l_ess, ess)
    flush(stdout)
    nbr_tot_sim = nbr_particles 
    current_time = time()
    old_epsilon = epsilon 
	mat_p = zeros(nbr_p, nbr_particles)
    wl_current = zeros(nbr_particles)
    l_nbr_sim = zeros(Int, nbr_particles) 
    while (epsilon > last_epsilon) && (current_time - begin_time <= duration_time) && (nbr_tot_sim <= bound_sim)
        if dir_results != "" && save_iterations
            step_dir = dir_results * "/step_$t/"
            if !isdir(step_dir) mkdir(step_dir) end
            writedlm(step_dir * "weights.csv", wl_old, ',')
            writedlm(step_dir * "mat_p.csv", mat_p_old, ',')
            writedlm(step_dir * "vec_dist.csv", vec_dist, ',')
        end
        t += 1
        begin_time_ite = time()
        @info "Step $t"
        # Set new epsilon
        epsilon = _compute_epsilon(vec_dist, alpha, old_epsilon, last_epsilon)
        @info "Current epsilon" epsilon
        @debug "5 first dist values" sort(vec_dist)[1:5]
        @debug mean(vec_dist), maximum(vec_dist), median(vec_dist), std(vec_dist)
        
        # Broadcast simulations
        mat_cov = nothing
        tree_mat_p = nothing
        if kernel_type == "mvnormal"
            mat_cov = 2 * cov(mat_p_old, ProbabilityWeights(wl_old), 2; corrected=false)
            @debug diag(mat_cov)
            if det(mat_cov) == 0.0
                @debug det(mat_cov), rank(mat_cov), effective_sample_size(wl_old)
                @error "Bad inv mat cov"
            end
        end
        if kernel_type == "knn_mvnormal"
            tree_mat_p = KDTree(mat_p_old)
        end
        Threads.@threads for i = eachindex(vec_dist)
            _update_param!(mat_p, vec_dist, l_nbr_sim, wl_current, i, pm, epsilon,
                           wl_old, mat_p_old, mat_cov, tree_mat_p, kernel_type, sym_var_aut;
                           l_obs = l_obs, func_dist = func_dist)
        end
        normalize!(wl_current, 1)
        step_nbr_sim = sum(l_nbr_sim)
        nbr_tot_sim += step_nbr_sim 
        ess = effective_sample_size(wl_current)
        l_ess = push!(l_ess, ess)
        @debug ess
        # Resampling
        if ess < NT
            @info "Resampling.."
            d_weights = Categorical(wl_current)
            ind_w = rand(d_weights, nbr_particles)
            mat_p = mat_p[:,ind_w]
            wl_current = ones(nbr_particles)
            normalize!(wl_current, 1)
            @info "End"
        end
        current_time = time()
        @info "After this step, time spent and number of simulations" steptime=(current_time-begin_time_ite) step_nbr_sim nbr_tot_sim
        mat_p_old = copy(mat_p)
        wl_old = copy(wl_current)
        fill!(l_nbr_sim, 0)
        flush(stdout)
        old_epsilon = epsilon
	end

    mat_cov = cov(mat_p_old, ProbabilityWeights(wl_old), 2; corrected=false)
    if dir_results != nothing
        writedlm(dir_results * "weights_end.csv", wl_old, ',')
        writedlm(dir_results * "mat_p_end.csv", mat_p_old, ',')
        writedlm(dir_results * "vec_dist_end.csv", vec_dist, ',')
        file_cfg = open(dir_results * "results_abc.out", "w")
        write(file_cfg, "\n")
        write(file_cfg, "About the results: \n")
        write(file_cfg, "Total number of simulations: $nbr_tot_sim\n")
        write(file_cfg, "Final tolerance : $(old_epsilon) \n")
        write(file_cfg, "Execution time: $(time() - begin_time)\n")
        write(file_cfg, "Number of jobs: $(nprocs())\n")
        write(file_cfg, "Number of threads: $(Threads.nthreads())\n")
        close(file_cfg)
    end
    r = ResultAbc(mat_p_old, mat_cov, nbr_tot_sim, time() - begin_time, vec_dist, old_epsilon, wl_old, l_ess)
    return r
end

function _distributed_abc_smc(pm::ParametricModel, nbr_particles::Int, tolerance::Float64, alpha::Float64, 
                              kernel_type::String, NT::Float64, duration_time::Float64, bound_sim::Int, 
                              dir_results::Union{Nothing,String}, sym_var_aut::VariableAutomaton, save_iterations::Bool,
                              init_mat_p::Union{Matrix{Float64},Nothing}, 
                              init_weights::Union{Vector{Float64},Nothing}, init_vec_dist::Union{Vector{Float64},Nothing};
                              l_obs::Union{Nothing,AbstractVector} = nothing, func_dist::Union{Nothing,Function} = nothing)
    @info "Distributed ABC PMC with $(nworkers()) processus and $(Threads.nthreads()) threads"
    begin_time = time()
    nbr_p = pm.df
    last_epsilon = tolerance
    # Init. Iteration 1
    t = 1
    epsilon = Inf
    mat_p_old = zeros(nbr_p, nbr_particles)
    vec_dist = zeros(nbr_particles)
    wl_old = zeros(nbr_particles)
    @info "Step 1 : Init"
    if init_mat_p == nothing
        _init_abc_automaton!(mat_p_old, vec_dist, pm, sym_var_aut; l_obs = l_obs, func_dist = func_dist)
        prior_pdf!(wl_old, pm, mat_p_old)
        normalize!(wl_old, 1)
    else
        mat_p_old = init_mat_p
        wl_old = init_weights
        vec_dist = init_vec_dist
    end
    ess = effective_sample_size(wl_old)
    l_ess = zeros(0)
    l_ess = push!(l_ess, ess)
    flush(stdout)
    nbr_tot_sim = nbr_particles 
    current_time = time()
    old_epsilon = epsilon 
	d_mat_p = dzeros(nbr_p, nbr_particles)
    d_vec_dist = dzeros(nbr_particles)
    d_wl_current = dzeros(nbr_particles)
    mat_p = zeros(0,0)
    wl_current = zeros(0)
    while (epsilon > last_epsilon) && (current_time - begin_time <= duration_time) && (nbr_tot_sim <= bound_sim)
        if dir_results != "" && save_iterations
            step_dir = dir_results * "/step_$t/"
            if !isdir(step_dir) mkdir(step_dir) end
            writedlm(step_dir * "weights.csv", wl_old, ',')
            writedlm(step_dir * "mat_p.csv", mat_p_old, ',')
            writedlm(step_dir * "vec_dist.csv", vec_dist, ',')
        end
        t += 1
        begin_time_ite = time()
        @info "Step $t"
        # Set new epsilon
        epsilon = _compute_epsilon(vec_dist, alpha, old_epsilon, last_epsilon)
        @info "Current epsilon" epsilon
        @debug "5 first dist values" sort(vec_dist)[1:5]
        @debug mean(vec_dist), maximum(vec_dist), median(vec_dist), std(vec_dist)
        
        # Broadcast simulations
        mat_cov = nothing
        tree_mat_p = nothing
        if kernel_type == "mvnormal"
            mat_cov = 2 * cov(mat_p_old, ProbabilityWeights(wl_old), 2; corrected=false)
            @debug diag(mat_cov)
            if det(mat_cov) == 0.0
                @debug det(mat_cov), rank(mat_cov), effective_sample_size(wl_old)
                @error "Bad inv mat cov"
            end
        end
        if kernel_type == "knn_mvnormal"
            tree_mat_p = KDTree(mat_p_old)
        end
        l_nbr_sim = zeros(Int, nworkers())
        @sync for id_w in workers()
            t_id_w = id_w - workers()[1] + 1
            @async l_nbr_sim[t_id_w] = 
                remotecall_fetch(() -> _task_worker!(d_mat_p, d_vec_dist, d_wl_current, 
                                                     pm, epsilon, wl_old, mat_p_old, mat_cov, tree_mat_p, 
                                                     kernel_type, sym_var_aut;
                                                     l_obs = l_obs, func_dist = func_dist), id_w)
        end
        wl_current = convert(Array, d_wl_current)
        normalize!(wl_current, 1)
        mat_p = convert(Array, d_mat_p)
        step_nbr_sim = sum(l_nbr_sim)
        nbr_tot_sim += step_nbr_sim 
        ess = effective_sample_size(wl_current)
        l_ess = push!(l_ess, ess)
        @debug ess
        # Resampling
        if ess < NT
            @info "Resampling.."
            d_weights = Categorical(wl_current)
            ind_w = rand(d_weights, nbr_particles)
            mat_p = mat_p[:,ind_w]
            wl_current = ones(nbr_particles)
            normalize!(wl_current, 1)
            @info "End"
        end
        current_time = time()
        @info "After this step, time spent and number of simulations" steptime=(current_time-begin_time_ite) step_nbr_sim nbr_tot_sim
        mat_p_old = mat_p
        wl_old = wl_current
        vec_dist = convert(Array, d_vec_dist)
        fill!(l_nbr_sim, 0)
        flush(stdout)
        old_epsilon = epsilon
	end
    
    mat_cov = cov(mat_p_old, ProbabilityWeights(wl_old), 2; corrected=false)
    save_mat_p_end = true
    if dir_results != nothing
        writedlm(dir_results * "weights_end.csv", wl_old, ',')
        writedlm(dir_results * "mat_p_end.csv", mat_p_old, ',')
        writedlm(dir_results * "vec_dist_end.csv", vec_dist, ',')
        file_cfg = open(dir_results * "results_abc.out", "w")
        write(file_cfg, "\n")
        write(file_cfg, "About the results: \n")
        write(file_cfg, "Total number of simulations: $nbr_tot_sim\n")
        write(file_cfg, "Final tolerance : $(old_epsilon) \n")
        write(file_cfg, "Execution time: $(time() - begin_time)\n")
        write(file_cfg, "Number of jobs: $(nprocs())\n")
        write(file_cfg, "Number of threads: $(Threads.nthreads())\n")
        close(file_cfg)
    end
    r = ResultAbc(mat_p_old, mat_cov, nbr_tot_sim, time() - begin_time, convert(Array, d_vec_dist), old_epsilon, wl_old, l_ess)
    return r
end

