
function _init_abc_automaton!(mat_p_old::Matrix{Float64}, vec_dist::Vector{Float64}, 
                              pm::ParametricModel, sym_var_aut::VariableAutomaton;
                              l_obs::Union{Nothing,AbstractVector} = nothing, 
                              func_dist::Union{Nothing,Function} = nothing)
    vec_p = zeros(pm.df)
    for i = eachindex(vec_dist)
        draw!(vec_p, pm)
        mat_p_old[:,i] = vec_p
        if l_obs == nothing
            S = volatile_simulate(pm, vec_p)
            vec_dist[i] = S[sym_var_aut]
        else
            l_sim = [simulate(pm, vec_p) for i = 1:length(l_obs)]
            vec_dist[i] = func_dist(l_sim, l_obs)
        end
    end
end

function _compute_epsilon(vec_dist::Vector{Float64}, α::Float64, 
                          old_eps::Float64, last_eps::Float64)
    u_vec_dist = unique(vec_dist)
    if length(u_vec_dist) == 2
        epsilon = quantile(u_vec_dist, α)
    else
        epsilon = quantile(vec_dist, α)
    end
    if old_eps == epsilon
        print("eps == old_eps,  we readjust")
        s_dist = sort(u_vec_dist)
        max_dist = 0.0
        for d in s_dist
            if (d < epsilon) && (max_dist < d)
                max_dist = d
            end
        end
        max_dist = (max_dist == 0) ? 0.9*epsilon : max_dist
        epsilon = max_dist
    end
    epsilon = (epsilon < last_eps) ? last_eps : epsilon
end

function _task_worker!(d_mat_p::DArray{Float64,2}, d_vec_dist::DArray{Float64,1},
                       d_wl_current::DArray{Float64,1},
                       pm::ParametricModel, epsilon::Float64,
                       wl_old::Vector{Float64}, mat_p_old::Matrix{Float64},
                       mat_cov::Matrix{Float64}, tree_mat_p::Union{Nothing,KDTree}, 
                       kernel_type::String, sym_var_aut::VariableAutomaton;
                       l_obs::Union{Nothing,AbstractVector} = nothing, 
                       func_dist::Union{Nothing,Function} = nothing)
    mat_p = localpart(d_mat_p)
    vec_dist = localpart(d_vec_dist)
    wl_current = localpart(d_wl_current)
    l_nbr_sim = zeros(Int, length(vec_dist))
    Threads.@threads for i = eachindex(vec_dist)
        _update_param!(mat_p, vec_dist, l_nbr_sim, wl_current, i, pm, epsilon,
                       wl_old, mat_p_old, mat_cov, tree_mat_p, kernel_type, sym_var_aut;
                       l_obs = l_obs, func_dist = func_dist)
    end
    return sum(l_nbr_sim)
end

function _draw_param_kernel!(vec_p_prime::Vector{Float64}, 
                             vec_p_star::SubArray{Float64,1}, 
                             mat_p_old::Matrix{Float64}, wl_old::Vector{Float64},
                             mat_cov::Matrix{Float64}, 
                             tree_mat_p::Union{KDTree,Nothing}, kernel_type::String)
    if kernel_type == "mvnormal"
        d_mvnormal = MvNormal(vec_p_star, mat_cov)
        rand!(d_mvnormal, vec_p_prime)
    elseif kernel_type == "knn_mvnormal"
        k = Int(round(0.25 * size(mat_p_old)[2]))
        idxs, dist = knn(tree_mat_p, vec_p_star, k, true)
        knn_mat_cov = 2 * cov(mat_p_old[:,idxs], ProbabilityWeights(wl_old[idxs]), 2; corrected=false)
        d_knn_mvnormal = Distributions.MvNormal(vec_p_star, knn_mat_cov)
        rand!(d_knn_mvnormal, vec_p_prime)
        return knn_mat_cov
    else
        error("Unknown specified kernel")
    end
end

function _update_param!(mat_p::Matrix{Float64}, vec_dist::Vector{Float64}, 
                        l_nbr_sim::Vector{Int}, wl_current::Vector{Float64}, 
                        idx::Int,
                        pm::ParametricModel, epsilon::Float64,
                        wl_old::Vector{Float64}, mat_p_old::Matrix{Float64},
                        mat_cov::Matrix{Float64}, tree_mat_p::Union{Nothing,KDTree}, 
                        kernel_type::String, sym_var_aut::VariableAutomaton;
                        l_obs::Union{Nothing,AbstractVector} = nothing, 
                        func_dist::Union{Nothing,Function} = nothing)
    d_weights = Categorical(wl_old)
    dist_sim = Inf
    nbr_sim = 0
    vec_p_prime = zeros(pm.df)
    while dist_sim > epsilon
        ind_p_star = rand(d_weights)
        vec_p_star = view(mat_p_old, :, ind_p_star)
        knn_mat_cov = _draw_param_kernel!(vec_p_prime, vec_p_star, mat_p_old, wl_old, mat_cov, tree_mat_p, kernel_type)
        if !insupport(pm, vec_p_prime)
            continue
        end
        if l_obs == nothing
            S = volatile_simulate(pm, vec_p_prime)
            dist_sim = S[sym_var_aut]
        else
            l_sim = [simulate(pm, vec_p_prime) for i = 1:length(l_obs)]
            dist_sim = func_dist(l_sim, l_obs)
        end
        nbr_sim += 1
    end
    
    if kernel_type == "mvnormal"
        mat_cov_kernel = mat_cov 
    elseif kernel_type == "knn_mvnormal"
        mat_cov_kernel = knn_mat_cov
    else
        error("Unknown kernel")
    end
    # Update
    mat_p[:,idx] = vec_p_prime
    vec_dist[idx] = dist_sim
    l_nbr_sim[idx] = nbr_sim
    _update_weight!(wl_current, idx, pm, wl_old, mat_p_old, vec_p_prime, mat_cov_kernel, kernel_type)
end

function _update_weight!(wl_current::Vector{Float64}, idx::Int,
                         pm::ParametricModel,
                         wl_old::Vector{Float64}, mat_p_old::Matrix{Float64}, 
                         vec_p_prime::Vector{Float64},
                         mat_cov_kernel::Matrix{Float64}, kernel_type::String) 
    denom = 0.0
    for j in 1:length(wl_old)
        #denom += wl_old[j] * Distributions.pdf(d_normal, inv_sqrt_mat_cov * (vec_p_current - mat_p_old[:,j]))::Float64 
        denom += wl_old[j] * pdf(MvNormal(mat_p_old[:,j], mat_cov_kernel), vec_p_prime)::Float64 
    end
    wl_current[idx] = prior_pdf(pm, vec_p_prime) / denom 
end

function effective_sample_size(wl::AbstractVector{Float64})
    n_eff = 0.0
    wls = sum(wl)
    if wls > 0.0
        n_eff = wls ^ 2 / sum(wl .^ 2)
    end
    return n_eff
end

