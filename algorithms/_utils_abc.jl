
import StatsBase: cov, ProbabilityWeights
import Statistics: quantile
import NearestNeighbors: knn, KDTree
import Distributions: MvNormal, Categorical
import Random: rand!

function _init_abc_automaton!(old_mat_p::Matrix{Float64}, vec_dist::Vector{Float64}, 
                              pm::ParametricModel, str_var_aut::String)
    vec_p = zeros(pm.df)
    for i = eachindex(vec_dist)
        draw!(vec_p, pm)
        old_mat_p[:,i] = vec_p
        σ = simulate(pm, vec_p)
        vec_dist[i] = σ.S[str_var_aut]
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

function _draw_param_kernel!(vec_p_prime::Vector{Float64}, 
                             vec_p_star::SubArray{Float64,1}, old_mat_p::Matrix{Float64}, wl_old::Vector{Float64},
                             mat_cov::Matrix{Float64}, tree_mat_p::Union{KDTree,Nothing}, kernel_type::String)
    if kernel_type == "mvnormal"
        d_mvnormal = MvNormal(vec_p_star, mat_cov)
        rand!(d_mvnormal, vec_p_prime)
    elseif kernel_type == "knn_mvnormal"
        k = Int(round(0.25 * size(old_mat_p)[2]))
        idxs, dist = knn(tree_mat_p, vec_p_star, k, true)
        knn_mat_cov = 2 * cov(old_mat_p[:,idxs], ProbabilityWeights(wl_old[idxs]), 2; corrected=false)
        d_knn_mvnormal = Distributions.MvNormal(vec_p_star, knn_mat_cov)
        rand!(d_knn_mvnormal, vec_p_prime)
        return knn_mat_cov
    else
        error("Unknown specified kernel")
    end
end

function _update_param!(mat_p::Matrix{Float64}, vec_dist::Vector{Float64}, 
                       l_nbr_sim::Vector{Int}, wl_current::Vector{Float64}, idx::Int,
                       pm::ParametricModel, epsilon::Float64,
                       wl_old::Vector{Float64}, old_mat_p::Matrix{Float64},
                       mat_cov::Matrix{Float64}, tree_mat_p::Union{Nothing,KDTree}, 
                       kernel_type::String, str_var_aut::String)
    d_weights = Categorical(wl_old)
    dist_sim = Inf
    nbr_sim = 0
    vec_p_prime = zeros(pm.df)
    while dist_sim > epsilon
        ind_p_star = rand(d_weights)
        vec_p_star = view(old_mat_p, :, ind_p_star)
        knn_mat_cov = _draw_param_kernel!(vec_p_prime, vec_p_star, old_mat_p, wl_old, mat_cov, tree_mat_p, kernel_type)
        if !insupport(pm, vec_p_prime)
            continue
        end
        σ = simulate(pm, vec_p_prime)
        dist_sim = σ.S[str_var_aut]
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
    _update_weight!(wl_current, idx, pm, wl_old, old_mat_p, vec_p_prime, mat_cov_kernel, kernel_type)
end

function _update_weight!(wl_current::Vector{Float64}, idx::Int,
                         pm::ParametricModel,
                         wl_old::Vector{Float64}, old_mat_p::Matrix{Float64}, 
                         vec_p_prime::Vector{Float64},
                         mat_cov_kernel::Matrix{Float64}, 
                         kernel_type::String) 
    denom = 0.0
    for j in 1:length(wl_old)
        #denom += wl_old[j] * Distributions.pdf(d_normal, inv_sqrt_mat_cov * (vec_p_current - old_mat_p[:,j]))::Float64 
        denom += wl_old[j] * pdf(MvNormal(old_mat_p[:,j], mat_cov_kernel), vec_p_prime)::Float64 
    end
    wl_current[idx] = prior_pdf(pm, vec_p_prime) / denom 
end

function effective_sample_size(wl::Vector{Float64})
    n_eff = 0.0
    wls = sum(wl)
    if wls > 0.0
        n_eff = wls ^ 2 / sum(wl .^ 2)
    end
    return n_eff
end

