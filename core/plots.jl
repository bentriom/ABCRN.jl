
import Plots: plot, plot!, scatter!, hline!, Shape, text
import Plots: current, palette, display, png, close

"""
    `plot(σ, var...; plot_transitions=false)`

Plot a simulated trajectory σ. var... is a tuple of stirng variables.
`plot(σ)` will plot all the variables simulated in σ 
whereas `plot(σ, :I, :R)` only plots the variables I and R of the trajectory (if it exists).
If `plot_transitions=true`, a marker that corresponds to a transition of the model will be plotted
at each break of the trajectory.
"""
function plot(σ::AbstractTrajectory, vars::VariableModel...; plot_transitions::Bool = false, A::Union{Nothing,LHA} = nothing, filename::String = "")
    # Setup 
    palette_tr = palette(:default)
    l_tr = unique(transitions(σ))
    map_tr_color(tr) = palette_tr[findfirst(x->x==tr, l_tr)]
    to_plot = vars
    if length(vars) ==  0
        to_plot = get_obs_var(σ)
    end
    
    # Plots
    p = plot(title = "Trajectory of $(σ.m.name) model", palette = :lightrainbow, legend = :outertopright, background_color_legend=:transparent, dpi = 480)
    for var in to_plot
        @assert var in get_obs_var(σ) "Variable $var is not observed." 
        plot!(p, times(σ), σ[var], 
              xlabel = "Time", ylabel = "Number of species",
              label = var,
              linetype=:steppost)
    end
    if plot_transitions
        for (i, var) in enumerate(to_plot)
            for tr in l_tr
                idx_tr = findall(x->x==tr, transitions(σ))
                label = (tr == nothing || i > 1) ? "" : tr
                alpha = (tr == nothing) ? 0.0 : 0.5
                scatter!(p, times(σ)[idx_tr], σ[var][idx_tr], label=label, 
                         markershape=:cross, markeralpha=alpha, 
                         markersize = 2,
                         markercolor=palette_tr[findfirst(x->x==tr, l_tr)])
            end
        end
    end
    if A != nothing
        plot!(A)
    end
    if filename == ""
        display(p)
    else
        png(p, filename)
    end
end

function plot!(A::LHA)
    if A.name in ["F property", "G property"]
        if haskey(A.constants, :x1) && haskey(A.constants, :x2) && haskey(A.constants, :t1) && haskey(A.constants, :t2) 
            x1, x2, t1, t2 = A.constants[:x1], A.constants[:x2], A.constants[:t1], A.constants[:t2] 
            plot!(Shape([(t1,x1), (t1,x2), (t2,x2), (t2,x1), (t1,x1)]), opacity = 0.5, label = "Region LHA")
        end
    end
    if A.name in ["G and F property"]    
        if haskey(A.constants, :x1) && haskey(A.constants, :x2) && haskey(A.constants, :t1) && haskey(A.constants, :t2) 
            x1, x2, t1, t2 = A.constants[:x1], A.constants[:x2], A.constants[:t1], A.constants[:t2] 
            plot!(Shape([(t1,x1), (t1,x2), (t2,x2), (t2,x1), (t1,x1)]), opacity = 0.5, label = "Region LHA")
        end
        if haskey(A.constants, :x3) && haskey(A.constants, :x4) && haskey(A.constants, :t3) && haskey(A.constants, :t4) 
            x3, x4, t3, t4 = A.constants[:x3], A.constants[:x4], A.constants[:t3], A.constants[:t4] 
            plot!(Shape([(t3,x3), (t3,x4), (t4,x4), (t4,x3), (t3,x3)]), opacity = 0.5, label = "Region LHA")
        end
    end
    if A.name == "Period"
        hline!([A.constants[:L], A.constants[:H]], label = "L, H", linestyle = :dot)
    end
end

# For tests purposes
function plot_periodic_trajectory(A::LHA, σ::SynchronizedTrajectory, sym_obs::Symbol; verbose = false, filename::String = "")
    @assert sym_obs in get_obs_var(σ) "Variable is not observed in the model"
    @assert A.name in ["Period"]
    A_new = LHA(A, (σ.m)._map_obs_var_idx)
    p_sim = (σ.m).p
    l_t = times(σ)
    l_tr = transitions(σ)
    Sn = init_state(A_new, σ[1], l_t[1])
    Snplus1 = copy(Sn)
    nbr_states = length_states(σ)
    locations_trajectory = Vector{Location}(undef, nbr_states)
    locations_trajectory[1] = Sn.loc
    idx_n = [1]
    values_n = [Sn[:n]]
    if verbose println("Init: ") end
    if verbose @show Sn end
    for n in 2:nbr_states
        next_state!(Snplus1, A_new, σ[n], l_t[n], l_tr[n], Sn, σ[n-1], p_sim; verbose = verbose)
        copyto!(Sn, Snplus1)
        locations_trajectory[n] = Sn.loc
        if Sn[:n] != values_n[end]
            push!(idx_n, n)
            push!(values_n, Sn[:n])
        end
        if Snplus1.loc in A_new.locations_final 
            break 
        end
    end
    p = plot(title = "Oscillatory trajectory of $(σ.m.name) model", 
             background_color_legend=:transparent, dpi = 480, legend = :outertopright) #legendfontsize, legend
    colors_loc = Dict(:l0 => :silver, :l0prime => :silver, :final => :black,
                      :low => :skyblue2, :mid => :orange, :high => :red)
    loc_to_color(loc::Symbol) = colors_loc[loc]
    for loc in A.locations
        idx_locations = findall(x->x==loc, locations_trajectory)
        label_state = (loc in [:low, :mid, :high]) ? "$sym_obs ($loc loc)" : ""
        scatter!(p, times(σ)[idx_locations], σ[sym_obs][idx_locations], color = colors_loc[loc], 
                    markersize = 1.0, markershape = :cross, 
                    label = label_state, xlabel = "Time", ylabel = "Species $sym_obs")
    end
    annot_n = [(times(σ)[idx_n[i]], σ[sym_obs][idx_n[i]] - 10, text("n = $(values_n[i])", 6, :bottom)) for i = eachindex(idx_n)]
    scatter!(p, times(σ)[idx_n], σ[sym_obs][idx_n], annotations = annot_n,
                             markershape = :utriangle, markersize = 3, label = "n")
    hline!(p, [A.constants[:L], A.constants[:H]], label = "L, H", color = :grey; linestyle = :dot)
    
    if filename == ""
        display(p)
    else
        png(p, filename)
    end
end

export plot, plot!, plot_periodic_trajectory

