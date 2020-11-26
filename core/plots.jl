
import Plots: plot, plot!, scatter!
import Plots: palette, display, png, close

function plot(σ::AbstractTrajectory, vars::String...; filename::String = "", plot_transitions = false)
    # Setup 
    palette_tr = palette(:default)
    l_tr = unique(transitions(σ))
    map_tr_color(tr) = palette_tr[findfirst(x->x==tr, l_tr)]
    to_plot = vars
    if length(vars) ==  0
        to_plot = get_obs_var(σ)
    end
    
    # Plots
    p = plot(title = "Trajectory", palette = :lightrainbow)
    for var in to_plot
        @assert var in get_obs_var(σ) 
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
    if filename == ""
        display(p)
    else
        png(p, filename)
        close(p)
    end
end

export plot

