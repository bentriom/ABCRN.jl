
import Plots: plot, plot!, scatter!
import Plots: current, palette, display, png, close

"""
    `plot(σ, var...; plot_transitions=false)`

Plot a simulated trajectory σ. var... is a tuple of stirng variables.
`plot(σ)` will plot all the variables simulated in σ 
whereas `plot(σ, "I", "R")` only plots the variables I and R of the trajectory (if it exists).
If `plot_transitions=true`, a marker that corresponds to a transition of the model will be plotted
at each break of the trajectory.
"""
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

