
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
    p = plot(title = "Trajectory of $(Symbol(typeof(σ.m)))", palette = :lightrainbow, legend = :outertopright, background_color_legend=:transparent, dpi = 480)
    for var in to_plot
        @assert var in get_obs_var(σ) "Variable $var is not observed." 
        plot!(p, times(σ), σ[var], 
              xlabel = "Time", ylabel = "Number of species",
              label = "$var",
              linetype=:steppost)
    end
    if plot_transitions
        for (i, var) in enumerate(to_plot)
            for tr in l_tr
                idx_tr = findall(x->x==tr, transitions(σ))
                label = (tr == nothing || i > 1) ? "" : "$tr"
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
        savefig(p, filename)
    end
end

function plot!(A::LHA; label::String = "")
    label_region = (label == "") ? "Region LHA" : label
    if (@isdefined(AutomatonF) && (typeof(A) <: AutomatonF)) || (@isdefined(AutomatonG) && (typeof(A) <: AutomatonG))
        x1, x2, t1, t2 = A.constants[:x1], A.constants[:x2], A.constants[:t1], A.constants[:t2] 
        plot!(Shape([(t1,x1), (t1,x2), (t2,x2), (t2,x1), (t1,x1)]), opacity = 0.5, label = label_region)
    end
    if @isdefined(AutomatonGandF) && typeof(A) <: AutomatonGandF
        x1, x2, t1, t2 = A.constants[:x1], A.constants[:x2], A.constants[:t1], A.constants[:t2] 
        plot!(Shape([(t1,x1), (t1,x2), (t2,x2), (t2,x1), (t1,x1)]), opacity = 0.5, label = label_region)
        x3, x4, t3, t4 = A.constants[:x3], A.constants[:x4], A.constants[:t3], A.constants[:t4] 
        plot!(Shape([(t3,x3), (t3,x4), (t4,x4), (t4,x3), (t3,x3)]), opacity = 0.5, label = label_region)
    end
    label_region = (label == "") ? "L, H" : label
    if @isdefined(PeriodAutomaton) && typeof(A) <: PeriodAutomaton
        hline!([A.constants[:L], A.constants[:H]], label = label_region, linestyle = :dot)
    end
end

# For tests purposes
function plot_periodic_trajectory(A::LHA, σ::SynchronizedTrajectory, sym_obs::Symbol; 
                                  verbose = false, annot_size = 6, show_tp = false, filename::String = "")
    @assert sym_obs in get_obs_var(σ) "Variable is not observed in the model"
    @assert typeof(A) <: PeriodAutomaton "The automaton is not a period automaton"
    @assert σ.m._g_idx == 1:σ.m.dim_state "All the variables must be observed. Call observe_all!(model) before."
    p_sim = (σ.m).p
    l_t = times(σ)
    l_tr = transitions(σ)
    S0 = init_state(A, σ[1], l_t[1])
    ptr_loc = [S0.loc]
    ptr_time = [S0.time]
    values = S0.values
    nbr_states = length_states(σ)
    locations_trajectory = Vector{Location}(undef, nbr_states)
    locations_trajectory[1] = S0.loc
    idx_n = [1]
    values_n = [MarkovProcesses.get_state_lha_value(A, values, :n)]
    values_tp = [MarkovProcesses.get_state_lha_value(A, values, :tp)]
    edge_candidates = Vector{EdgePeriodAutomaton}(undef, 2)
    if verbose println("Init: ") end
    for k in 2:nbr_states
        tp_k = MarkovProcesses.get_state_lha_value(A, values, :tp)
        next_state!(A, ptr_loc, values, ptr_time, σ[k], l_t[k], l_tr[k], σ[k-1], p_sim, edge_candidates; verbose = verbose)
        n_kplus1 = MarkovProcesses.get_state_lha_value(A, values, :n)
        if n_kplus1 != values_n[end]
            push!(idx_n, k)
            push!(values_n, n_kplus1)
            push!(values_tp, tp_k)
        end
        locations_trajectory[k] = ptr_loc[1]
        if ptr_loc[1] in A.locations_final 
            break
        end
    end
    p = plot(title = "Oscillatory trajectory of $(Symbol(typeof(σ.m)))", palette = :lightrainbow,
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
    annot_n = [(times(σ)[idx_n[i]], σ[sym_obs][idx_n[i]] - 20, text("n = $(values_n[i])", annot_size, :top)) for i = eachindex(idx_n)]
    annot_tp = [(times(σ)[idx_n[i]], σ[sym_obs][idx_n[i]] - 20, text("tp = $(round(values_tp[i], digits = 5))", annot_size, :bottom)) for i = eachindex(idx_n)]
    annots = (show_tp) ? vcat(annot_n, annot_tp) : annot_n
    scatter!(p, times(σ)[idx_n], σ[sym_obs][idx_n], annotations = annots,
                             markershape = :utriangle, markersize = 3, label = "n")
    hline!(p, [A.constants[:L], A.constants[:H]], label = "L, H", color = :grey; linestyle = :dot)
    if filename == ""
        display(p)
    else
        savefig(p, filename)
    end
end

export plot, plot!, plot_periodic_trajectory

