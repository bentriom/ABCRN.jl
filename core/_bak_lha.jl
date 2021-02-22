
length_var(A::LHA) = length(getfield(A, :map_var_automaton_idx))
get_value(S::StateLHA, x::Vector{Int}, var::VariableModel) = x[getfield(getfield(S, :A), :map_var_model_idx)[var]]
get_value(A::LHA, x::Vector{Int}, var::VariableModel) = x[getfield(A, :map_var_model_idx)[var]]

copy(S::StateLHA) = StateLHA(getfield(S, :A), getfield(S, :loc), copy(getfield(S, :values)), getfield(S, :time))
# Not overring getproperty, setproperty to avoid a conversion Symbol => String for the dict key

# From the variable automaton var symbol this function get the index in S.values
get_idx_var_automaton(S::StateLHA, var::VariableAutomaton) = getfield(getfield(S, :A), :map_var_automaton_idx)[var]
get_value(S::StateLHA, idx_var::Int) = getfield(S, :values)[idx_var]

getindex(S::StateLHA, var::VariableAutomaton) = getindex(getfield(S, :values), get_idx_var_automaton(S, var))
setindex!(S::StateLHA, val::Float64, var::VariableAutomaton) = setindex!(getfield(S, :values), val, get_idx_var_automaton(S, var))

getindex(S::StateLHA, l_var::Vector{VariableAutomaton}) = [S[var] for var in l_var]
setindex!(S::StateLHA, val::Int, var::VariableAutomaton) = S[var] = convert(Float64, val)
setindex!(S::StateLHA, val::Bool, var::VariableAutomaton) = S[var] = convert(Float64, val)

function Base.show(io::IO, A::LHA)
    print(io, "$(Symbol(typeof(A))) automaton (LHA)\n")
    print(io, "- initial locations : $(join(A.locations_init,','))\n")
    print(io, "- final locations : $(join(A.locations_final,','))\n")
    print(io, "- labeling prop Λ :\n")
    for loc  in keys(A.Λ)
        print(io, "* $loc: $(Symbol(A.Λ[loc]))\n")
    end
    print(io, "- variables :\n")
    for (var, idx) in A.map_var_automaton_idx
        print(io, "* $var (index = $idx in variables space)\n")
    end
    print(io, "- flow :\n")
    for (loc, flow_loc) in A.flow
        print(io, "* $loc: $flow_loc\n")
    end
    print(io, "- edges :\n")
    for from_loc in keys(A.map_edges)
        for to_loc in keys(A.map_edges[from_loc])
            edges = A.map_edges[from_loc][to_loc]
            print(io, "* $from_loc => $to_loc ($(length(edges))): $(join(edges,','))\n")
        end
    end
    print(io, "- constants :\n")
    for (name_constant, val_constant) in A.constants
        print(io, "* $name_constant: $val_constant\n")
    end
    print(io, "- transitions : $(join(A.transitions,','))\n")
end

function Base.show(io::IO, S::StateLHA)
    print(io, "State of LHA\n")
    print(io, "- location: $(S.loc)\n")
    print(io, "- variables:\n")
    for (var, idx) in (S.A).map_var_automaton_idx
        print(io, "* $var = $(S.values[idx]) (idx = $idx)\n")
    end
    print(io, "- time: $(S.time)")
end

function Base.show(io::IO, E::Edge)
    print(io, "($(typeof(E)): ")
    print(io, (E.transitions == nothing) ? "Asynchronous #)" : ("Synchronized with " * join(E.transitions,',') * ")"))
end

function Base.copyto!(Sdest::StateLHA, Ssrc::StateLHA)
    Sdest.A = Ssrc.A
    Sdest.loc = Ssrc.loc
    for i = eachindex(Sdest.values)
        @inbounds Sdest.values[i] = Ssrc.values[i]
    end
    Sdest.time = Ssrc.time
end

# In future check_consistency(LHA), check if constant has the name
# of one of the LHA fields

isaccepted(S::StateLHA) = (getfield(S, :loc) in getfield(getfield(S, :A), :locations_final))

# Methods for synchronize / read the trajectory
function init_state(A::LHA, x0::Vector{Int}, t0::Float64)
    S0 = StateLHA(A, :init, zeros(length_var(A)), t0)
    for loc in getfield(A, :locations_init)
        if A.Λ[loc](x0) 
            S0.loc = loc
            break
        end
    end
    return S0
end

function generate_code_next_state(lha_name::Symbol, edge_type::Symbol, 
                                  check_constraints::Symbol, update_state!::Symbol)
    
    return quote 
        # A push! method implementend by myself because of preallocation of edge_candidates
        function _push_edge!(edge_candidates::Vector{<:$(edge_type)}, edge::$(edge_type), nbr_candidates::Int)
            if nbr_candidates < 2
                @inbounds edge_candidates[nbr_candidates+1] = edge
            else
                push!(edge_candidates, edge)
            end
        end

        function _find_edge_candidates!(edge_candidates::Vector{$(edge_type)},
                                        edges_from_current_loc::Dict{Location,Vector{$(edge_type)}},
                                        Λ::Dict{Location,Function},
                                        S_time::Float64, S_values::Vector{Float64},
                                        x::Vector{Int}, p::Vector{Float64},
                                        only_asynchronous::Bool)
            nbr_candidates = 0
            for target_loc in keys(edges_from_current_loc)
                if !Λ[target_loc](x) continue end
                for edge in edges_from_current_loc[target_loc]
                    if $(check_constraints)(edge, S_time, S_values, x, p)
                        if getfield(edge, :transitions) == nothing
                            _push_edge!(edge_candidates, edge, nbr_candidates)
                            nbr_candidates += 1
                            return nbr_candidates
                        else
                            if !only_asynchronous
                                _push_edge!(edge_candidates, edge, nbr_candidates)
                                nbr_candidates += 1
                            end
                        end
                    end
                end
            end
            return nbr_candidates
        end

        function _get_edge_index(edge_candidates::Vector{$(edge_type)}, nbr_candidates::Int,
                                 detected_event::Bool, tr_nplus1::Transition)
            ind_edge = 0
            bool_event = detected_event
            for i = 1:nbr_candidates
                edge = edge_candidates[i]
                # Asynchronous edge detection: we fire it
                if getfield(edge, :transitions) == nothing
                    return (i, detected_event)
                end
                # Synchronous detection
                if !detected_event && tr_nplus1 != nothing
                    if (getfield(edge, :transitions)[1] == :ALL) || 
                        (tr_nplus1 in getfield(edge, :transitions))
                        ind_edge = i
                        bool_event = true
                    end
                end
            end
            return (ind_edge, bool_event)
        end

        function next_state!(Snplus1::StateLHA, A::$(lha_name),
                             xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                             Sn::StateLHA, xn::Vector{Int}, p::Vector{Float64},
                             edge_candidates::Vector{$(edge_type)}; verbose::Bool = false)
            # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
            detected_event::Bool = false
            turns = 0
            current_values = getfield(Snplus1, :values)
            current_time = getfield(Snplus1, :time)
            current_loc = getfield(Snplus1, :loc)
            Λ = getfield(A, :Λ)
            flow = getfield(A, :flow)
            map_edges = getfield(A, :map_edges)

            if verbose 
                println("##### Begin next_state!")
                @show xnplus1, tnplus1, tr_nplus1
                @show Sn 
            end
            # In terms of values not reference, Snplus1 == Sn
            # First, we check the asynchronous transitions
            while true
                turns += 1
                #edge_candidates = empty!(edge_candidates) 
                edges_from_current_loc = map_edges[current_loc]
                # Save all edges that satisfies transition predicate (asynchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                        current_time, current_values, xn, p, true)
                # Search the one we must chose, here the event is nothing because 
                # we're not processing yet the next event
                ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, nothing)
                # Update the state with the chosen one (if it exists)
                # Should be xn here
                #first_round = false
                if ind_edge > 0
                    firing_edge = edge_candidates[ind_edge]
                    current_loc = $(update_state!)(firing_edge, current_time, current_values, xn, p)
                else
                    if verbose println("No edge fired") end
                    break 
                end
                if verbose
                    @show turns
                    @show edge_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show current_loc
                    @show current_time
                    @show current_values
                    if turns == 500
                        @warn "We've reached 500 turns"
                    end
                end
                # For debug
                #=
                if turns > 100
                println("Number of turns in next_state! is suspicious")
                @show first_round, detected_event
                @show length(edge_candidates)
                @show tnplus1, tr_nplus1, xnplus1
                @show edge_candidates
                error("Unpredicted behavior automaton")
                end
                =#
            end
            if verbose 
                println("Time flies with the flow...")
            end
            # Now time flies according to the flow
            for i in eachindex(current_values)
                @inbounds coeff_deriv = flow[current_loc][i]
                if coeff_deriv > 0
                    @inbounds current_values[i] += coeff_deriv*(tnplus1 - current_time)
                end
            end
            current_time = tnplus1
            if verbose 
                @show current_loc
                @show current_time
                @show current_values
            end
            # Now firing an edge according to the event 
            while true
                turns += 1
                edges_from_current_loc = map_edges[current_loc]
                # Save all edges that satisfies transition predicate (synchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                        current_time, current_values, xnplus1, p, false)
                # Search the one we must chose
                ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, tr_nplus1)
                # Update the state with the chosen one (if it exists)
                if ind_edge > 0
                    firing_edge = edge_candidates[ind_edge]
                    current_loc = $(update_state!)(firing_edge, current_time, current_values, xnplus1, p)
                end
                if ind_edge == 0 || detected_event
                    if verbose 
                        if detected_event    
                            println("Synchronized with $(tr_nplus1)") 
                            @show turns
                            @show edge_candidates
                            @show ind_edge, detected_event, nbr_candidates
                            @show detected_event
                            @show current_loc
                            @show current_time
                            @show current_values
                        else
                            println("No edge fired")
                        end
                    end
                    break 
                end
                if verbose
                    @show turns
                    @show edge_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show detected_event
                    @show current_loc
                    @show current_time
                    @show current_values
                    if turns == 500
                        @warn "We've reached 500 turns"
                    end
                end
                # For debug
                #=
                if turns > 100
                println("Number of turns in next_state! is suspicious")
                @show detected_event
                @show length(edge_candidates)
                @show tnplus1, tr_nplus1, xnplus1
                @show edge_candidates
                error("Unpredicted behavior automaton")
                end
                =#
            end
            setfield!(Snplus1, :loc, current_loc)
            setfield!(Snplus1, :time, current_time)
            if verbose 
                @show Snplus1
                println("##### End next_state!") 
            end
        end
    end
end

# For tests purposes
function read_trajectory(A::LHA, σ::AbstractTrajectory; verbose = false)
    proba_model = get_proba_model(σ.m)
    @assert proba_model.dim_state == proba_model.dim_obs_state # Model should be entirely obserbed 
    A_new = LHA(A, proba_model._map_obs_var_idx)
    p_sim = proba_model.p
    l_t = times(σ)
    l_tr = transitions(σ)
    Sn = init_state(A_new, σ[1], l_t[1])
    Snplus1 = copy(Sn)
    edge_candidates = Vector{Edge}(undef, 2)
    if verbose println("Init: ") end
    if verbose @show Sn end
    for n in 2:length_states(σ)
        next_state!(Snplus1, A_new, σ[n], l_t[n], l_tr[n], Sn, σ[n-1], p_sim, edge_candidates; verbose = verbose)
        copyto!(Sn, Snplus1)
        if Snplus1.loc in A_new.locations_final 
            break 
        end
    end
    return Sn
end

