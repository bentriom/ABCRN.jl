
length_var(A::LHA) = length(getfield(A, :map_var_automaton_idx))
get_value(S::StateLHA, x::Vector{Int}, var::VariableModel) = x[getfield(getfield(S, :A), :map_var_model_idx)[var]]
get_value(A::LHA, x::Vector{Int}, var::VariableModel) = x[getfield(A, :map_var_model_idx)[var]]

copy(S::StateLHA) = StateLHA(getfield(S, :A), getfield(S, :loc), copy(getfield(S, :values)), getfield(S, :time))
# Not overring getproperty, setproperty to avoid a conversion Symbol => String for the dict key

# From the variable automaton var symbol this function get the index in S.values
get_idx_var_automaton(S::StateLHA, var::VariableAutomaton) = getfield(getfield(S, :A), :map_var_automaton_idx)[var]
get_value(S::StateLHA, idx_var::Int) = getfield(S, :values)[idx_var]
get_state_lha_value(A::LHA, values::Vector{Float64}, var::VariableAutomaton) = values[A.map_var_automaton_idx[var]]

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
        Sdest.values[i] = Ssrc.values[i]
    end
    Sdest.time = Ssrc.time
end

# In future check_consistency(LHA), check if constant has the name
# of one of the LHA fields

isaccepted(S::StateLHA) = (getfield(S, :loc) in getfield(getfield(S, :A), :locations_final))
isaccepted(loc::Symbol, A::LHA) = loc in A.locations_final

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

function generate_code_next_state_with_dicts_lha(lha_name::Symbol, edge_type::Symbol)
    
    return quote 
        # A push! method implementend by myself because of preallocation of edge_candidates
        function _push_edge!(edge_id_candidates::Vector{Int}, target_loc_candidates::Vector{Symbol}, 
                             edge_id::Int, target_loc::Symbol, nbr_candidates::Int)
            if nbr_candidates < length(edge_id_candidates)
                 edge_id_candidates[nbr_candidates+1] = edge_id
                 target_loc_candidates[nbr_candidates+1] = target_loc
            else
                push!(edge_id_candidates, edge_id)
                push!(target_loc_candidates, target_loc)
            end
        end

        function _find_edge_candidates!(edge_id_candidates::Vector{Int}, target_loc_candidates::Vector{Symbol},
                                        dict_transitions_from_current_loc::Dict{Location,Vector{TransitionSet}},
                                        dict_check_constraints_from_current_loc::Dict{Location,Vector{Function}},
                                        Λ::Dict{Location,Function},
                                        S_time::Float64, S_values::Vector{Float64},
                                        x::Vector{Int}, p::Vector{Float64},
                                        only_asynchronous::Bool)
            nbr_candidates = 0
            for target_loc in keys(dict_transitions_from_current_loc)
                if !Λ[target_loc](x) continue end
                for i = eachindex(dict_transitions_from_current_loc[target_loc])
                    check_constraints_edge = dict_check_constraints_from_current_loc[target_loc][i]
                    if check_constraints_edge(S_time, S_values, x, p)
                        transitions = dict_transitions_from_current_loc[target_loc][i]
                        if transitions == nothing
                            _push_edge!(edge_id_candidates, target_loc_candidates, i, target_loc, nbr_candidates)
                            nbr_candidates += 1
                            return nbr_candidates
                        else
                            if !only_asynchronous
                                _push_edge!(edge_id_candidates, target_loc_candidates, i, target_loc, nbr_candidates)
                                nbr_candidates += 1
                            end
                        end
                    end
                end
            end
            return nbr_candidates
        end

        function _get_edge_index(edge_id_candidates::Vector{Int}, target_loc_candidates::Vector{Symbol}, nbr_candidates::Int,
                                 dict_transitions_from_current_loc::Dict{Location,Vector{TransitionSet}},
                                 detected_event::Bool, tr_nplus1::Transition)
            ind_edge = 0
            bool_event = detected_event
            for i = 1:nbr_candidates
                target_loc = target_loc_candidates[i]
                edge_id = edge_id_candidates[i]
                transitions = dict_transitions_from_current_loc[target_loc][edge_id]
                # Asynchronous edge detection: we fire it
                if transitions == nothing
                    return (i, detected_event)
                end
                # Synchronous detection
                if !detected_event && tr_nplus1 != nothing
                    if (transitions[1] == :ALL) || (tr_nplus1 in transitions)
                        ind_edge = i
                        bool_event = true
                    end
                end
            end
            return (ind_edge, bool_event)
        end

        function next_state!(A::$(lha_name),
                             ptr_loc_state::Vector{Symbol}, values_state::Vector{Float64}, ptr_time_state::Vector{Float64},
                             xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                             xn::Vector{Int}, p::Vector{Float64},
                             edge_id_candidates::Vector{Int}, target_loc_candidates::Vector{Symbol}; verbose::Bool = false)
            # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
            detected_event::Bool = false
            turns = 0
            Λ = getfield(A, :Λ)
            flow = getfield(A, :flow)
            map_edges = A.map_edges
            map_edges_transitions = A.map_edges_transitions
            map_edges_check_constraints = A.map_edges_check_constraints
            map_edges_update_state = A.map_edges_update_state
            if verbose 
                println("##### Begin next_state!")
                @show xnplus1, tnplus1, tr_nplus1
            end
            # First, we check the asynchronous transitions
            while true
                turns += 1
                if verbose @show turns end
                #edge_candidates = empty!(edge_candidates) 
                dict_transitions_from_current_loc = map_edges_transitions[ptr_loc_state[1]]
                dict_check_constraints_from_current_loc = map_edges_check_constraints[ptr_loc_state[1]]
                # Save all edges that satisfies transition predicate (asynchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_id_candidates, target_loc_candidates, 
                                                        dict_transitions_from_current_loc, dict_check_constraints_from_current_loc,
                                                        Λ, ptr_time_state[1], values_state, xn, p, true)
                # Search the one we must chose, here the event is nothing because 
                # we're not processing yet the next event
                ind_edge, detected_event = _get_edge_index(edge_id_candidates, target_loc_candidates, nbr_candidates, 
                                                           dict_transitions_from_current_loc,
                                                           detected_event, nothing)
                # Update the state with the chosen one (if it exists)
                # Should be xn here
                #first_round = false
                if ind_edge > 0
                    edge_target_loc = target_loc_candidates[ind_edge]
                    edge_id = edge_id_candidates[ind_edge]
                    firing_update_state! = map_edges_update_state[ptr_loc_state[1]][edge_target_loc][edge_id]
                    ptr_loc_state[1] = firing_update_state!(ptr_time_state[1], values_state, xn, p)
                else
                    if verbose 
                        println("No edge fired:")
                        @show ind_edge, detected_event, nbr_candidates
                    end
                    break 
                end
                if verbose
                    println("Edge fired:")
                    @show edge_id_candidates, target_loc_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show ptr_loc_state[1]
                    @show ptr_time_state[1]
                    @show values_state
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
            for i in eachindex(values_state)
                 coeff_deriv = flow[ptr_loc_state[1]][i]
                if coeff_deriv > 0
                      values_state[i] += coeff_deriv*(tnplus1 - ptr_time_state[1])
                end
            end
            ptr_time_state[1] = tnplus1
            if verbose 
                @show ptr_loc_state[1]
                @show ptr_time_state[1]
                @show values_state
            end
            # Now firing an edge according to the event 
            while true
                turns += 1
                if verbose @show turns end
                edges_from_current_loc = map_edges[ptr_loc_state[1]]
                dict_transitions_from_current_loc = map_edges_transitions[ptr_loc_state[1]]
                dict_check_constraints_from_current_loc = map_edges_check_constraints[ptr_loc_state[1]]
                # Save all edges that satisfies transition predicate (synchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_id_candidates, target_loc_candidates,
                                                        dict_transitions_from_current_loc, dict_check_constraints_from_current_loc,
                                                        Λ, ptr_time_state[1], values_state, xnplus1, p, false)
                # Search the one we must chose
                ind_edge, detected_event = _get_edge_index(edge_id_candidates, target_loc_candidates, nbr_candidates, 
                                                           dict_transitions_from_current_loc,
                                                           detected_event, tr_nplus1)
                # Update the state with the chosen one (if it exists)
                if ind_edge > 0
                    edge_target_loc = target_loc_candidates[ind_edge]
                    edge_id = edge_id_candidates[ind_edge]
                    firing_update_state! = map_edges_update_state[ptr_loc_state[1]][edge_target_loc][edge_id]
                    ptr_loc_state[1] = firing_update_state!(ptr_time_state[1], values_state, xnplus1, p)
                end
                if ind_edge == 0 || detected_event
                    if verbose 
                        if detected_event    
                            println("Synchronized with $(tr_nplus1)") 
                            @show edge_id_candidates, target_loc_candidates
                            @show ind_edge, detected_event, nbr_candidates
                            @show detected_event
                            @show ptr_loc_state[1]
                            @show ptr_time_state[1]
                            @show values_state
                        else
                            println("No edge fired")
                        end
                    end
                    break 
                end
                if verbose
                    @show edge_id_candidates, target_loc_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show detected_event
                    @show ptr_loc_state[1]
                    @show ptr_time_state[1]
                    @show values_state
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
            if verbose 
                println("##### End next_state!") 
            end
        end
    end
end

############################################################################################################

function generate_code_next_state(lha_name::Symbol, edge_type::Symbol)
    
    return quote 
        # A push! method implementend by myself because of preallocation of edge_candidates
        function _push_edge!(edge_candidates::Vector{<:$(edge_type)}, edge::$(edge_type), nbr_candidates::Int)
            if nbr_candidates < length(edge_candidates)
                 edge_candidates[nbr_candidates+1] = edge
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
                    if edge.check_constraints(S_time, S_values, x, p)
                        if edge.transitions == nothing
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
                if edge.transitions == nothing
                    return (i, detected_event)
                end
                # Synchronous detection
                if !detected_event && tr_nplus1 != nothing
                    if (edge.transitions[1] == :ALL) || (tr_nplus1 in edge.transitions)
                        ind_edge = i
                        bool_event = true
                    end
                end
            end
            return (ind_edge, bool_event)
        end

        function next_state!(A::$(lha_name),
                             ptr_loc_state::Vector{Symbol}, values_state::Vector{Float64}, ptr_time_state::Vector{Float64},
                             xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                             xn::Vector{Int}, p::Vector{Float64},
                             edge_candidates::Vector{$(edge_type)}; verbose::Bool = false)
            # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
            detected_event::Bool = false
            turns = 0
            Λ = getfield(A, :Λ)
            flow = getfield(A, :flow)
            map_edges = getfield(A, :map_edges)

            if verbose 
                println("##### Begin next_state!")
                @show xnplus1, tnplus1, tr_nplus1
            end
            # First, we check the asynchronous transitions
            while true
                turns += 1
                #edge_candidates = empty!(edge_candidates) 
                edges_from_current_loc = map_edges[ptr_loc_state[1]]
                # Save all edges that satisfies transition predicate (asynchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                        ptr_time_state[1], values_state, xn, p, true)
                # Search the one we must chose, here the event is nothing because 
                # we're not processing yet the next event
                ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, nothing)
                # Update the state with the chosen one (if it exists)
                # Should be xn here
                #first_round = false
                if ind_edge > 0
                    firing_edge = edge_candidates[ind_edge]
                    ptr_loc_state[1] = firing_edge.update_state!(ptr_time_state[1], values_state, xn, p)
                else
                    if verbose println("No edge fired") end
                    break 
                end
                if verbose
                    @show turns
                    @show edge_candidates
                    @show ind_edge, detected_event, nbr_candidates
                    @show ptr_loc_state[1]
                    @show ptr_time_state[1]
                    @show values_state
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
            for i in eachindex(values_state)
                 coeff_deriv = flow[ptr_loc_state[1]][i]
                if coeff_deriv > 0
                      values_state[i] += coeff_deriv*(tnplus1 - ptr_time_state[1])
                end
            end
            ptr_time_state[1] = tnplus1
            if verbose 
                @show ptr_loc_state[1]
                @show ptr_time_state[1]
                @show values_state
            end
            # Now firing an edge according to the event 
            while true
                turns += 1
                edges_from_current_loc = map_edges[ptr_loc_state[1]]
                # Save all edges that satisfies transition predicate (synchronous ones)
                nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, Λ, 
                                                        ptr_time_state[1], values_state, xnplus1, p, false)
                # Search the one we must chose
                ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, tr_nplus1)
                # Update the state with the chosen one (if it exists)
                if ind_edge > 0
                    firing_edge = edge_candidates[ind_edge]
                    ptr_loc_state[1] = firing_edge.update_state!(ptr_time_state[1], values_state, xnplus1, p)
                end
                if ind_edge == 0 || detected_event
                    if verbose 
                        if detected_event    
                            println("Synchronized with $(tr_nplus1)") 
                            @show turns
                            @show edge_candidates
                            @show ind_edge, detected_event, nbr_candidates
                            @show detected_event
                            @show ptr_loc_state[1]
                            @show ptr_time_state[1]
                            @show values_state
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
                    @show ptr_loc_state[1]
                    @show ptr_time_state[1]
                    @show values_state
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
            if verbose 
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
    S = init_state(A, σ[1], l_t[1])
    ptr_loc = [S.loc]
    ptr_time = [S.time]
    values = S.values
    lha_name = Symbol(typeof(A))
    edge_type = getfield(Main, Symbol("Edge$(lha_name)"))
    edge_candidates = Vector{edge_type}(undef, 2)
    if verbose println("Init: ") end
    if verbose @show S end
    func_next_state! = getfield(Main, :next_state!) # It's defined after the compilation of the package
    for k in 2:length_states(σ)
        func_next_state!(A, ptr_loc, values, ptr_time, σ[k], l_t[k], l_tr[k], σ[k-1], p_sim, edge_candidates; verbose = verbose)
        if ptr_loc[1] in A_new.locations_final 
            break 
        end
    end
    setfield!(S, :loc, ptr_loc[1])
    setfield!(S, :time, ptr_time[1])
    return S
end

