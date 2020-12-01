
length_var(A::LHA) = length(A.map_var_automaton_idx)
get_value(A::LHA, x::Vector{Int}, var::String) = x[A.map_var_model_idx[var]]

copy(S::StateLHA) = StateLHA(S.A, S.loc, S.l_var, S.time)
# Not overring getproperty, setproperty to avoid a conversion Symbol => String for the dict key
getindex(S::StateLHA, var::VariableAutomaton) = (S.l_var)[(S.A).map_var_automaton_idx[var]]
setindex!(S::StateLHA, val::Float64, var::VariableAutomaton) = (S.l_var)[(S.A).map_var_automaton_idx[var]] = val
setindex!(S::StateLHA, val::Int, var::VariableAutomaton) = (S.l_var)[(S.A).map_var_automaton_idx[var]] = convert(Float64, val)
setindex!(S::StateLHA, val::Bool, var::VariableAutomaton) = (S.l_var)[(S.A).map_var_automaton_idx[var]] = convert(Float64, val)

function Base.show(io::IO, S::StateLHA)
    print(io, "State of LHA\n")
    print(io, "- location: $(S.loc)\n")
    print(io, "- time: $(S.time)\n")
    print(io, "- variables:\n")
    for (var, idx) in (S.A).map_var_automaton_idx
        print(io, "* $var = $(S.l_var[idx]) (idx = $idx)\n")
    end
end

function Base.copyto!(Sdest::StateLHA, Ssrc::StateLHA)
    Sdest.A = Ssrc.A
    Sdest.loc = Ssrc.loc
    for i = eachindex(Sdest.l_var)
        Sdest.l_var[i] = Ssrc.l_var[i]
    end
    Sdest.time = Ssrc.time
end

isaccepted(S::StateLHA) = (S.loc in (S.A).l_loc_final)

# Methods for synchronize / read the trajectory
function init_state(A::LHA, x0::Vector{Int}, t0::Float64)
    S0 = StateLHA(A, "", zeros(length_var(A)), t0)
    for loc in A.l_loc_init
        if A.Λ[loc](A,S0) 
            S0.loc = loc
            break
        end
    end
    return S0
end

function _find_edge_candidates!(edge_candidates::Vector{Edge}, current_loc::Location, 
                                A::LHA, Snplus1::StateLHA, only_asynchronous::Bool)
    for loc in A.l_loc
        tuple_edges = (current_loc, loc)
        if haskey(A.map_edges, tuple_edges)
            for edge in A.map_edges[tuple_edges]
                if edge.check_constraints(A, Snplus1)
                    if edge.transitions[1] == nothing
                        pushfirst!(edge_candidates, edge)
                    end
                    if !only_asynchronous && edge.transitions[1] != nothing
                        push!(edge_candidates, edge)
                    end
                end
            end
        end
    end
end

function _get_edge_index(edge_candidates::Vector{Edge},
                         detected_event::Bool, tr_nplus1::Transition)
    ind_edge = 0
    bool_event = detected_event
    for i in eachindex(edge_candidates)
        edge = edge_candidates[i]
        # Asynchronous detection
        if edge.transitions[1] == nothing
            return (i, bool_event)
        end
        # Synchronous detection
        if !detected_event && tr_nplus1 != nothing
            if (length(edge.transitions) == 1 && edge.transitions[1] == "ALL") || 
               (tr_nplus1 in edge.transitions)
                return (i, true)
            end
        end
    end
    return (ind_edge, bool_event)
end

function next_state!(Snplus1::StateLHA, A::LHA, 
                    xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                    Sn::StateLHA; verbose::Bool = false)
    # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
    edge_candidates = Edge[]
    first_round::Bool = true
    detected_event::Bool = (tr_nplus1 == nothing) ? true : false
    turns = 0
       
    if verbose 
        println("#####")
        @show xnplus1, tnplus1, tr_nplus1
        @show Sn 
        @show Snplus1 
    end
   
    # In terms of values not reference, Snplus1 == Sn
    # First, we check the asynchronous transitions
    while first_round || length(edge_candidates) > 0
        turns += 1
        edge_candidates = empty!(edge_candidates)
        current_loc = Snplus1.loc
        # Save all edges that satisfies transition predicate (asynchronous ones)
        _find_edge_candidates!(edge_candidates, current_loc, A, Snplus1, true)
        # Search the one we must chose, here the event is nothing because 
        # we're not processing yet the next event
        ind_edge, detected_event = _get_edge_index(edge_candidates, detected_event, nothing)
        # Update the state with the chosen one (if it exists)
        if ind_edge > 0
            edge_candidates[ind_edge].update_state!(A, Snplus1, xnplus1)
        end
        first_round = false
        if verbose
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event
            println("After update")
            @show Snplus1
        end
        if (ind_edge == 0) 
            break 
        end
        # For debug
        if turns > 100
            println("Number of turns in next_state! is suspicious")
            @show first_round, detected_event
            @show length(edge_candidates)
            @show tnplus1, tr_nplus1, xnplus1
            @show edge_candidates
            for edge in edge_candidates
                @show edge.check_constraints(A, Snplus1)
            end
            error("Unpredicted behavior automaton")
        end
    end
    if verbose 
        @show Snplus1 
        println("Time flies with the flow...")
    end
    # Now time flies according to the flow
    for i in eachindex(Snplus1.l_var)
        coeff_deriv = (A.l_flow[Snplus1.loc])[i]
        if coeff_deriv > 0
            Snplus1.l_var[i] += coeff_deriv*(tnplus1 - Snplus1.time)
        end
    end
    Snplus1.time = tnplus1
    if verbose 
        @show Snplus1 
    end
    # Now firing an edge according to the event 
    first_round = true
    while first_round || length(edge_candidates) > 0
        turns += 1
        edge_candidates = empty!(edge_candidates)
        current_loc = Snplus1.loc
        # Save all edges that satisfies transition predicate (synchronous ones)
        _find_edge_candidates!(edge_candidates, current_loc, A, Snplus1, false)
        # Search the one we must chose
        ind_edge, detected_event = _get_edge_index(edge_candidates, detected_event, tr_nplus1)
        # Update the state with the chosen one (if it exists)
        if verbose 
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event
        end
        if ind_edge > 0
            edge_candidates[ind_edge].update_state!(A, Snplus1, xnplus1)
        end
        first_round = false
        if verbose
            println("After update")
            @show Snplus1
        end
        if (ind_edge == 0 || detected_event) 
            break 
        end
        # For debug
        if turns > 100
            println("Number of turns in next_state! is suspicious")
            @show detected_event
            @show length(edge_candidates)
            @show tnplus1, tr_nplus1, xnplus1
            @show edge_candidates
            for edge in edge_candidates
                @show edge.check_constraints(A, Snplus1)
            end
            error("Unpredicted behavior automaton")
        end
    end
end

# For tests purposes
function read_trajectory(A::LHA, σ::Trajectory; verbose = false)
    @assert σ.m.d == σ.m.dobs # Model should be entirely obserbed 
    A_new = LHA(A, (σ.m)._map_obs_var_idx)
    l_t = times(σ)
    l_tr = transitions(σ)
    Sn = init_state(A_new, σ[1], l_t[1])
    Snplus1 = copy(Sn)
    if verbose println("Init: ") end
    if verbose @show Sn end
    for n in 1:length_states(σ)
        next_state!(Snplus1, A_new, σ[n], l_t[n], l_tr[n], Sn; verbose = verbose)
        copyto!(Sn, Snplus1)
        if Snplus1.loc in A_new.l_loc_final 
            break 
        end
    end
    return Sn
end

