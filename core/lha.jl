
length_var(A::LHA) = length(getfield(A, :map_var_automaton_idx))
get_value(S::StateLHA, x::Vector{Int}, var::VariableModel) = x[getfield(getfield(S, :A), :map_var_model_idx)[var]]
get_value(A::LHA, x::Vector{Int}, var::VariableModel) = x[getfield(A, :map_var_model_idx)[var]]

copy(S::StateLHA) = StateLHA(getfield(S, :A), getfield(S, :loc), getfield(S, :values), getfield(S, :time))
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
    print(io, "$(A.name) automaton (LHA)\n")
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
    print(io, "(Edge: ")
    print(io, (E.transitions[1] == nothing) ? "Asynchronous #, " : ("Synchronized with " * join(E.transitions,',') * ", "))
    print(io, Symbol(E.check_constraints))
    print(io, ", ")
    print(io, Symbol(E.update_state!))
    print(io, ")")
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

# A push! method implementend by myself because I know in most cases the edge candidates
# are not greater than 2
function _push_edge!(edge_candidates::Vector{Edge}, edge::Edge, nbr_candidates::Int)
    if nbr_candidates < 2
        edge_candidates[nbr_candidates+1] = edge
    else
        push!(edge_candidates, edge)
    end
end

function _find_edge_candidates!(edge_candidates::Vector{Edge},
                                edges_from_current_loc::Dict{Location,Vector{Edge}},
                                Λ::Dict{Location,Function},
                                Snplus1::StateLHA,
                                x::Vector{Int}, p::Vector{Float64},
                                only_asynchronous::Bool)
    nbr_candidates = 0
    for target_loc in keys(edges_from_current_loc)
        for edge in edges_from_current_loc[target_loc]
            if Λ[target_loc](x) && getfield(edge, :check_constraints)(Snplus1, x, p)
                if getfield(edge, :transitions)[1] == nothing
                    _push_edge!(edge_candidates, edge, nbr_candidates)
                    nbr_candidates += 1
                    return nbr_candidates
                end
                if !only_asynchronous && getfield(edge, :transitions)[1] != nothing
                    _push_edge!(edge_candidates, edge, nbr_candidates)
                    nbr_candidates += 1
                end
            end
        end
    end
    return nbr_candidates
end

function _get_edge_index(edge_candidates::Vector{Edge}, nbr_candidates::Int,
                         detected_event::Bool, tr_nplus1::Transition)
    ind_edge = 0
    bool_event = detected_event
    for i = 1:nbr_candidates
        edge = edge_candidates[i]
        # Asynchronous edge detection: we fire it
        if getfield(edge, :transitions)[1] == nothing 
            return (i, bool_event)
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

function next_state!(Snplus1::StateLHA, A::LHA, 
                     xnplus1::Vector{Int}, tnplus1::Float64, tr_nplus1::Transition, 
                     Sn::StateLHA, xn::Vector{Int}, p::Vector{Float64}; verbose::Bool = false)
    # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
    edge_candidates = Vector{Edge}(undef, 2)
    first_round::Bool = true
    detected_event::Bool = false
    turns = 0
       
    if verbose 
        println("##### Begin next_state!")
        @show xnplus1, tnplus1, tr_nplus1
        @show Sn 
        @show Snplus1 
    end
    # In terms of values not reference, Snplus1 == Sn
    # First, we check the asynchronous transitions
    while first_round || length(edge_candidates) > 0
        turns += 1
        #edge_candidates = empty!(edge_candidates)
        current_loc = getfield(Snplus1, :loc)
        edges_from_current_loc = getfield(A, :map_edges)[current_loc]
        # Save all edges that satisfies transition predicate (asynchronous ones)
        nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, getfield(A, :Λ), Sn, xn, p, true) #Snplus1=>Sn
        # Search the one we must chose, here the event is nothing because 
        # we're not processing yet the next event
        ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, nothing)
        # Update the state with the chosen one (if it exists)
        # Should be xn here
        if ind_edge > 0
            getfield(edge_candidates[ind_edge], :update_state!)(Snplus1, xn, p)
        end
        first_round = false
        if verbose
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event, nbr_candidates
            println("After update")
            @show Snplus1
        end
        if (ind_edge == 0) 
            break 
        end
        # For debug
        #=
        if turns > 100
            println("Number of turns in next_state! is suspicious")
            @show first_round, detected_event
            @show length(edge_candidates)
            @show tnplus1, tr_nplus1, xnplus1
            @show edge_candidates
            for edge in edge_candidates
                @show getfield(edge, :check_constraints)(Snplus1, xn, p)
            end
            error("Unpredicted behavior automaton")
        end
        =#
    end
    if verbose 
        @show Snplus1 
        println("Time flies with the flow...")
    end
    # Now time flies according to the flow
    for i in eachindex(Snplus1.values)
        coeff_deriv = (A.flow[Snplus1.loc])[i]
        if coeff_deriv > 0
            Snplus1.values[i] += coeff_deriv*(tnplus1 - Snplus1.time)
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
        current_loc = getfield(Snplus1, :loc)
        edges_from_current_loc = getfield(A, :map_edges)[current_loc]
        # Save all edges that satisfies transition predicate (synchronous ones)
        nbr_candidates = _find_edge_candidates!(edge_candidates, edges_from_current_loc, getfield(A, :Λ), Snplus1, xnplus1, p, false)
        # Search the one we must chose
        ind_edge, detected_event = _get_edge_index(edge_candidates, nbr_candidates, detected_event, tr_nplus1)
        # Update the state with the chosen one (if it exists)
        if verbose 
            @show turns
            @show edge_candidates
            @show ind_edge, detected_event, nbr_candidates
        end
        if ind_edge > 0
            getfield(edge_candidates[ind_edge], :update_state!)(Snplus1, xnplus1, p)
        end
        first_round = false
        if verbose
            println("After update")
            @show detected_event
            @show Snplus1
        end
        if (ind_edge == 0 || detected_event) 
            break 
        end
        # For debug
        #=
        if turns > 100
            println("Number of turns in next_state! is suspicious")
            @show detected_event
            @show length(edge_candidates)
            @show tnplus1, tr_nplus1, xnplus1
            @show edge_candidates
            for edge in edge_candidates
                @show getfield(edge, :check_constraints)(Snplus1, x, p)
            end
            error("Unpredicted behavior automaton")
        end
        =#
    end
    if verbose println("##### End next_state!") end
end

# For tests purposes
function read_trajectory(A::LHA, σ::Trajectory; verbose = false)
    @assert (σ.m).dim_state == σ.m.dim_obs_state # Model should be entirely obserbed 
    A_new = LHA(A, (σ.m)._map_obs_var_idx)
    p_sim = (σ.m).p
    l_t = times(σ)
    l_tr = transitions(σ)
    Sn = init_state(A_new, σ[1], l_t[1])
    Snplus1 = copy(Sn)
    if verbose println("Init: ") end
    if verbose @show Sn end
    for n in 2:length_states(σ)
        next_state!(Snplus1, A_new, σ[n], l_t[n], l_tr[n], Sn, σ[n-1], p_sim; verbose = verbose)
        copyto!(Sn, Snplus1)
        if Snplus1.loc in A_new.locations_final 
            break 
        end
    end
    return Sn
end

