
load_automaton(automaton::String) = include(get_module_path() * "/automata/$(automaton).jl")

length_var(A::LHA) = length(A.map_var_automaton_idx)
get_value(A::LHA, x::SubArray{Int,1}, var::String) = x[A.map_var_model_idx[var]]

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
        print(io, "* $var = $(S.l_var[idx])\n")
    end
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

function next_state!(Snplus1::StateLHA, A::LHA, 
                    xnplus1::SubArray{Int,1}, tnplus1::Float64, tr_nplus1::Transition, 
                    Sn::StateLHA; verbose = false)
    for i in eachindex(Snplus1.l_var)
        Snplus1.l_var[i] += (A.l_flow[Sn.loc])[i]*(tnplus1 - Sn.time) 
    end
    Snplus1.time = tnplus1

    # En fait d'apres observation de Cosmos, après qu'on ait lu la transition on devrait stop.
    edge_candidates = Edge[]
    tuple_candidates = Tuple{Location, Location}[]
    first_round = true
    detected_event = (tr_nplus1 == nothing) ? true : false
    turns = 1
    while first_round || !detected_event || length(edge_candidates) > 0
        edge_candidates = Edge[]
        tuple_candidates = Tuple{Location,Location}[]
        current_loc = Snplus1.loc
        if verbose
            @show turns
        end
        # Save all edges that satisfies transition predicate (synchronous or autonomous)
        for loc in A.l_loc
            tuple_edges = (current_loc, loc)
            if haskey(A.map_edges, tuple_edges)
                for edge in A.map_edges[tuple_edges]
                    if edge.check_constraints(A, Snplus1)
                        push!(edge_candidates, edge)
                        push!(tuple_candidates, tuple_edges)
                    end
                end
            end
        end
        # Search the one we must chose
        ind_edge = 0
        for (i, edge) in enumerate(edge_candidates)
            if edge.transitions[1] == nothing
                ind_edge = i
                break
            end
            if first_round || !detected_event
                if (length(edge.transitions) == 1 && tr_nplus1 != nothing && edge.transitions[1] == "ALL") || 
                    (tr_nplus1 in edge.transitions)
                    detected_event = true
                    ind_edge = i
                end
            end
        end
        # Update the state with the chosen one (if it exists)
        if ind_edge > 0
            edge_candidates[ind_edge].update_state!(A, Snplus1, xnplus1)
            # Should add something like if edges_candidates[ind_edge].transition != nohting break end ??
        end
        if verbose
            @show first_round, detected_event
            @show tnplus1, tr_nplus1, xnplus1
            @show ind_edge
            @show edge_candidates
            @show tuple_candidates
            @show Snplus1
        end
        if ind_edge == 0 break end
        # For debug
        if turns > 10
            println("Turns, Bad behavior of region2 automaton")
            @show first_round, detected_event
            @show length(edge_candidates)
            @show tnplus1, tr_nplus1, xnplus1
            @show edge_candidates
            for edge in edge_candidates
                @show edge.check_constraints(A, Snplus1)
            end
            error("Unpredicted behavior automaton F v2")
        end
        turns += 1
        first_round = false
    end
end

# For debug purposes
function read_trajectory(A::LHA, σ::Trajectory; verbose = false)
    A_new = LHA(A, (σ.m)._map_obs_var_idx)
    l_t = times(σ)
    l_tr = transitions(σ)
    mat_x = zeros(Int, length_states(σ), σ.m.d)
    for (i,var) in enumerate(σ.m.g)
        mat_x[:,i] = σ[var]
    end
    Sn = init_state(A_new, σ[1], l_t[1])
    Snplus1 = copy(Sn)
    if verbose println("Init: ") end
    if verbose @show Sn end
    for n in 1:length_states(σ)
        xn = view(mat_x, n, :)
        next_state!(Snplus1, A_new, xn, l_t[n], l_tr[n], Sn; verbose = verbose)
        if Snplus1.loc in A_new.l_loc_final 
            break 
        end
        Sn = Snplus1 
    end
    return Sn
end

