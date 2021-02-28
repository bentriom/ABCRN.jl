
using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

MAKE_SECOND_AUTOMATON_TESTS = false

load_model("SIR")
load_automaton("euclidean_distance_automaton")
if MAKE_SECOND_AUTOMATON_TESTS
    load_automaton("euclidean_distance_automaton_2")
end

tml_obs = 0:0.5:200
set_time_bound!(SIR, 200.0)
y_obs = vectorize(simulate(SIR), :I, tml_obs)
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
sync_SIR = SIR * aut1

test_all = true
nbr_sim = 100

for i = 1:nbr_sim
    let σ, test
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
            @show σ
        end
        global test_all = test_all && test
    end
end

if MAKE_SECOND_AUTOMATON_TESTS
    aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs, :I)
    sync_SIR = SIR * aut2
    for i = 1:nbr_sim
        let σ, test2
            σ = simulate(sync_SIR)
            test2 = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
            if !test2
                @show euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
                @show σ
            end
            global test_all = test_all && test2
        end
    end
end

tml_obs = 0:20.0:200
set_time_bound!(SIR, 200.0)
y_obs = vectorize(simulate(SIR), :I, tml_obs)
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
sync_SIR = SIR * aut1

test_all = true
nbr_sim = 100

for i = 1:nbr_sim
    let σ, test
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
            @show σ
        end
        global test_all = test_all && test
    end
end

if MAKE_SECOND_AUTOMATON_TESTS
    aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs, :I)
    sync_SIR = SIR * aut2
    for i = 1:nbr_sim
        let σ, test2
            σ = simulate(sync_SIR)
            test2 = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
            if !test2
                @show euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
                @show σ
            end
            global test_all = test_all && test2
        end
    end
end

return test_all

