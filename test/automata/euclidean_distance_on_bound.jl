
using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

MAKE_SECOND_AUTOMATON_TESTS = false

load_model("SIR")
load_model("repressilator")
load_automaton("euclidean_distance_automaton")
load_automaton("abc_euclidean_distance_automaton")
if MAKE_SECOND_AUTOMATON_TESTS
    load_automaton("euclidean_distance_automaton_2")
end

tml_obs = 0:0.5:200
set_time_bound!(SIR, 200.0)
y_obs_sir = vectorize(simulate(SIR), :I, tml_obs)
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs_sir, :I)
sync_SIR = SIR * aut1
set_time_bound!(repressilator, 200.0)
y_obs_repr = vectorize(simulate(repressilator), :P1, tml_obs)
aut1 = create_euclidean_distance_automaton(repressilator, tml_obs, y_obs_repr, :P1)
aut1_abc = create_abc_euclidean_distance_automaton(repressilator, tml_obs, y_obs_repr, :P1)
sync_repressilator = repressilator * aut1
sync_abc_repressilator = repressilator * aut1_abc

test_all = true
nbr_sim = 10

for i = 1:nbr_sim
    let σ, test
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs_sir) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :I, tml_obs, y_obs_sir), σ.state_lha_end[:d]
            @show σ
        end
        σ = simulate(sync_repressilator)
        test = test && euclidean_distance(σ, :P1, tml_obs, y_obs_repr) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :P1, tml_obs, y_obs_repr), σ.state_lha_end[:d]
            #@show σ
        end
        σ = simulate(sync_abc_repressilator)
        test = test && euclidean_distance(σ, :P1, tml_obs, y_obs_repr) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :P1, tml_obs, y_obs_repr), σ.state_lha_end[:d]
            #@show σ
        end
        global test_all = test_all && test
    end
end

if MAKE_SECOND_AUTOMATON_TESTS
    aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs_sir, :I)
    sync_SIR = SIR * aut2
    for i = 1:nbr_sim
        let σ, test2
            σ = simulate(sync_SIR)
            test2 = euclidean_distance(σ, :I, tml_obs, y_obs_sir) == σ.state_lha_end[:d]
            if !test2
                @show euclidean_distance(σ, :I, tml_obs, y_obs_sir), σ.state_lha_end[:d]
                @show σ
            end
            global test_all = test_all && test2
        end
    end
end

tml_obs = 0:20.0:200
set_time_bound!(SIR, 200.0)
y_obs_sir = vectorize(simulate(SIR), :I, tml_obs)
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs_sir, :I)
sync_SIR = SIR * aut1
set_time_bound!(repressilator, 200.0)
y_obs_repr = vectorize(simulate(repressilator), :P1, tml_obs)
aut1 = create_euclidean_distance_automaton(repressilator, tml_obs, y_obs_repr, :P1)
aut1_abc = create_abc_euclidean_distance_automaton(repressilator, tml_obs, y_obs_repr, :P1)
sync_repressilator = repressilator * aut1
sync_abc_repressilator = repressilator * aut1_abc

test_all = true
nbr_sim = 10

for i = 1:nbr_sim
    let σ, test
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs_sir) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :I, tml_obs, y_obs_sir), σ.state_lha_end[:d]
            @show σ
        end
        σ = simulate(sync_repressilator)
        test = test && euclidean_distance(σ, :P1, tml_obs, y_obs_repr) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :P1, tml_obs, y_obs_repr), σ.state_lha_end[:d]
            #@show σ
        end
        σ = simulate(sync_abc_repressilator)
        test = test && euclidean_distance(σ, :P1, tml_obs, y_obs_repr) == σ.state_lha_end[:d]
        if !test
            @show tml_obs
            @show euclidean_distance(σ, :P1, tml_obs, y_obs_repr), σ.state_lha_end[:d]
            #@show σ
        end
        global test_all = test_all && test
    end
end

if MAKE_SECOND_AUTOMATON_TESTS
    aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs_sir, :I)
    sync_SIR = SIR * aut2
    for i = 1:nbr_sim
        let σ, test2
            σ = simulate(sync_SIR)
            test2 = euclidean_distance(σ, :I, tml_obs, y_obs_sir) == σ.state_lha_end[:d]
            if !test2
                @show euclidean_distance(σ, :I, tml_obs, y_obs_sir), σ.state_lha_end[:d]
                @show σ
            end
            global test_all = test_all && test2
        end
    end
end

return test_all

