
using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

MAKE_SECOND_AUTOMATON_TESTS = false

load_model("SIR")
tml_obs = 0:0.5:200
set_time_bound!(SIR, 200.0)
y_obs = vectorize(simulate(SIR), :I, tml_obs)

load_automaton("euclidean_distance_automaton")
aut1 = create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
sync_SIR = SIR * aut1
σ = simulate(sync_SIR)
test = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]

if MAKE_SECOND_AUTOMATON_TESTS
    load_automaton("euclidean_distance_automaton_2")
    aut2 = create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs, :I)
    sync_SIR = SIR * aut2
    σ = simulate(sync_SIR)
    test2 = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
else
    test2 = true
end

test_all = test && test2

return test_all

