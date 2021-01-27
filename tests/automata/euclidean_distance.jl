
using MarkovProcesses
import LinearAlgebra: dot
import Distributions: Uniform

load_automaton("euclidean_distance_automaton")
load_model("SIR")
load_model("ER")

test_all = true

# SIR model
nbr_sim = 10
for i = 1:nbr_sim
    set_param!(SIR, [:ki, :kr], [rand(Uniform(5E-5, 3E-3)), rand(Uniform(5E-3, 0.2))]) 
    let tml_obs, y_obs, sync_SIR, σ, test
        tml_obs = rand(Uniform(0.0, 5.0)):1.0:rand(Uniform(50.0, 100.0))
        y_obs = vectorize(simulate(SIR), :I, tml_obs)
        sync_SIR = SIR * create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
        #@show test, euclidean_distance(σ, tml_obs, y_obs, :I), σ.state_lha_end[:d]
        global test_all = test_all && test
    end
end

# ER model
for i = 1:nbr_sim
    let tml_obs, y_obs, sync_SIR, σ, test
        set_param!(ER, :k3, rand(Uniform(0.0, 100.0))) 
        tml_obs = rand(Uniform(0.0, 0.2)):1.0:rand(Uniform(0.5,10.0))
        y_obs = vectorize(simulate(ER), :P, tml_obs)
        sync_ER = ER * create_euclidean_distance_automaton(ER, tml_obs, y_obs, :P)
        σ = simulate(sync_ER)
        test = euclidean_distance(σ, :P, tml_obs, y_obs) == σ.state_lha_end[:d]
        #@show test, euclidean_distance(σ, tml_obs, y_obs, :P), σ.state_lha_end[:d]
        global test_all = test_all && test
    end
end

return test_all

