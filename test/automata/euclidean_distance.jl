
using BiochemNetABC
import LinearAlgebra: dot
import Distributions: Uniform

load_automaton("euclidean_distance_automaton")
load_automaton("euclidean_distance_automaton_2")
load_automaton("abc_euclidean_distance_automaton")
load_model("SIR")
load_model("ER")
observe_all!(SIR)
observe_all!(ER)

MAKE_SECOND_AUTOMATON_TESTS = false
test_all = true

# SIR model
nbr_sim = 10
for i = 1:nbr_sim
    set_param!(SIR, [:ki, :kr], [rand(Uniform(5E-5, 3E-3)), rand(Uniform(5E-3, 0.2))]) 
    let tml_obs, y_obs, sync_SIR, σ, test, test2
        tml_obs = rand(Uniform(0.0, 5.0)):1.0:rand(Uniform(50.0, 100.0))
        y_obs = vectorize(simulate(SIR), :I, tml_obs)
        sync_SIR = SIR * create_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
        σ = simulate(sync_SIR)
        test = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test
            @show test, euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
            global err = σ
            global tml = tml_obs
            global y = y_obs
            global sync_model = sync_SIR
            break
        end
        if MAKE_SECOND_AUTOMATON_TESTS
            sync_SIR = SIR * create_euclidean_distance_automaton_2(SIR, tml_obs, y_obs, :I)
            σ = simulate(sync_SIR)
            test2 = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
            if !test2
                @show test2, euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
                global err = σ
                global tml = tml_obs
                global y = y_obs
                global sync_model = sync_SIR
                break
            end
        else
            test2 = true
        end
        sync_SIR = SIR * create_abc_euclidean_distance_automaton(SIR, tml_obs, y_obs, :I)
        σ = simulate(sync_SIR)
        test3 = euclidean_distance(σ, :I, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test3
            @show test3, euclidean_distance(σ, :I, tml_obs, y_obs), σ.state_lha_end[:d]
            global err = σ
            global tml = tml_obs
            global y = y_obs
            global sync_model = sync_SIR
            break
        end
        global test_all = test_all && test && test2 && test3
    end
end

# ER model
for i = 1:nbr_sim
    let tml_obs, y_obs, sync_SIR, σ, test, test2
        set_param!(ER, :k3, rand(Uniform(0.0, 100.0))) 
        tml_obs = rand(Uniform(0.0, 0.2)):1.0:rand(Uniform(0.5,10.0))
        y_obs = vectorize(simulate(ER), :P, tml_obs)
        sync_ER = ER * create_euclidean_distance_automaton(ER, tml_obs, y_obs, :P)
        σ = simulate(sync_ER)
        test = euclidean_distance(σ, :P, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test
            @show test, euclidean_distance(σ, :P, tml_obs, y_obs), σ.state_lha_end[:d]
            global err = σ
            global tml = tml_obs
            global y = y_obs
            global sync_model = sync_ER
            break
        end
        if MAKE_SECOND_AUTOMATON_TESTS
            sync_ER = ER * create_euclidean_distance_automaton_2(ER, tml_obs, y_obs, :P)
            σ = simulate(sync_ER)
            test2 = euclidean_distance(σ, :P, tml_obs, y_obs) == σ.state_lha_end[:d]
            if !test2
                @show test2, euclidean_distance(σ, :P, tml_obs, y_obs), σ.state_lha_end[:d]
                global err = σ
                global tml = tml_obs
                global y = y_obs
                global sync_model = sync_ER
                break
            end
        else
            test2 = true
        end
        sync_ER = ER * create_abc_euclidean_distance_automaton(ER, tml_obs, y_obs, :P)
        σ = simulate(sync_ER)
        test3 = euclidean_distance(σ, :P, tml_obs, y_obs) == σ.state_lha_end[:d]
        if !test3
            @show test3, euclidean_distance(σ, :P, tml_obs, y_obs), σ.state_lha_end[:d]
            global err = σ
            global tml = tml_obs
            global y = y_obs
            global sync_model = sync_ER
            break
        end
        global test_all = test_all && test && test2 && test3
    end
end

return test_all

