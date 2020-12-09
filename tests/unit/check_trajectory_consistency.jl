
using MarkovProcesses

load_model("SIR")
load_model("ER")
load_automaton("automaton_F")

test_all = true
nb_sim = 1000

function show_traj(io::IO, σ::AbstractTrajectory, m::Model)
    println(io, "length(σ.values[1]), length(σ.times), length(σ.transitions)")
    println(io, "$(length(σ.values[1])), $(length(σ.times)), $(length(σ.transitions))")
    println(io, "isbounded(m), isbounded(σ)")
    println(io, "$(isbounded(m)), $(isbounded(σ))")
    println(io, σ.values)
    println(io, "$(σ.values)")
    println(io, σ.times)
    println(io, "$(σ.times)")
    println(io, "σ.transitions")
    println(io, "$(σ.transitions)")
end

σ_dump = simulate(SIR)
σ2_dump = simulate(ER)

for i = 1:nb_sim
    local σ = simulate(SIR)
    local σ2 = simulate(ER)
    try
        global test_all = test_all && check_consistency(σ) && check_consistency(σ2) &&
                          !isbounded(σ) && !isbounded(σ2)
    catch err
        show_traj(stdout, σ, SIR)
        show_traj(stdout, σ2, ER)
        global σ_dump = σ
        global σ2_dump = σ2
        throw(err)
    end
    if !test_all
        show_traj(stdout, σ, SIR)
        show_traj(stdout, σ2, ER)
        global σ_dump = σ
        global σ2_dump = σ2
        error("Ouch")
    end
end

new_SIR = deepcopy(SIR)
sync_SIR = new_SIR * create_automaton_F(new_SIR, 0.0, Inf, 100.0, 110.0, :I)
new_ER = deepcopy(ER)
sync_ER = new_ER * create_automaton_F(new_ER, 0.0, 100.0, 4.0, 5.0, :P)
for i = 1:nb_sim
    local σ = simulate(sync_SIR)
    local σ2 = simulate(sync_ER)
    try
        global test_all = test_all && check_consistency(σ) && check_consistency(σ2) &&
                          !isbounded(σ) && !isbounded(σ2)
    catch err
        show_traj(stdout, σ, sync_SIR)
        show_traj(stdout, σ2, sync_ER)
        global σ_dump = σ
        global σ2_dump = σ2
        throw(err)
    end
    if !test_all
        show_traj(stdout, σ, sync_SIR)
        show_traj(stdout, σ2, sync_ER)
        global σ_dump = σ
        global σ2_dump = σ2
        error("Ouch")
    end
end

SIR.time_bound = 1.0
ER.time_bound = 0.01
for i = 1:nb_sim
    local σ = simulate(SIR)
    local σ2 = simulate(ER)
    try
        global test_all = test_all && check_consistency(σ) && check_consistency(σ2) &&
                          isbounded(σ) && isbounded(σ2)
    catch err
        show_traj(stdout, σ, SIR)
        show_traj(stdout, σ2, ER)
        global σ_dump = σ
        global σ2_dump = σ2
        throw(err)
    end
    if !test_all
        show_traj(stdout, σ, SIR)
        show_traj(stdout, σ2, ER)
        global σ_dump = σ
        global σ2_dump = σ2
        error("Ouch")
    end
end

SIR.time_bound = Inf
ER.time_bound = Inf
new_x0_SIR = [95, 0, 5]
new_x0_ER = [0, 0, 0, 100]
SIR.x0 = new_x0_SIR
ER.x0 = new_x0_ER
σ = simulate(SIR)
σ2 = simulate(ER)
test_all = test_all && SIR.isabsorbing(SIR.p,new_x0_SIR) && ER.isabsorbing(ER.p,new_x0_ER) &&
                    length_states(σ) == 1 && length_states(σ2) == 1 &&
                    check_consistency(σ) && check_consistency(σ2)

SIR.time_bound = 1.0
ER.time_bound = 0.01
σ = simulate(SIR)
σ2 = simulate(ER)
test_all = test_all && SIR.isabsorbing(SIR.p,new_x0_SIR) && ER.isabsorbing(ER.p,new_x0_ER) &&
                    length_states(σ) == 2 && length_states(σ2) == 2 &&
                    check_consistency(σ) && check_consistency(σ2)

return test_all

