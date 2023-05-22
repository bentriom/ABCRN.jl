
using ABCRN

# Settings with two SIR models
load_model("SIR")
SIR2 = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end "SIR2"
set_x0!(SIR, [:S, :I, :R], [95, 5, 0])
set_x0!(SIR2, [:S, :I, :R], [95, 5, 0])
set_param!(SIR, [:ki, :kr], [0.012, 0.05])
set_param!(SIR2, [:ki, :kr], [0.012, 0.05])
# Settings ER model
load_model("ER")
set_x0!(ER, [:E, :S, :ES, :P], [100, 100, 0, 0])
set_param!(ER, [:k1, :k2, :k3], [1.0, 1.0, 1.0])
# Settings repressilator model
load_model("repressilator")
set_observed_var!(repressilator, [:mRNA1, :mRNA2, :mRNA3, :P1, :P2, :P3])
set_x0!(repressilator, [:mRNA1, :mRNA2, :mRNA3], fill(0, 3))
set_x0!(repressilator, [:P1, :P2, :P3], [5, 0, 15])
set_param!(repressilator, :n, 2.0)
set_param!(repressilator, [:α, :α0, :β, :n], [400.0, 0.0, 2.0, 2.0])
set_time_bound!(repressilator, 200.0)

test_all = true

# First test : infinite simulation
set_time_bound!(ER, Inf)
set_time_bound!(SIR, Inf)
set_time_bound!(SIR2, Inf)
for i = 1:10
    global test_all = test_all && simulate(ER).P[end] == 100
    global test_all = test_all && simulate(SIR).I[end] == 0
    global test_all = test_all && simulate(SIR2).I[end] == 0
end

# Second test: change the simulation stop criteria
get_idx(m::ContinuousTimeModel, var::Symbol) = m.map_var_idx[var]
idx_I_SIR, idx_I_SIR2, idx_P_ER, idx_P1_repressilator = 
get_idx(SIR, :I), get_idx(SIR2, :I), get_idx(ER, :P), get_idx(repressilator, :P1)

@eval new_isabsorbing_SIR(p::Vector{Float64}, x::Vector{Int}) = x[$(idx_I_SIR)] == 2
@eval new_isabsorbing_SIR2(p::Vector{Float64}, x::Vector{Int}) = x[$(idx_I_SIR2)] == 3
@eval new_isabsorbing_ER(p::Vector{Float64}, x::Vector{Int}) = x[$(idx_P_ER)] == 50
@eval new_isabsorbing_repressilator(p, x) = x[$(idx_P1_repressilator)] == 100

change_simulation_stop_criteria(SIR, :new_isabsorbing_SIR)
change_simulation_stop_criteria(SIR2, :new_isabsorbing_SIR2)
change_simulation_stop_criteria(ER, :new_isabsorbing_ER)
change_simulation_stop_criteria(repressilator, :new_isabsorbing_repressilator)

for i = 1:10
    global test_all = test_all && simulate(SIR).I[end] == 2
    global test_all = test_all && simulate(SIR2).I[end] == 3
    global test_all = test_all && simulate(ER).P[end] == 50
    global test_all = test_all && simulate(repressilator).P1[end] == 100
end

return test_all

