
using ABCRN 

test_all = true

model_SIR = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end "SIR"
set_x0!(model_SIR, [95,5,0])
set_param!(model_SIR, [0.012, 0.05])
simulate(model_SIR)

test_all = test_all && keys(model_SIR.map_var_idx) == Set([:S, :I, :R]) &&
                       keys(model_SIR.map_param_idx) == Set([:ki, :kr])


model_unnamed_SIR = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end
set_x0!(model_unnamed_SIR, [95,5,0])
set_param!(model_unnamed_SIR, [0.012, 0.05])
simulate(model_unnamed_SIR)

test_all = test_all && keys(model_unnamed_SIR.map_var_idx) == Set([:S, :I, :R]) &&
                       keys(model_unnamed_SIR.map_param_idx) == Set([:ki, :kr])

model_ER = @network_model begin
    R1: (E+S => ES, k1*E*S)
    R2: (ES => E+S, k2*ES)
    R3: (ES => E+P, k3*ES)
end "ER"
set_x0!(model_ER, [100,100,0,0])
set_param!(model_ER, [1.0,1.0,1.0])
simulate(model_ER)

test_all = test_all && keys(model_ER.map_var_idx) == Set([:E, :S, :ES, :P]) &&
                       keys(model_ER.map_param_idx) == Set([:k1, :k2, :k3])

model_birth_death = @network_model begin
    Birth: (X => 2X, λ*X)
    Death: (X => ∅, μ*X)
end "Birth-death process"

test_all = test_all && keys(model_birth_death.map_var_idx) == Set([:X]) &&
                       keys(model_birth_death.map_param_idx) == Set([:λ, :μ])

model_birth_death_2 = @network_model begin
    Birth: (X => 2X, λ*X)
    Death: (X => 0, μ*X)
end "Birth-death process"

test_all = test_all && keys(model_birth_death_2.map_var_idx) == Set([:X]) &&
                       keys(model_birth_death_2.map_param_idx) == Set([:λ, :μ])

return test_all

#=
@network_model "test1" begin
    R1: (S+I => 2I+Z,  ki*S*I)
    R2: (I => R, kr*I)
end

@network_model "test2" begin
    R1: (S+I => 2I,  ki*S*I)
    R2: (I => R, kr*I)
end

@network_model "test3" begin
    R1: (S+I => I,  ki)
    R2: (2I => R, kr*I*c)
end
=#

