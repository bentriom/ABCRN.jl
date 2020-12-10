
using MarkovProcesses 

model_SIR = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end "SIR"
set_x0!(model_SIR, [95,5,0])
set_param!(model_SIR, [0.012, 0.05])
simulate(model_SIR)

model_unnamed_SIR = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end
set_x0!(model_unnamed_SIR, [95,5,0])
set_param!(model_unnamed_SIR, [0.012, 0.05])
simulate(model_unnamed_SIR)

model_ER = @network_model begin
    R1: (E+S => ES, k1*E*S)
    R2: (ES => E+S, k2*ES)
    R3: (ES => E+P, k3*ES)
end "ER"
set_x0!(model_ER, [100,100,0,0])
set_param!(model_ER, [1.0,1.0,1.0])
simulate(model_ER)

return true

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

