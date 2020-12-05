
using BenchmarkTools
using MarkovProcesses
using Catalyst
using DiffEqJump

load_model("ER")
set_param!(ER, "k1", 0.2)
set_param!(ER, "k2", 40.0)
ER.buffer_size = 100
ER.estim_min_states = 8000

b_pkg = @benchmark simulate(ER)
@show minimum(b_pkg), mean(b_pkg), maximum(b_pkg)

rs = @reaction_network begin
  c1, S + E --> SE
  c2, SE --> S + E
  c3, SE --> P + E
end c1 c2 c3
p = (0.2,40.0,1.0)   # [c1,c2,c3]
tspan = (0., 100.)
u0 = [100., 100., 0., 0.]  # [S,E,SE,P]

# solve JumpProblem
dprob = DiscreteProblem(rs, u0, tspan, p)
jprob = JumpProblem(rs, dprob, Direct())
jsol = solve(jprob, SSAStepper())

b_catalyst = @benchmark solve(jprob, SSAStepper())
@show minimum(b_catalyst), mean(b_catalyst), maximum(b_catalyst)

#plot(jsol,lw=2,title="Gillespie: Michaelis-Menten Enzyme Kinetics")

