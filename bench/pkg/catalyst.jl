
using BenchmarkTools
using BiochemNetABC
using Catalyst
using DiffEqJump

path_latex = "./"

# Bench 1: A long simulation of ER

load_model("ER")
set_param!(ER, :k1, 0.2)
set_param!(ER, :k2, 40.0)
ER.buffer_size = 100
ER.estim_min_states = 8000

bench1_pkg = @benchmark simulate(ER)
@btime simulate(ER)
@show minimum(bench1_pkg), mean(bench1_pkg), maximum(bench1_pkg)

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

bench1_catalyst = @benchmark solve(jprob, SSAStepper())
@btime solve(jprob, SSAStepper())
@show minimum(bench1_catalyst), mean(bench1_catalyst), maximum(bench1_catalyst)

str_latex_bench1 = "
\\begin{tabular}{|c|c|c|c|c|}
    \\hline
    Bench 1 & Mean time (ms) & Max. time (ms) &
Min. time (ms) & \\begin{tabular}[c]{@{}c@{}}Mean\\\\Memory (KB)\\end{tabular} \\\\
    \\hline
    Package & $(round(time(mean(bench1_pkg))/1E6, digits=2))     & $(round(time(maximum(bench1_pkg))/1E6, digits=2)) &
    $(round(time(minimum(bench1_pkg))/1E6, digits=2))  & $(round(memory(mean(bench1_pkg))/1024, digits=2)) \\\\
    \\hline
    Catalyst.jl  & $(round(time(mean(bench1_catalyst))/1E6, digits=2))  & $(round(time(maximum(bench1_catalyst))/1E6, digits=2)) &
    $(round(time(minimum(bench1_catalyst))/1E6, digits=2)) & $(round(memory(mean(bench1_catalyst))/1024, digits=2)) \\\\
    \\hline
\\end{tabular}"

# Bench 2: Creation of model + a long simulation of ER

bench2_pkg = @benchmark begin
    ER = @network_model begin
        R1: (E+S => ES, k1*E*S)
        R2: (ES => E+S, k2*ES)
        R3: (ES => E+P, k3*ES)
    end "ER"
    set_param!(ER, :k1, 0.2)
    set_param!(ER, :k2, 40.0)
    ER.buffer_size = 100
    ER.estim_min_states = 8000
    simulate(ER)
end
@show minimum(bench2_pkg), mean(bench2_pkg), maximum(bench2_pkg)

bench2_catalyst = @benchmark begin
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
end
@show minimum(bench2_catalyst), mean(bench2_catalyst), maximum(bench2_catalyst)

str_latex_bench2 = "
\\begin{tabular}{|c|c|c|c|c|}
    \\hline
    Bench 1 & Mean time (ms) & Max. time (ms) &
Min. time (ms) & \\begin{tabular}[c]{@{}c@{}}Mean\\\\Memory (KB)\\end{tabular} \\\\
    \\hline
    Package & $(round(time(mean(bench2_pkg))/1E6, digits=2))     & $(round(time(maximum(bench2_pkg))/1E6, digits=2)) &
    $(round(time(minimum(bench2_pkg))/1E6, digits=2))  & $(round(memory(mean(bench2_pkg))/1024, digits=2)) \\\\
    \\hline
    Catalyst.jl  & $(round(time(mean(bench2_catalyst))/1E6, digits=2))  & $(round(time(maximum(bench2_catalyst))/1E6, digits=2)) &
    $(round(time(minimum(bench2_catalyst))/1E6, digits=2)) & $(round(memory(mean(bench2_catalyst))/1024, digits=2)) \\\\
    \\hline
\\end{tabular}"


open(path_latex * "bench1.tex", "w+") do io
    write(io, str_latex_bench1)
end;
open(path_latex * "bench2.tex", "w+") do io
    write(io, str_latex_bench2)
end;

