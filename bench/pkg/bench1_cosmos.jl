
@everywhere using MarkovProcesses
using Statistics
using DelimitedFiles
import Distributed: nworkers

absolute_path = get_module_path() * "/tests/cosmos/"
path_latex = "./"

dict_automata_lha = Dict("A_G" => (absolute_path * "distance_G/dist_G_ER.lha"))
dict_models_cosmos = Dict("ER" => absolute_path * "models/ER.gspn")

# Configuration 
k1, k2, k3 = 1.5, 1.0, 1.0
new_p = [k1, k2, k3]
width = 0.1
level = 0.99
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8

bench1_cosmos_time = zeros(0)
bench1_cosmos_mem = zeros(0)
bench1_pkg_time = zeros(0)
bench1_pkg_mem = zeros(0)
nb_sim = 0

# Bench 1 : Estim of d in
total_runs = 5
for nb_run = 1:total_runs
    command = `Cosmos $(absolute_path * "models/ER.gspn") 
    $(absolute_path * "distance_G/dist_G_ER.lha") --njob $(nworkers()) 
    --const k_1=$(k1),k_2=$(k2),k_3=$(k3),x1=$x1,x2=$x2,t1=$t1,t2=$t2 
    --level $(level) --width $(width) 
    --verbose 0` 
    run(pipeline(command, stderr=devnull))
    dict_values = cosmos_get_values("Result_dist_G_ER.res") 
    global nb_sim = convert(Int, dict_values["Total paths"][1])
    nb_accepted = dict_values["Accepted paths"][1]
    time_cosmos = dict_values["Total CPU time"][1]
    mem_cosmos = dict_values["Total Memory used"][1] * 10^6
    #@show time_cosmos, mem_cosmos
    push!(bench1_cosmos_time, time_cosmos)
    push!(bench1_cosmos_mem, mem_cosmos)
    rm("Result_dist_G_ER.res")
    rm("Result.res")

    @everywhere load_model("ER")
    observe_all!(ER)
    @everywhere load_automaton("automaton_G")
    A_G = create_automaton_G(ER, x1, x2, t1, t2, :E) 
    sync_ER = ER*A_G
    set_param!(ER, new_p)
    distribute_mean_value_lha(sync_ER, :d, 2)
    dist_pkg = @timed mean_value_lha(sync_ER, :d, nb_sim)
    push!(bench1_pkg_time, dist_pkg.time)
    push!(bench1_pkg_mem, dist_pkg.bytes)
    #@show dist_pkg.value, dict_values["Estimated value"][1]
    if (nb_run % 10 == 0) && (nb_run != 0)  println("$(div(total_runs,nb_run))0%"); end
end
println()

str_latex = "
\\begin{tabular}{|c|c|c|c|c|c|}
    \\hline
    Bench 1 & Mean time (s) & Max. time (s) &
Min. time (s) & \\begin{tabular}[c]{@{}c@{}}Mean\\\\Memory (MB)\\end{tabular} & Sim. \\\\
    \\hline
    Package & $(round(mean(bench1_pkg_time), digits=2))     & $(round(maximum(bench1_pkg_time), digits=2)) &
    $(round(minimum(bench1_pkg_time), digits=2))  & $(round(mean(bench1_pkg_mem)/(1024^2), digits=2)) & $nb_sim \\\\
    \\hline
    Cosmos  & $(round(mean(bench1_cosmos_time), digits=2))  & $(round(maximum(bench1_cosmos_time), digits=2)) &
    $(round(minimum(bench1_cosmos_time), digits=2)) & $(round(mean(bench1_cosmos_mem)/(1024^2), digits=2)) & $nb_sim \\\\
    \\hline
\\end{tabular}"

str_end_file = nworkers() > 1 ? "_distributed_$(nworkers())" : ""
open(path_latex * "bench1$(str_end_file).tex", "w+") do io
    write(io, str_latex)
end;
writedlm(path_latex * "values_bench1_cosmos$(str_end_file).csv", [bench1_cosmos_time bench1_cosmos_mem], ',')
writedlm(path_latex * "values_bench1_pkg$(str_end_file).csv", [bench1_pkg_time bench1_pkg_mem], ',')

@show mean(bench1_pkg_time), mean(bench1_cosmos_time)
@show mean(bench1_pkg_mem), mean(bench1_cosmos_mem)

