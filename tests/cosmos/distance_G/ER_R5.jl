
using MarkovProcesses
absolute_path = get_module_path() * "/tests/cosmos/"

# Remarque: la valeur estimée par Cosmos varie plus que de 0.01

# Values x1, x2  t1, t2
str_model = "ER"
load_model(str_model)
model = ER
ER.buffer_size = 100
load_automaton("automaton_G")

test_all = true
width = 0.1
level = 0.95
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
A_G = create_automaton_G(model, x1, x2, t1, t2, "P")  

l_k1 = 0.0:0.2:1.4
l_k1 = 0.2:0.2
l_k2 = 0:20:100
l_k2 = 40:40
nb_k1 = length(l_k1)
nb_k2 = length(l_k2)
mat_dist_cosmos = zeros(nb_k1,nb_k2)
mat_dist_pkg = zeros(nb_k1,nb_k2)
mat_full_k1 = zeros(nb_k1,nb_k2)
mat_full_k2 = zeros(nb_k1,nb_k2)
for i = 1:nb_k1
    for j = 1:nb_k2
        # Cosmos estimation
        k1 = l_k1[i]
        k2 = l_k2[j]
        command = `Cosmos $(absolute_path * "models/" * str_model * ".gspn") 
        $(absolute_path * "distance_G/dist_G_"  * str_model * ".lha") --njob $(ENV["JULIA_NUM_THREADS"]) 
        --const k_1=$(k1),k2=$(k2),x1=$x1,x2=$x2,t1=$t1,t2=$t2 
        --level $(level) --width $(width) 
        --verbose 2` 
        #run(pipeline(command, stderr=devnull))
        @timev run(pipeline(command))
        dict_values = cosmos_get_values("Result_dist_G_$(str_model).res")
        mat_dist_cosmos[i,j] = dict_values["Estimated value"]
        nb_sim = dict_values["Total paths"]
        nb_accepted = dict_values["Accepted paths"]
        nb_sim = convert(Int, nb_sim)
        # MarkovProcesses estimation
        set_param!(ER, "k1", convert(Float64, k1))
        set_param!(ER, "k2", convert(Float64, k2))
        sync_ER = ER*A_G
        mat_dist_ij_threads = zeros(Threads.nthreads()) 
        @timev Threads.@threads for j = 1:nb_sim
            σ = simulate(sync_ER)
            mat_dist_ij_threads[Threads.threadid()] += (σ.S)["d"]
        end
        mat_dist_pkg[i,j] = sum(mat_dist_ij_threads) / nb_sim
        global test_all = test_all && isapprox(mat_dist_cosmos[i,j], mat_dist_pkg[i,j]; atol = width)
    end
end

@show mat_dist_pkg
@show mat_dist_cosmos

rm("Result_dist_G_$(str_model).res")
rm("Result.res")

return test_all

