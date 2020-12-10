
@everywhere begin 
    using MarkovProcesses
    import Distributed: nworkers
    absolute_path = get_module_path() * "/tests/cosmos/"
    # Values x1, x2  t1, t2
    str_model = "ER"
    load_model(str_model)
    model = ER
    observe_all!(ER)
    ER.buffer_size = 100
    load_automaton("automaton_G")
    width = 0.2
    level = 0.95
    x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
    A_G = create_automaton_G(model, x1, x2, t1, t2, :E)  
    l_k1 = 0.0:0.5:1.5
    #l_k1 = 0.2:0.2
    l_k2 = 0:40:100
    #l_k2 = 40:40
end

test_all = true
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
        $(absolute_path * "distance_G/dist_G_"  * str_model * ".lha") --njob $(nworkers()) 
        --const k_1=$(k1),k_2=$(k2),x1=$x1,x2=$x2,t1=$t1,t2=$t2 
        --level $(level) --width $(width) 
        --verbose 0` 
        run(pipeline(command, stderr=devnull))
        dict_values = cosmos_get_values("Result_dist_G_$(str_model).res")
        mat_dist_cosmos[i,j] = dict_values["Estimated value"][1]
        nb_sim = dict_values["Total paths"][1]
        nb_accepted = dict_values["Accepted paths"][1]
        nb_sim = convert(Int, nb_sim)
        # MarkovProcesses estimation
        set_param!(ER, :k1, convert(Float64, k1))
        set_param!(ER, :k2, convert(Float64, k2))
        sync_ER = ER*A_G
        mat_dist_pkg[i,j], nb_accepts_pkg = distribute_mean_value_lha(sync_ER, :d, nb_sim; with_accepts = true)
        #@info "About accepts" nb_sim nb_accepted nb_accepts_pkg
        test = isapprox(mat_dist_cosmos[i,j], mat_dist_pkg[i,j]; atol = width*1.01)
        test2 = nb_accepts_pkg == (nb_sim / nb_accepted)
        if !test
            @info "Distances too different" (k1,k2) mat_dist_pkg[i,j] mat_dist_cosmos[i,j]
        end
        global test_all = test_all && test && test2 
    end
end

@info "Distances R5 pkg" mat_dist_pkg
@info "Distances R5 Cosmos" mat_dist_cosmos

rm("Result_dist_G_$(str_model).res")
rm("Result.res")

return test_all

