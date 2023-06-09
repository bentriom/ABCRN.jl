
@everywhere begin
    using BiochemNetABC
    import Distributed: nworkers
    absolute_path = get_module_path() * "/test/cosmos/"
    # Values x1, x2  t1, t2
    dict_exp = Dict(
                    :R1 => [50, 75, 0.025, 0.05],
                    :R2 => [50, 75, 0.05, 0.075],
                    :R3 => [25, 50, 0.05, 0.075]
                   )
    str_model = "ER"
    if str_model == "ER"
        load_model(str_model)
        model = ER
        observe_all!(ER)
    end
    load_automaton("automaton_F")
    test_all = true
    width = 0.2
    level = 0.95
    l_exp = [:R1, :R2, :R3]
end
#l_exp = [:R2]
for exp in l_exp
    if !(exp in keys(dict_exp))
        println("Unrecognized experiment: <<$exp>>")
        continue
    end
    val_exp = dict_exp[exp]
    l_k3 = 0:10:100
    nb_param = length(l_k3)
    l_dist_cosmos = zeros(nb_param)
    l_dist_pkg = zeros(nb_param)
    for i in 1:nb_param
        let test, test2, nb_sim, sync_ER, k3, x1, x2, t1, t2, A_F
            x1, x2, t1, t2 = val_exp[1], val_exp[2], val_exp[3], val_exp[4]
            A_F = create_automaton_F(model, x1, x2, t1, t2, :P)  
            # Cosmos estimation
            k3 = l_k3[i]
            command = `Cosmos $(absolute_path * "models/" * str_model * ".gspn") 
            $(absolute_path * "distance_F/dist_F_"  * str_model * ".lha") --njob $(nworkers()) 
            --const k_3=$(k3),x1=$x1,x2=$x2,t1=$t1,t2=$t2 
            --level $(level) --width $(width) 
            --verbose 0` 
            #run(command)
            run(pipeline(command, stderr=devnull))
            dict_values = cosmos_get_values("Result_dist_F_$(str_model).res")
            l_dist_cosmos[i] = dict_values["Estimated value"][1]
            nb_sim = dict_values["Total paths"][1]
            nb_accepted = dict_values["Accepted paths"][1]
            nb_sim = convert(Int, nb_sim)
            # BiochemNetABC estimation
            set_param!(ER, :k3, convert(Float64, k3))
            sync_ER = ER*A_F
            l_dist_pkg[i], nb_accepts_pkg = distribute_mean_value_lha(sync_ER, :d, nb_sim; with_accepts = true)
            #@info "About accepts" nb_sim nb_accepted nb_accepts_pkg
            test = isapprox(l_dist_cosmos[i], l_dist_pkg[i]; atol = width*1.01) || 
            (l_dist_cosmos[i] == 9997999 && l_dist_pkg[i] == Inf)
            test2 = nb_accepts_pkg == (nb_sim / nb_accepted)
            global test_all = test_all && test && test2
            if !test
                @show l_dist_pkg[i], l_dist_cosmos[i]
                @show exp
                @show k3
                @show nb_sim
            end
        end
    end
    @info exp l_dist_pkg l_dist_cosmos
end

rm("Result_dist_F_$(str_model).res")
rm("Result.res")

return test_all

