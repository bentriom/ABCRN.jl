
using BenchmarkTools
using MarkovProcesses

load_model("ER")

b1 = @benchmark begin 
    Threads.@threads for i = 1:1000
        simulate(ER)
    end
end

b2 = @benchmark begin 
    for i = 1:1000
        Threads.@spawn simulate(ER)
    end
end


b3 = @benchmark begin 
    for i = 1:1000
        simulate(ER)
    end
end

@show minimum(b1), mean(b1), maximum(b1)
@show minimum(b2), mean(b2), maximum(b2)
@show minimum(b3), mean(b3), maximum(b3)

