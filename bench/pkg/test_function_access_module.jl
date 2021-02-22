
module MyModule

export Edge, load_file, test

struct Edge
    a::Int
    func::Function
end

function test(e::Edge, t::Float64, x::Vector{Float64})
    for i = 1:1E6
        e.func(t, x)
    end
end

load_file(str::String) = include(str)

end

