
using BenchmarkTools
include("test_function_access_module.jl")
using .MyModule

load_file("test_e1.jl")
example_func2(t::Float64, x::Vector{Float64}) = t+x[2]
e2 = Edge(2, example_func2)

t, x = 1.0, 2ones(2)

@btime test(e1, t, x)
@btime test(e2, t, x)

