
example_func(t::Float64, x::Vector{Float64}) = t+x[1]
e1 = Edge(1, example_func)

export example_func, e1

