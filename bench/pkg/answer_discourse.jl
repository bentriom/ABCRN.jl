
import Statistics: mean
using BenchmarkTools

EdgeTransition = Union{Nothing,Vector{Symbol}}
struct EdgeStruct3{F1 <: Function, F2 <: Function}
    tr::EdgeTransition
    func1::F1
    func2::F2
end

function f(t::Float64, values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64})
    return (t <= 0.025 && (values[1] < 50.0 || values[1] > 75.0))
end

t = 0.1
values = [100.0, Inf, 0.0]
x = [99, 99, 1, 0] 
p = [1.0, 1.0]

edge_struct_3 = EdgeStruct3(nothing, f, mean)
@btime edge_struct_3.func1(t, values, x, p)

