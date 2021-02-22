
import Statistics: mean
using BenchmarkTools
using Profile

EdgeTransition = Union{Nothing,Vector{Symbol}}
struct EdgeStruct
    tr::EdgeTransition
    func1::Function
    func2::Function
end
struct EdgeStruct2
    tr::EdgeTransition
    func1::Symbol
    func2::Symbol
end
struct EdgeStruct3{F1 <: Function, F2 <: Function}
    tr::EdgeTransition
    func1::F1
    func2::F2
end
EdgeTuple = Tuple{EdgeTransition,Function,Function}
EdgeTuple2 = Tuple{EdgeTransition,Symbol,Symbol}

function f(t::Float64, values::Vector{Float64}, x::Vector{Int}, p::Vector{Float64})
    return (t <= 0.025 && (values[1] < 50.0 || values[1] > 75.0))
end

t = 0.1
values = [100.0, Inf, 0.0]
x = [99, 99, 1, 0] 
p = [1.0, 1.0]

edge_struct_1 = EdgeStruct(nothing, f, mean)
edge_struct_2 = EdgeStruct2(nothing, :f, :mean)
edge_struct_3 = EdgeStruct3(nothing, f, mean)
edge_tuple_1 = (nothing, getfield(Main, :f), getfield(Main, :mean))
edge_tuple_2 = (nothing, :f, :mean)
@assert typeof(edge_struct_1) <: EdgeStruct && typeof(edge_struct_2) <: EdgeStruct2 &&
        typeof(edge_tuple_1) <: EdgeTuple && typeof(edge_tuple_2) <: EdgeTuple2

println("Time execution of f")
@btime f(t, values, x, p)

println("Time execution of f with edges")
println("- Structs")
@btime edge_struct_1.func1(t, values, x, p)
@btime getfield(Main, edge_struct_2.func1)(t, values, x, p)
@btime edge_struct_3.func1(t, values, x, p)

println("- Tuples")
@btime edge_tuple_1[2](t, values, x, p)
@btime getfield(Main, edge_tuple_2[2])(t, values, x, p)

println("Time access of variables")
println("- Structs")
@btime edge_struct_1.func1
@btime edge_struct_2.func1
@btime edge_struct_3.func1
println("- Tuples")
@btime edge_tuple_1[2]
@btime edge_tuple_2[2]

