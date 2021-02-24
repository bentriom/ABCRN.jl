
using Profile
using StaticArrays
using BenchmarkTools

println("Cost of 20 boolean tests:")
f_test(a::Int) = a == 2
@btime for i = 1:20
    f_test(3)
end

# Multiple dispatch
abstract type SuperType end
for i = 1:20
    include_string(Main, "struct Type$i <: SuperType v::Symbol end")
    include_string(Main, "f(item::Type$i, a::Int) = a == $(i)")
end
vec_types = SuperType[]
for i = 1:20
    rand_number_type = rand(1:20)
    @eval item = $(Symbol("Type$(rand_number_type)"))(:test)
    push!(vec_types, item)
end
println("Multiple dispatch:")
function read(vec_types::Vector{SuperType}, a::Int) 
    for item in vec_types
        f(item, a)
    end
end
@btime read(vec_types, 3)

# Parametric struct
struct AnotherSuperType{T<:Function}
    v::Symbol 
    f::T
end
vec_others_types = AnotherSuperType[]
for i = 1:20
    rand_number_type = rand(1:20)
    include_string(Main, "f_other_$(i)(a::Int) = a == $(i)")
    sym_func = Symbol("f_other_$i")
    @eval push!(vec_others_types, AnotherSuperType(:test, $(sym_func)))
end
println("Parametrized struct:")
function read(vec_others_types::Vector{AnotherSuperType}, a::Int)
    for item in vec_others_types
        item.f(a)
    end
end
@btime read(vec_others_types, 3)

# With two vectors
println("Vectors:")
vec_tr = Union{Nothing,Symbol}[]
vec_func = Function[]
for i = 1:20
    rand_number_type = rand(1:20)
    include_string(Main, "f_vec_$(i)(a::Int) = a == $(i)")
    sym_func = Symbol("f_vec_$i")
    push!(vec_tr, :test)
    @eval push!(vec_func, $(sym_func))
end
function read(vec_tr::Vector{Union{Nothing,Symbol}}, vec_func::Vector{Function}, a::Int) 
    for i = eachindex(vec_tr)
        vec_func[i](a)
    end
end
@btime read(vec_tr, vec_func, 3)

exit()
println("Transitions second:")
vec_tr = Union{Nothing,Symbol}[]
vec_func = Function[]
for i = 1:20
    rand_number_type = rand(1:20)
    name_func = Symbol("f_vec2_$(i)")
    str_func = "$(name_func)(a::Int) = a == $(rand_number_type)"
    include_string(Main, str_func)
    push!(vec_tr, :test)
    @eval push!(vec_func, $(name_func))
end
@btime read(vec_tr, vec_func, 3)


# How to store and read efficiently abstract types ?
d = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()
for i = 1:20
    sym_a = Symbol("a$i")
    sym_b = Symbol("a$i")
    d[sym_a] = Dict{Symbol,Vector{Float64}}()
    d[sym_a][sym_b] = zeros(4)
end
function f(dict::Dict)
    for key in keys(dict)
        for key2 in keys(dict[key])
            for i = eachindex(dict[key][key2])
                dict[key][key2][i]
            end
        end
    end
end
@btime f(d)

