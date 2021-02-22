
using Distributed

function f(a::Int, b::Int)
    name = Symbol("Edge_$(a)_$(b)_$(rand(1:10))")
    expr = quote
        struct $(name) tr::Int end
        check_constraints(edge::$(name), u::Int) = u >= 1.0
    end
    @everywhere eval($expr)
    eval(quote edge = $(name)(2.0) end)
end


