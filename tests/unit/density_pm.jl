
using MarkovProcesses
import Distributions: rand, pdf, product_distribution
load_model("SIR")

tol = 0.0

dist = Uniform(0.0, 1.0)
pm = ParametricModel(SIR, ("kr", dist)) 
test1 = !insupport(pm, [2.0])

dist = Normal()
pm = ParametricModel(SIR, ("kr", dist)) 
test2 = isapprox(prior_pdf(pm, [0.05]), pdf(dist, 0.05); atol = tol)


mat_u = [[rand()] for i = 1:10]
vec_u = [mat_u[i][1] for i = 1:10]
_prior_pdf(x::Vector{Float64}) = prior_pdf(pm, x)
test3 = isapprox(_prior_pdf.(mat_u), pdf.(dist, vec_u); atol = tol)

dist, dist2 = Normal(), Normal(1.0, 2.0)
prod_dist = product_distribution([dist, dist2])
pm = ParametricModel(SIR, ("ki", dist), ("kr", dist2)) 
mat_u = rand(2,10)
vec_u = [mat_u[:,i] for i = 1:10]
vec_res = zeros(10)
prior_pdf!(vec_res, pm, mat_u)
test4 = isapprox(vec_res, [pdf(prod_dist, vec_u[i]) for i = 1:10]; atol = tol)

return test1 && test2 && test3 && test4

