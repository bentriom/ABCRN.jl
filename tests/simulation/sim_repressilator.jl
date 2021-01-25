
using MarkovProcesses

load_model("repressilator")
set_param!(repressilator, [:α, :n, :β, :α0], [1.0, 2.0, 5.0, 1000.0])
set_x0!(repressilator, [:m1, :m2, :m3, :p1, :p2, :p3], [0, 2, 0, 1, 2, 3])

return true

