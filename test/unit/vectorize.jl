
using MarkovProcesses
load_model("ER")

test_all = true

# Unbounded
tml = 0:0.1:10.0
σ = simulate(ER)
vec_traj = vectorize(σ, :P, tml)
idx = 0
for i = eachindex(vec_traj)
    local test = convert(Int, vec_traj[i]) == get_var_from_time(σ, :P, tml[i])
    global test_all = test_all && test
    if !test
        @show tml[i], convert(Int, vec_traj[i]), get_var_from_time(σ, :P, tml[i])
    end
end

# Bounded
set_time_bound!(ER, 2.0)
tml = 0:0.1:2.0
σ = simulate(ER)
vec_traj = vectorize(σ, :P, tml)
for i = eachindex(vec_traj)
    local test = convert(Int, vec_traj[i]) == get_var_from_time(σ, :P, tml[i])
    global test_all = test_all && test
    if !test
        @show tml[i], convert(Int, vec_traj[i]), get_var_from_time(σ, :P, tml[i])
    end
end

set_time_bound!(ER, 1.0)
tml = 0:0.02:0.9
σ = simulate(ER)
vec_traj = vectorize(σ, :P, tml)
for i = eachindex(vec_traj)
    local test = convert(Int, vec_traj[i]) == get_var_from_time(σ, :P, tml[i])
    global test_all = test_all && test
    if !test
        @show tml[i], convert(Int, vec_traj[i]), get_var_from_time(σ, :P, tml[i])
    end
end

set_time_bound!(ER, 2.0)
tml = 0.5:0.1:2.0
σ = simulate(ER)
vec_traj = vectorize(σ, :P, tml)
for i = eachindex(vec_traj)
    local test = convert(Int, vec_traj[i]) == get_var_from_time(σ, :P, tml[i])
    global test_all = test_all && test
    if !test
        @show tml[i], convert(Int, vec_traj[i]), get_var_from_time(σ, :P, tml[i])
    end
end

return test_all

