
using Distributed
using MarkovProcesses
# A bunch of code to define the package on the workers
# because they are created after the execution of Julia
addprocs(2)
module_path = get_module_path()
@everywhere module_path = $module_path
@everywhere push!(LOAD_PATH, "$(module_path)/core")
@everywhere using MarkovProcesses

test_all = true

# SIR created by macro
model_SIR = @network_model begin
    R1: (S+I => 2I, ki*S*I)
    R2: (I => R, kr*I)
end "SIR"
set_x0!(model_SIR, [95,5,0])
set_param!(model_SIR, [0.012, 0.05])
# Check if we send the model to the workers, 
# then the simulation works
@everywhere simulate($model_SIR)
sym_f! = Symbol(model_SIR.f!)
sym_isabs = Symbol(model_SIR.isabsorbing)
test_all = test_all && fetch(@spawnat 2 isdefined(Main, sym_f!)) && fetch(@spawnat 2 isdefined(Main, sym_isabs))
test_all = test_all && fetch(@spawnat 3 isdefined(Main, sym_f!)) && fetch(@spawnat 3 isdefined(Main, sym_isabs))
# Check the model object is not defined on workers
test_all = test_all && !fetch(@spawnat 2 isdefined(Main, :model_SIR)) && !fetch(@spawnat 2 isdefined(Main, :model_SIR))
test_all = test_all && !fetch(@spawnat 3 isdefined(Main, :model_SIR)) && !fetch(@spawnat 3 isdefined(Main, :model_SIR))

# SIR pkg
load_model("SIR")
# Check if we send the model to the workers, 
# then the simulation works
@everywhere simulate($SIR)
sym_f! = Symbol(SIR.f!)
sym_isabs = Symbol(SIR.isabsorbing)
test_all = test_all && fetch(@spawnat 2 isdefined(Main, sym_f!)) && fetch(@spawnat 2 isdefined(Main, sym_isabs))
test_all = test_all && fetch(@spawnat 3 isdefined(Main, sym_f!)) && fetch(@spawnat 3 isdefined(Main, sym_isabs))
# Check the model object is not defined on workers
test_all = test_all && !fetch(@spawnat 2 isdefined(Main, :SIR)) && !fetch(@spawnat 2 isdefined(Main, :SIR))
test_all = test_all && !fetch(@spawnat 3 isdefined(Main, :SIR)) && !fetch(@spawnat 3 isdefined(Main, :SIR))

# ER pkg
load_model("ER")
# Check if we send the model to the workers, 
# then the simulation works
@everywhere simulate($ER)
sym_f! = Symbol(ER.f!)
sym_isabs = Symbol(ER.isabsorbing)
test_all = test_all && fetch(@spawnat 2 isdefined(Main, sym_f!)) && fetch(@spawnat 2 isdefined(Main, sym_isabs))
test_all = test_all && fetch(@spawnat 3 isdefined(Main, sym_f!)) && fetch(@spawnat 3 isdefined(Main, sym_isabs))
# Check the model object is not defined on workers
test_all = test_all && !fetch(@spawnat 2 isdefined(Main, :ER)) && !fetch(@spawnat 2 isdefined(Main, :ER))
test_all = test_all && !fetch(@spawnat 3 isdefined(Main, :ER)) && !fetch(@spawnat 3 isdefined(Main, :ER))

rmprocs(2)

return test_all

