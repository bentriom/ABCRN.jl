
using Distributed
using ABCRN
addprocs(2)
@everywhere using ABCRN

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
# Check the creation of model has defined the necessary functions on all workers
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

# Synchronized model
load_model("ER")
observe_all!(ER)
load_automaton("automaton_G")
x1, x2, t1, t2 = 50.0, 100.0, 0.0, 0.8
A_G_R5 = create_automaton_G(ER, x1, x2, t1, t2, :E) 
sync_ER = ER * A_G_R5
@everywhere simulate($sync_ER)
test_all = test_all && !fetch(@spawnat 2 isdefined(Main, :sync_ER)) && !fetch(@spawnat 2 isdefined(Main, :sync_ER))
test_all = test_all && !fetch(@spawnat 3 isdefined(Main, :sync_ER)) && !fetch(@spawnat 3 isdefined(Main, :sync_ER))

rmprocs(2)

return test_all

