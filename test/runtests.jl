
using Distributed
using Test
using ABCRN
using Logging
Logging.disable_logging(Logging.Warn)

include("run_unit.jl")
include("run_simulation.jl")
include("run_dist_lp.jl")
include("run_automata.jl")
include("run_abc_smc.jl")
include("run_abc_model_choice.jl")

