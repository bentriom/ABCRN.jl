
# Getting Started

## Installation

1. Launch Julia's REPL (for example by entering julia in your shell)
2. Enter Pkg's REPL by typing `]`
3. Enter
   ```julia 
   pkg> add https://gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl/
   ```

## Context - Mathematical framework

In this package, we are focused on Continuous-Time Markov Chains (CTMC, also called Markov jump processes), that can be described by Chemical Reaction Networks. The future state only depends on the current state. It is defined by two properties:

- ``\forall t, s \in \mathbb{R}_{\geq 0}, \mathbb{P}(S_t | S_v, v \in [0,s]) = \mathbb{P}(S_t| S_s)`` (Memoryless/Markov property)
- ``\forall t,v \in \mathbb{R}_{\geq 0}, t > v, \mathbb{P}(S_t | S_v) = \mathbb{P}(S_{t-v} | S_0)`` (Time-homogeneity).

A Chemical Reaction Network (CRN) is a formalism that describes biological phenomena. An example is the Susceptible-Infected-Removed model (SIR):

``
R_1: S + I \xrightarrow{k_i} 2I
``
``
R_2: I \xrightarrow{k_r} R
``

This CRN has two reactions that models two phenomena: infection ($R_1$) or immunisation ($R_2$). Each reaction is parametrized by a kinetic rate ($k_i$ or $k_r$). The stochastic dynamics of a CRN are described by CTMCs.

## Models

In the package, models are objects and their types all derived from the abstract type `Model`. Let's load the SIR model, which is pre-written within the package.

```julia
julia> load_model("SIR")
create_SIR (generic function with 1 method)

julia> println(SIR)
SIRModel <: ContinuousTimeModel model
- variables :
* I (index = 2 in state space)
* R (index = 3 in state space)
* S (index = 1 in state space)
- parameters :
* ki (index = 1 in parameter space)
* kr (index = 2 in parameter space)
- transitions : R1,R2
- observed variables :
* I (index = 1 in observed state space, index = 2 in state space)
p = [0.0012, 0.05]
x0 = [95, 5, 0]
t0 = 0.0
time bound = Inf
```

`load_model` has created a variable called SIR, which contains all the information for simulating the SIR model described above. It also contains a parameter vector and an initial point.

```julia
julia> @show SIR.p 
       @show SIR.x0
SIR.p = [0.0012, 0.05]
SIR.x0 = [95, 5, 0]
```

You can change the parameters or the initial state with the functions `set_param!` and `set_x0!`.

```julia
julia> set_param!(SIR, :ki, 0.015)
julia> @show SIR.p
SIR.p = [0.015, 0.05]

julia> set_param!(SIR, [0.02, 0.07])
julia> @show SIR.p
SIR.p = [0.02, 0.07]

julia> set_x0!(SIR, :S, 93)
julia> @show SIR.x0
SIR.x0 = [93, 5, 0]
```

## Trajectories

The simulation of the model is done by the function simulate.

```julia
julia> σ = simulate(SIR)
```

`simulate` returns a trajectory, which type is derived from `AbstractTrajectory`. It can be either an object of type `Trajectory` for models `::ContinuousTimeModel` or `SynchronizedTrajectory` for models that includes an automaton (but this is the subject of another section). It is easy to access the values of a trajectory.

```julia
julia> @show σ[3] # the third state of the trajectory
       @show length_states(σ) # number of states
       @show σ[:I, 4] # Fourth value of the variable I
       @show σ.I[4] # Fourth value of the variable I
       @show get_state_from_time(σ, 2.3)

σ[3] = [7]
length_states(σ) = 196
σ[:I, 4] = 8
σ.I[4] = 8
get_state_from_time(σ, 2.3) = [79]
```

The SIR object includes an observation model symbolized by the vector `SIR.g`. Even if the variables of the model are `[:S, :I, :R]`, only I will be observed.

```julia
julia> @show SIR.map_var_idx
       @show SIR.g
       @show size(σ.values), length(σ[:I]) # Only one column which corresponds to the I variable
SIR.map_var_idx = Dict(:I => 2,:R => 3,:S => 1)
SIR.g = [:I]
(size(σ.values), length(σ[:I])) = ((1,), 196)
```

The SIR model is by default unbounded, i.e. each trajectory is simulated until it reaches an absorbing state.

```julia
julia> @show isbounded(SIR)
       @show isbounded(σ)
isbounded(SIR) = false
isbounded(σ) = false
```

We can bound the SIR's trajectories until time 120 by running `set_time_bound!(SIR, 120.0)`.

