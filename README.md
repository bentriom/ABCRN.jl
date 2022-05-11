
MarkovProcesses.jl 
==================

A Julia package for efficient simulation, statistical inference and verification of Continuous Time Markov Chains.

It implements:

* A core of simulation for Markov Processes.
* A simple interface for Biochemical Networks / Stochastic Petri Nets.
* Synchronized simulation with Linear Hybrid Automata.
* Approximate Bayesian Computation (a likelihood-free inference method)
* Automaton-ABC: a statistical method for verification of parametric CTMC (cite paper)

## Install

This package is not yet accessible via the Julia package manager. For the install of the package:

1. Clone this git repository on your computer.
2. Add the "src" directory of this repository to your `LOAD_PATH`. This can be done by two different ways:
    * Add the Julia line code 
    ```julia
    import Distributed: @everywhere
    @everywhere push!(LOAD_PATH, /path/to/markovprocesses.jl/src")
    ```
    on your Julia startup file which is often located in `~/.julia/config/startup.jl` in Unix systems. 
    "/path/to/markovprocesses.jl/core" is the path to the core directory of this git repository.
    
    * If you don't want to add this in your startup file, you can add these lines in your Julia script before `using MarkovProcesses`.

## Getting started 

A few notebooks are available in examples/notebooks for a quick presentation of the different features of the package.

## Tests

Execution tests and statistical tests are available. It can be run by:

`julia test/run_all.jl`

> :warning: The statistical tests run by `test/run_cosmos.jl` needs [Cosmos](http://cosmos.lacl.fr/) in your PATH environment variable.

## Benchmarks

Benchmarks have been made to test the performance of the package compared to well-known efficient other packages such as `DifferentialEquations.jl`.

## Info

This package was written during my PhD thesis. The mathematical fundations and the package archtecture are presented in it.

