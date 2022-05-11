
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
    * Launch Julia's REPL (for example by entering `julia` in your terminal)
    * Enter Pkg's REPL by typing `]`
    * Type
    ```julia
    add https://gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl/
    ```

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

