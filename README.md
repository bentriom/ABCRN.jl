
ABCRN.jl 
==================

<!--
![Pipeline status](https://gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl/badges/master/pipeline.svg)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://2017bentrioum.pages.centralesupelec.fr/markovprocesses.jl/)
-->

A Julia package for efficient simulation, statistical inference and verification of Continuous Time Markov Chains
modeled by Chemical Reaction Networks with Approximate Bayesian Computation methods.

It implements:

* A core for simulation of Markov Processes.
* A simple interface for Biochemical Networks / Stochastic Petri Nets.
* Synchronized simulation with Linear Hybrid Automata.
* Approximate Bayesian Computation (a likelihood-free inference method)
* Automaton-ABC: a statistical method for verification of parametric CTMC

<!--
## Resources

* **Documentation**: <https://2017bentrioum.pages.centralesupelec.fr/markovprocesses.jl/public>
-->

## Install

This package is not yet registered in the Julia's packages. For the install of the package:

1. Launch Julia's REPL (for example by entering `julia` in the shell)
2. Enter Pkg's REPL by typing `]`
3. Enter
   ```julia
   pkg> add https://github.com/bentriom/ABCRN.jl
   ```

## Getting started 

A few notebooks are available in examples/notebooks for a quick presentation of the different features of the package.

## Tests

Execution and statistical tests are available. It can be run by

`julia test/runtests.jl`

or in Pkg's REPL:

```julia
pkg> test ABCRN
```

> :warning: The statistical tests run by `test/run_cosmos.jl` needs [Cosmos](http://cosmos.lacl.fr/) in your PATH environment variable.

## Benchmarks

Benchmarks have been made to test the performance of the package compared to well-known efficient other packages such as `DifferentialEquations.jl`.

## Info

This package was written during my PhD thesis. The mathematical fundations and the package archtecture are presented in it.

