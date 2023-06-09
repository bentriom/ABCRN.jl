
BiochemNetABC.jl 
==================

<!--
![Pipeline status](https://gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl/badges/master/pipeline.svg)
-->
[![Documentation](https://github.com/bentriom/BiochemNetABC.jl/actions/workflows/doc.yml/badge.svg)](https://bentriom.github.io/BiochemNetABC.jl/)

A Julia package for efficient simulation, statistical inference and verification of 
Markov Processes/Continuous Time Markov Chains (CTMC) modeled by Biochemical networks/Chemical Reaction Networks (CRN) 
with Approximate Bayesian Computation methods (ABC).

It implements:

* A core for simulation of Markov Processes/CTMC,
* A simple interface for Biochemical Networks/Stochastic Petri Nets,
* Synchronized simulation with Linear Hybrid Automata.
* Approximate Bayesian Computation (a likelihood-free inference method),
* Automaton-ABC: a statistical method for verification of parametric CTMC.

<!--
## Resources

* **Documentation**: <https://2017bentrioum.pages.centralesupelec.fr/markovprocesses.jl/public>
-->

## Install

This package is not yet registered in the Julia's General registry. For the install of the package:

1. Launch Julia's REPL (for example by entering `julia` in the shell)
2. Enter Pkg's REPL by typing `]`
3. Enter
   ```julia
   pkg> add https://github.com/bentriom/BiochemNetABC.jl
   ```

## Getting started 

A few notebooks are available in examples/notebooks for a quick presentation of the packages' features.

## Tests

Execution and statistical tests can be run through:

`julia test/runtests.jl`

or in Pkg's REPL:

```julia
pkg> test BiochemNetABC
```

> :warning: The statistical tests run by `test/run_cosmos.jl` needs [Cosmos](http://cosmos.lacl.fr/) in your PATH environment variable.

## Benchmarks

Benchmarks have been made to test the package performance compared to well-known efficient other packages such as `DifferentialEquations.jl`.

## Info

The mathematical fundations and the package architecture are presented in my [PhD thesis](https://theses.hal.science/tel-03621447).

