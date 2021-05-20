
using Documenter, MarkovProcesses

makedocs(
    sitename = "MarkovProcesses.jl",
    modules = [MarkovProcesses],
    pages = [
        "Home" => "index.md",
        "Approximate Bayesian Computation" => "abc.md"
    ],
    format = Documenter.HTML(prettyurls = false)
)

