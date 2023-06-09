
using Documenter, BiochemNetABC

makedocs(
    sitename = "BiochemNetABC.jl",
    modules = [BiochemNetABC],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "starting.md",
        "Create a model" => "create_model.md",
        "API" => Any[
            "Model" => "api/model.md",
            "Trajectory" => "api/trajectory.md",
            "Approximate Bayesian Computation" => "api/abc.md",
            "Plots" => "api/plots.md"
        ]
    ],
    format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/bentriom/BiochemNetABC.jl.git",
    devbranch = "master",
)

