
using Documenter, ABCRN

makedocs(
    sitename = "ABCRN.jl",
    modules = [ABCRN],
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
    repo = "github.com/bentriom/ABCRN.jl.git",
    devbranch = "master",
)

