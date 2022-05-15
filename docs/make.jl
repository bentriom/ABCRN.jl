
using Documenter, MarkovProcesses

makedocs(
    sitename = "MarkovProcesses.jl",
    modules = [MarkovProcesses],
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
config_gitlab = Documenter.GitLab()

deploydocs(
    repo = "gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl.git",
    devbranch = "ci_tests",
    branch = "pages"
)

