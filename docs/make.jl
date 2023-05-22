
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
config_gitlab = Documenter.GitLab()

deploydocs(
    repo = "gitlab-research.centralesupelec.fr/2017bentrioum/markovprocesses.jl.git",
    deploy_config = config_gitlab,
    devbranch = "ci_tests",
    devurl = "dev",
    branch = "pages"
)

