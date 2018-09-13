#Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[]) # JuliaLang/julia/pull/28625
ENV["GKSwstype"] = "100"

using Documenter, CopernicusUtils

#Before running this, make sure that Plots, Tensors, Distances,JLD2, Printf, Random, OrdinaryDiffEq and Clustering packages are
#installed and added to your current environment (]add )

makedocs(
    format=:html,
    sitename="CopernicusUtils.jl",
    pages = Any[
        "Home" => "index.md"
        "Example Usage" => "example.md"
    ]
    )

deploydocs(
    repo = "github.com/CoherentStructures/CopernicusUtils.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
