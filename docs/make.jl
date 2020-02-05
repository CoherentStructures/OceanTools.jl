if Base.HOME_PROJECT[] != nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

ENV["GKSwstype"] = "100"

using Documenter, CopernicusUtils

makedocs(
    sitename="CopernicusUtils.jl",
    pages = Any[
        "Home" => "index.md"
        "Example usage" => "example.md"
    ]
    )

deploydocs(
    repo = "github.com/CoherentStructures/CopernicusUtils.jl.git"
)
