if Base.HOME_PROJECT[] != nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

ENV["GKSwstype"] = "100"

using Documenter, OceanTools

makedocs(
    sitename="OceanTools.jl",
    pages = Any[
        "Home" => "index.md"
        "Loading NetCDF datafiles" => "loading_data.md"
        "Interpolation" => "interpolation.md"
        "Example" => "example.md"
    ]
    )

deploydocs(
    repo = "github.com/CoherentStructures/OceanTools.jl.git"
)
