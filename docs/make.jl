push!(LOAD_PATH,"../src/")
using Dispersions
using Documenter

makedocs(
    modules = [Dispersions],
    sitename = "Dispersions.jl",
    authors="Julian Stobbe <Atomtomate@gmx.de> and contributors",
    repo="https://github.com/Atomtomate/SeriesAcceleration.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://Atomtomate.github.io/Dispersions.jl",
        assets=String[],
    ),
    pages = ["Home" => "index.md"
             "Interface" => "interface.md"
            ]
)

deploydocs(
    branch="gh-pages",
    devbranch = "master",
    devurl = "stable",
    repo = "github.com/Atomtomate/Dispersions.jl.git"
)
