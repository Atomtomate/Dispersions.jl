using Documenter
push!(LOAD_PATH,"../src/")
include("../src/Dispersions.jl")
using .Dispersions

makedocs(
    sitename = "Dispersions.jl",
    format = Documenter.HTML(),
    modules = [Dispersions],
    pages = ["index.md"
             "Interface" => "interface.md"
            ]
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Atomtomate/Dispersions.jl.git",
    devbranch = "master",
    devurl = "stable"
)
