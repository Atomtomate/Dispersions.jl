using Documenter
using Dispersions

makedocs(
    sitename = "Dispersions",
    format = Documenter.HTML(),
    modules = [Dispersions]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
