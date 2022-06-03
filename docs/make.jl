using Documenter
using MDTools

makedocs(
    sitename = "MDTools",
    format = Documenter.HTML(),
    modules = [MDTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
