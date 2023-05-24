using Documenter
using TelecomUtils

makedocs(
    sitename = "TelecomUtils",
    format = Documenter.HTML(),
    modules = [TelecomUtils]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/disberd/TelecomUtils.jl.git"
)
