using Documenter
using MCIntegrate
using Pkg

Pkg.activate(".") # activate the docs environment

makedocs(
    sitename="MCIntegrate.jl Documentation",
    pages = [
        "Home" => "index.md",
        "Usage Guide" => "guide.md"
    ],
    highlightsig=true
)