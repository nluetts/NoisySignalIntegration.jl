using Documenter
using MCIntegrate

makedocs(
    sitename="MCIntegrate.jl Documentation",
    pages = [
        "Home" => "index.md",
        "Usage Guide" => "guide.md"
    ],
    highlightsig=true
)