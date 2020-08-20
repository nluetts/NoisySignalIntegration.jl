using Pkg

Pkg.activate(@__DIR__) # activate the docs environment

using Documenter
using MCIntegrate

# include doctests from docstrings
DocMeta.setdocmeta!(MCIntegrate, :DocTestSetup, :(using MCIntegrate); recursive=true)

makedocs(
    sitename="NoisySignalIntegration.jl Documentation",
    modules=[MCIntegrate],
    pages = [
        "Home" => "index.md",
        "Usage Guide" => "guide.md",
        "API reference" => "API.md"
    ],
    highlightsig=true
)