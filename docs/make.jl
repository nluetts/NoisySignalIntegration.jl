using Pkg

Pkg.activate(@__DIR__) # activate the docs environment

using Documenter
using MCIntegrate

# include doctests from docstrings
DocMeta.setdocmeta!(MCIntegrate, :DocTestSetup, :(using MCIntegrate); recursive=true)

makedocs(
    sitename="NoisySignalIntegration.jl",
    modules=[MCIntegrate],
    pages = [
        "Home" => "index.md",
        "Package Overview" => "overview.md",
        "Usage Guide" => "guide.md",
        "Examples" => "examples.md",
        "Internals" => "internals.md",
        "API reference" => "API.md"
    ],
    highlightsig=true
)